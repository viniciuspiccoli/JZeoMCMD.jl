#!/usr/bin/env julia
"""
reload_adsorbate.jl  —  v1.1

Replace adsorbate molecules in an existing LAMMPS data file with
new positions from a RASPA3 JSON restart. Framework + all coefficients
are preserved exactly. Writes complete class2 coefficients.

Usage:
  julia reload_adsorbate.jl <prev_npt.lmp> <raspa_restart.json> <output.lmp>

Dependencies (same directory):
  read_lammps_data.jl, build_loaded_zeolite.jl
"""

# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════

# for silicous zeolites
Base.@kwdef struct ReloadConfig
    adsorbate_types::Vector{Int} = [3, 4, 5, 6]
    framework_types::Vector{Int} = [1, 2]

    si_type::Int = 1
    o_type::Int = 2

    eth_atoms_per_mol::Int = 4
    eth_atom_names::Vector{String} = ["CH3", "CH2", "O_eth", "H_eth"]
    eth_types::Dict{String,Int} = Dict(
        "CH3"=>3, "CH2"=>4, "O_eth"=>5, "H_eth"=>6)
    eth_charges::Dict{String,Float64} = Dict(
        "CH3"=>0.0, "CH2"=>0.265, "O_eth"=>-0.700, "H_eth"=>0.435)
    eth_masses::Dict{String,Float64} = Dict(
        "CH3"=>15.035, "CH2"=>14.027, "O_eth"=>15.999, "H_eth"=>1.008)

    eth_bond_defs::Vector{Tuple{Int,Int,Int}} = [(2,1,2),(3,2,3),(4,3,4)]
    eth_angle_defs::Vector{Tuple{Int,Int,Int,Int}} = [(3,1,2,3),(4,2,3,4)]
    eth_dihedral_defs::Vector{Tuple{Int,Int,Int,Int,Int}} = [(2,1,2,3,4)]

    box_tolerance::Float64 = 0.5
    check_xyz::String = "reloaded_ethanol_check.xyz"
end



#=
for aluminosilicate

Base.@kwdef struct ReloadConfig
    adsorbate_types::Vector{Int} = [3, 4, 5, 6]
    framework_types::Vector{Int} = [1, 2]

    si_type::Int = 1
    o_type::Int = 2

    eth_atoms_per_mol::Int = 4
    eth_atom_names::Vector{String} = ["CH3", "CH2", "O_eth", "H_eth"]
    eth_types::Dict{String,Int} = Dict(
        "CH3"=>3, "CH2"=>4, "O_eth"=>5, "H_eth"=>6)
    eth_charges::Dict{String,Float64} = Dict(
        "CH3"=>0.0, "CH2"=>0.265, "O_eth"=>-0.700, "H_eth"=>0.435)
    eth_masses::Dict{String,Float64} = Dict(
        "CH3"=>15.035, "CH2"=>14.027, "O_eth"=>15.999, "H_eth"=>1.008)

    eth_bond_defs::Vector{Tuple{Int,Int,Int}} = [(2,1,2),(3,2,3),(4,3,4)]
    eth_angle_defs::Vector{Tuple{Int,Int,Int,Int}} = [(3,1,2,3),(4,2,3,4)]
    eth_dihedral_defs::Vector{Tuple{Int,Int,Int,Int,Int}} = [(2,1,2,3,4)]

    box_tolerance::Float64 = 0.5
    check_xyz::String = "reloaded_ethanol_check.xyz"
end

=#


# ═══════════════════════════════════════════════════════════════════
# Step 1: Strip adsorbate from existing data
# ═══════════════════════════════════════════════════════════════════

function strip_adsorbate(data, cfg::ReloadConfig)
    ads_set = Set(cfg.adsorbate_types)
    natoms = size(data.coords, 1)

    fw_indices = [j for j in 1:natoms if !(data.atom_labels[j] in ads_set)]
    nfw = length(fw_indices)
    n_removed = natoms - nfw

    # old atom ID → new atom ID
    id_map = Dict{Int,Int}()
    for (new_idx, old_idx) in enumerate(fw_indices)
        id_map[data.atom_ids[old_idx]] = new_idx
    end

    # Filter atoms
    data.coords          = data.coords[fw_indices, :]
    data.atom_ids        = collect(1:nfw)
    data.atom_labels     = data.atom_labels[fw_indices]
    data.atom_charges    = data.atom_charges[fw_indices]
    data.molecule_labels = data.molecule_labels[fw_indices]
    if !isempty(data.image_flags) && size(data.image_flags,1) == natoms
        data.image_flags = data.image_flags[fw_indices, :]
    else
        data.image_flags = zeros(Int, nfw, 3)
    end

    # Filter + remap bonds
    if size(data.bonds, 1) > 0
        keep = [all(haskey(id_map, data.bonds[k,c]) for c in 1:2)
                for k in 1:size(data.bonds,1)]
        ki = findall(keep)
        data.bond_labels = data.bond_labels[ki]
        nb = data.bonds[ki, :]
        for k in 1:size(nb,1), c in 1:2
            nb[k,c] = id_map[nb[k,c]]
        end
        data.bonds = nb
    end

    # Filter + remap angles
    if size(data.angles, 1) > 0
        keep = [all(haskey(id_map, data.angles[k,c]) for c in 1:3)
                for k in 1:size(data.angles,1)]
        ki = findall(keep)
        data.angle_labels = data.angle_labels[ki]
        na = data.angles[ki, :]
        for k in 1:size(na,1), c in 1:3
            na[k,c] = id_map[na[k,c]]
        end
        data.angles = na
    end

    # Filter + remap dihedrals
    if size(data.dihedrals, 1) > 0
        keep = [all(haskey(id_map, data.dihedrals[k,c]) for c in 1:4)
                for k in 1:size(data.dihedrals,1)]
        ki = findall(keep)
        data.dihedral_labels = data.dihedral_labels[ki]
        nd = data.dihedrals[ki, :]
        for k in 1:size(nd,1), c in 1:4
            nd[k,c] = id_map[nd[k,c]]
        end
        data.dihedrals = nd
    end

    return data, n_removed
end

# ═══════════════════════════════════════════════════════════════════
# Step 2: Read new ethanol from RASPA3 JSON
# ═══════════════════════════════════════════════════════════════════

function read_raspa3_molecules(json_path::String, box_dims, cfg::ReloadConfig)
    restart = JSON.parsefile(json_path)

   # component_key = ""
   # for k in keys(restart)
   #     if restart[k] isa AbstractVector && !isempty(restart[k]) &&
   #        restart[k][1] isa AbstractVector
   #         component_key = k; break
   #     end
   # end
   # isempty(component_key) && error("No adsorbate found in JSON")
   # println("  Component: \"$component_key\"")  

    # NEW:
    component_key = ""
    for k in keys(restart)
        if restart[k] isa AbstractVector
            component_key = k
            break
        end
    end
    if isempty(component_key)
        println("  WARNING: No adsorbate component in JSON")
        return []
    end
    println("  Component: \"$component_key\"")
    ##########
   
    all_pos = restart[component_key]
    apm = cfg.eth_atoms_per_mol
    nmols_total = length(all_pos) ÷ apm

    L  = [box_dims[d,2]-box_dims[d,1] for d in 1:3]
    lo = [box_dims[d,1] for d in 1:3]

    molecules = Vector{Vector{Tuple{String,Float64,Float64,Float64}}}()
    n_outside = 0

    for mi in 1:nmols_total
        raw = Matrix{Float64}(undef, apm, 3)
        for k in 1:apm
            p = all_pos[(mi-1)*apm + k]
            raw[k,:] = Float64.([p[1], p[2], p[3]])
        end

        # Wrap as rigid unit (COM-based)
        com = vec(mean(raw, dims=1))
        com_w = [lo[d] + mod(com[d]-lo[d], L[d]) for d in 1:3]
        shift = com_w .- com
        wrapped = raw .+ shift'

        ok = all(sqrt(sum((wrapped[k,:] .- wrapped[k+1,:]).^2)) < 3.0
                 for k in 1:(apm-1))
        tol = cfg.box_tolerance
        inside = all(lo[d]-tol <= com_w[d] <= lo[d]+L[d]+tol for d in 1:3)

        if !inside || !ok
            n_outside += 1; continue
        end

        mol = [(cfg.eth_atom_names[k], wrapped[k,1], wrapped[k,2], wrapped[k,3])
               for k in 1:apm]
        push!(molecules, mol)
    end

    println("  RASPA3: $nmols_total → $(length(molecules)) kept")
    n_outside > 0 && println("  Removed $n_outside outside/broken")
    return molecules
end

# ═══════════════════════════════════════════════════════════════════
# Step 3: Append new ethanol
# ═══════════════════════════════════════════════════════════════════

function append_adsorbate!(data, molecules, cfg::ReloadConfig)
    nmols = length(molecules)
    nmols == 0 && (println("  No molecules."); return data)

    nfw = size(data.coords, 1)
    apm = cfg.eth_atoms_per_mol
    neth = nmols * apm
    ntot = nfw + neth
    max_mol = maximum(data.molecule_labels; init=0)

    # Extend atoms
    oc = data.coords;     data.coords = zeros(ntot,3);   data.coords[1:nfw,:] = oc
    oi = data.atom_ids;   data.atom_ids = zeros(Int,ntot); data.atom_ids[1:nfw] = oi
    ot = data.atom_labels; data.atom_labels = zeros(Int,ntot); data.atom_labels[1:nfw] = ot
    oq = data.atom_charges; data.atom_charges = zeros(ntot); data.atom_charges[1:nfw] = oq
    om = data.molecule_labels; data.molecule_labels = zeros(Int,ntot); data.molecule_labels[1:nfw] = om
    og = data.image_flags; data.image_flags = zeros(Int,ntot,3)
    size(og,1)==nfw && (data.image_flags[1:nfw,:] = og)

    for (mi, mol) in enumerate(molecules)
        for (k, (name, x, y, z)) in enumerate(mol)
            idx = nfw + (mi-1)*apm + k
            data.atom_ids[idx] = idx
            data.atom_labels[idx] = cfg.eth_types[name]
            data.atom_charges[idx] = cfg.eth_charges[name]
            data.molecule_labels[idx] = max_mol + mi
            data.coords[idx,:] = [x, y, z]
        end
    end

    # Append ethanol bonds
    nfb = size(data.bonds,1); neb = nmols*length(cfg.eth_bond_defs)
    nb = zeros(Int,nfb+neb,2); nb[1:nfb,:] = data.bonds
    nbl = zeros(Int,nfb+neb);  nbl[1:nfb] = data.bond_labels
    for mi in 1:nmols, (bi,(bt,a1,a2)) in enumerate(cfg.eth_bond_defs)
        i = nfb+(mi-1)*length(cfg.eth_bond_defs)+bi
        base = nfw+(mi-1)*apm
        nbl[i]=bt; nb[i,:]=[base+a1,base+a2]
    end
    data.bonds=nb; data.bond_labels=nbl

    # Append ethanol angles
    nfa = size(data.angles,1); nea = nmols*length(cfg.eth_angle_defs)
    na = zeros(Int,nfa+nea,3); na[1:nfa,:] = data.angles
    nal = zeros(Int,nfa+nea);  nal[1:nfa] = data.angle_labels
    for mi in 1:nmols, (ai,(at,a1,a2,a3)) in enumerate(cfg.eth_angle_defs)
        i = nfa+(mi-1)*length(cfg.eth_angle_defs)+ai
        base = nfw+(mi-1)*apm
        nal[i]=at; na[i,:]=[base+a1,base+a2,base+a3]
    end
    data.angles=na; data.angle_labels=nal

    # Append ethanol dihedrals
    nfd = size(data.dihedrals,1); ned = nmols*length(cfg.eth_dihedral_defs)
    nd = zeros(Int,nfd+ned,4); nd[1:nfd,:] = data.dihedrals
    ndl = zeros(Int,nfd+ned);  ndl[1:nfd] = data.dihedral_labels
    for mi in 1:nmols, (di,(dt,a1,a2,a3,a4)) in enumerate(cfg.eth_dihedral_defs)
        i = nfd+(mi-1)*length(cfg.eth_dihedral_defs)+di
        base = nfw+(mi-1)*apm
        ndl[i]=dt; nd[i,:]=[base+a1,base+a2,base+a3,base+a4]
    end
    data.dihedrals=nd; data.dihedral_labels=ndl

    data.nbond_types = max(data.nbond_types, 4)
    data.nangle_types = max(data.nangle_types, 4)
    data.ndihedral_types = max(data.ndihedral_types, 2)
    for (name,t) in cfg.eth_types; data.masses[t] = cfg.eth_masses[name]; end

    return data
end

# ═══════════════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════════════

"""
    reload_adsorbate(prev_data, json_restart, output)

Read previous NPT data → strip old ethanol → read new from RASPA3 →
append → write complete data file with all class2 coefficients.
"""
function reload_adsorbate(prev_data::String, json_restart::String,
                           output::String;
                           cfg::ReloadConfig = ReloadConfig())
    println("╔═══════════════════════════════════════════════╗")
    println("║  reload_adsorbate.jl v1.1                    ║")
    println("╚═══════════════════════════════════════════════╝\n")

    println("═══ Reading $prev_data ═══")
    data = read_lammps_data(prev_data; verbose=false)
    println("  Atoms: $(size(data.coords,1))")

    println("\n═══ Stripping old adsorbate ═══")
    data, n_removed = strip_adsorbate(data, cfg)
    nfw = size(data.coords, 1)
    println("  Removed $n_removed adsorbate atoms → $nfw framework atoms")

    println("\n═══ Reading new ethanol from $json_restart ═══")
    molecules = read_raspa3_molecules(json_restart, data.box_dimensions, cfg)

    println("\n═══ Appending $(length(molecules)) molecules ═══")
    append_adsorbate!(data, molecules, cfg)
    ntot = size(data.coords, 1)
    @printf("  Final: %d atoms (%d fw + %d eth), charge = %.6f e\n",
            ntot, nfw, ntot-nfw, sum(data.atom_charges))

    # Write ethanol check XYZ
    if length(molecules) > 0
        open(cfg.check_xyz, "w") do io
            neth = sum(length(m) for m in molecules)
            println(io, neth)
            println(io, "Reloaded ethanol ($(length(molecules)) mols)")
            for mol in molecules, (name,x,y,z) in mol
                elem = name=="H_eth" ? "H" : (name=="O_eth" ? "O" : "C")
                @printf(io, "%s %.6f %.6f %.6f\n", elem, x, y, z)
            end
        end
        println("  Check: $(cfg.check_xyz)")
    end

    # Write using write_complete_data from build_loaded_zeolite.jl
    # which includes all class2 coefficients
    #println("\n═══ Writing $output (with all coefficients) ═══")
    #zcfg = ZeoliteConfig(output_data=output)
    #generate_impropers!(data; si_type=cfg.si_type, o_type=cfg.o_type)
    #write_complete_data(output, data, zcfg)

    println("\n═══ Writing $output (with all coefficients) ═══")
    if length(cfg.framework_types) <= 2
        # Silica path: generate Si impropers + embed class2 coefficients
        zcfg = ZeoliteConfig(output_data=output)
        generate_impropers!(data; si_type=cfg.si_type, o_type=cfg.o_type)
        write_complete_data(output, data, zcfg)
    else
        # Aluminosilicate path: framework impropers preserved from cycle 1,
        # coefficients come from .ff include file
        write_lammps_data(output, data;
        comment="Reloaded alumino framework + ethanol (reload_adsorbate.jl)")
    end

end

# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

function main_reload(args=ARGS)
    if length(args) < 3
        println("""
        Usage: julia reload_adsorbate.jl <prev_npt.lmp> <raspa.json> <output.lmp>

        Example:
          julia reload_adsorbate.jl loaded_npt_final.lmp \\
                                    restart_300_1e+06.s0.json \\
                                    cycle2_loaded.lmp
        """)
        return
    end
    reload_adsorbate(args[1], args[2], args[3])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_reload()
end
