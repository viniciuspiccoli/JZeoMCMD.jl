#!/usr/bin/env julia
"""
build_loaded_zeolite.jl  —  v1.1

Pipeline to build a LAMMPS data file for a zeolite loaded with adsorbate
molecules, starting from:
  1. Ovito-exported framework data file (coordinates only)
  2. RASPA3 JSON restart file (adsorbate positions)
  3. Hill-Sauer force field (table + class2 bonded)

Supports: all-silica and aluminosilicate frameworks.

Usage:
  julia build_loaded_zeolite.jl <ovito.data> [raspa_restart.json] [output.lmp]

"""

#using Printf
#using Statistics

# ═══════════════════════════════════════════════════════════════════
# Load dependencies
# ═══════════════════════════════════════════════════════════════════
#include(joinpath(@__DIR__, "read_lammps_data.jl"))
#include(joinpath(@__DIR__, "add_zeolite_topology.jl"))

#local _JSON
#try
#    global _JSON = Base.require(Base.PkgId(
#        Base.UUID("682c06a0-de6a-54ab-a142-c8b1cf79cde6"), "JSON"))
#catch
#    error("JSON.jl required. Run: using Pkg; Pkg.add(\"JSON\")")
#end

# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════

Base.@kwdef mutable struct ZeoliteConfig
    # ── Input ──
    ovito_data::String    = "test.data"
    raspa_restart::String = ""
    table_file::String    = "hillsauer_nb.table"

    # ── Output ──
    output_data::String   = "loaded_zeolite.lmp"
    output_input::String  = "run_loaded.in"
    fw_check_xyz::String  = "framework_check.xyz"
    eth_check_xyz::String = "ethanol_check.xyz"

    # ── Type remapping: Ovito → LAMMPS ──
    # Ovito exports O=1 Si=2; working LAMMPS uses Si=1 O=2
    ovito_type_remap::Dict{Int,Int} = Dict(1 => 2, 2 => 1)

    # ── Framework types (after remapping) ──
    si_type::Int = 1
    o_type::Int  = 2
    al_type::Int = 0   # 0 = not present; set to 3 for aluminosilicate
    h_acid_type::Int = 0  # Brønsted H; set to e.g. 7 if present

    # ── Topology cutoffs ──
    si_o_cutoff::Float64 = 1.85
    al_o_cutoff::Float64 = 2.00
    o_h_cutoff::Float64  = 1.05

    # ── Ethanol (TraPPE) ──
    eth_atoms_per_mol::Int = 4
    eth_atom_names::Vector{String} = ["CH3", "CH2", "O_eth", "H_eth"]
    eth_types::Dict{String,Int} = Dict(
        "CH3"=>3, "CH2"=>4, "O_eth"=>5, "H_eth"=>6)
    eth_charges::Dict{String,Float64} = Dict(
        "CH3"=>0.0, "CH2"=>0.265, "O_eth"=>-0.700, "H_eth"=>0.435)
    eth_masses::Dict{String,Float64} = Dict(
        "CH3"=>15.035, "CH2"=>14.027, "O_eth"=>15.999, "H_eth"=>1.008)

    # ── Ethanol bonded topology definitions ── 
    # (bond_type, atom1_in_mol, atom2_in_mol)
    eth_bond_defs::Vector{Tuple{Int,Int,Int}} = [(2,1,2),(3,2,3),(4,3,4)]
    # (angle_type, a1, a2, a3)
    eth_angle_defs::Vector{Tuple{Int,Int,Int,Int}} = [(3,1,2,3),(4,2,3,4)]
    # (dihedral_type, a1, a2, a3, a4)
    eth_dihedral_defs::Vector{Tuple{Int,Int,Int,Int,Int}} = [(2,1,2,3,4)]

    # ── Pair style ──
    pair_cutoff::Float64 = 12.0
    coul_cutoff::Float64 = 12.0

    # ── Wrapping tolerance (Å) — molecules outside box±tol are dropped ──
    box_tolerance::Float64 = 0.5
end

# ═══════════════════════════════════════════════════════════════════
# Step 1: Read and remap Ovito framework
# ═══════════════════════════════════════════════════════════════════

function read_and_remap_framework(cfg::ZeoliteConfig)
    println("═══ Step 1: Reading Ovito framework ═══")
    data = read_lammps_data(cfg.ovito_data; verbose=false)
    natoms = size(data.coords, 1)

    # Remap types
    if !isempty(cfg.ovito_type_remap)
        println("  Remapping atom types: $(cfg.ovito_type_remap)")
        for j in 1:natoms
            old = data.atom_labels[j]
            haskey(cfg.ovito_type_remap, old) && (data.atom_labels[j] = cfg.ovito_type_remap[old])
        end
    end

    data.masses = Dict(cfg.si_type => 28.0855, cfg.o_type => 15.9994)
    cfg.al_type > 0 && (data.masses[cfg.al_type] = 26.9815)
    cfg.h_acid_type > 0 && (data.masses[cfg.h_acid_type] = 1.008)

    # Charges via bond increments (Si: +0.5236, O: -0.2618)
    for j in 1:natoms
        t = data.atom_labels[j]
        if t == cfg.si_type;       data.atom_charges[j] =  0.5236
        elseif t == cfg.o_type;    data.atom_charges[j] = -0.2618
        elseif t == cfg.al_type;   data.atom_charges[j] =  0.5236  # placeholder
        elseif t == cfg.h_acid_type; data.atom_charges[j] = 0.0
        end
    end

    nsi = count(==(cfg.si_type), data.atom_labels)
    no  = count(==(cfg.o_type),  data.atom_labels)
    nal = cfg.al_type > 0 ? count(==(cfg.al_type), data.atom_labels) : 0
    println("  Atoms: $natoms ($nsi Si, $no O" *
            (nal > 0 ? ", $nal Al" : "") * ")")

    # Stoichiometry check
    expected_o = 2 * (nsi + nal)
    no != expected_o && @warn "O count ($no) ≠ 2×T ($expected_o)"

    # Write XYZ
    open(cfg.fw_check_xyz, "w") do io
        println(io, natoms)
        println(io, "Framework (Si=$(cfg.si_type), O=$(cfg.o_type))")
        for j in 1:natoms
            t = data.atom_labels[j]
            elem = t == cfg.si_type ? "Si" : (t == cfg.al_type ? "Al" :
                   (t == cfg.o_type ? "O" : "H"))
            @printf(io, "%s %.6f %.6f %.6f\n", elem,
                    data.coords[j,1], data.coords[j,2], data.coords[j,3])
        end
    end
    println("  Wrote $(cfg.fw_check_xyz)")
    return data
end

# ═══════════════════════════════════════════════════════════════════
# Step 2: Add framework topology
# ═══════════════════════════════════════════════════════════════════

function add_framework_topology!(data, cfg::ZeoliteConfig)
    println("\n═══ Step 2: Building framework topology ═══")
    add_zeolite_topology!(data;
        si_type     = cfg.si_type,
        o_type      = cfg.o_type,
        al_type     = cfg.al_type > 0 ? cfg.al_type : nothing,
        h_type      = cfg.h_acid_type > 0 ? cfg.h_acid_type : nothing,
        si_o_cutoff = cfg.si_o_cutoff,
        al_o_cutoff = cfg.al_o_cutoff,
        o_h_cutoff  = cfg.o_h_cutoff)
    return data
end

# ═══════════════════════════════════════════════════════════════════
# Step 3: Read RASPA3 JSON restart → extract ethanol
# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════
# FIXED: read_raspa3_ethanol — wrap molecule as a unit, not per-atom
# ═══════════════════════════════════════════════════════════════════

"""
    read_raspa3_ethanol(cfg, box_dims) -> Vector of molecules

RASPA3 JSON format:  { "ethanol": [ [x,y,z], [x,y,z], ... ] }
Flat array, every `atoms_per_mol` (4) entries = 1 molecule.
Order: CH3, CH2, O_eth, H_eth.
Coordinates may be unwrapped → wrap molecule COM, keep internal geometry intact.
"""
function read_raspa3_ethanol(cfg::ZeoliteConfig, box_dims)
    println("\n═══ Step 3: Reading RASPA3 restart ═══")
    isempty(cfg.raspa_restart) && (println("  No restart file — skipping."); return [])

    restart = JSON.parsefile(cfg.raspa_restart)

    # Find adsorbate key
    #component_key = ""
    #for k in keys(restart)
    #    if restart[k] isa AbstractVector && !isempty(restart[k]) &&
    #       restart[k][1] isa AbstractVector
    #        component_key = k
    #        break
    #    end
    #end
    #isempty(component_key) && error("No adsorbate positions found in JSON")
    #println("  Component: \"$component_key\"")
    #all_positions = restart[component_key]

  

    # find adsorbate key new method - this is done to avoid crash in cases where there is no 
    # adsornate. Also, this is to allow me to track bugs and errors in the simulations
    # that might comes from bad ff parameters, mistakes in the input files, etc.

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
    all_positions = restart[component_key]
    if isempty(all_positions)
        println("  WARNING: 0 molecules — GCMC produced empty loading")
        return []
    end

    npos = length(all_positions)
    apm = cfg.eth_atoms_per_mol   # atoms per molecule

    npos % apm != 0 &&
        error("$npos positions not divisible by $apm atoms/mol")
    nmols_total = npos ÷ apm
    println("  Raw: $nmols_total molecules ($npos positions)")

    # Box
    Lx = box_dims[1,2] - box_dims[1,1]
    Ly = box_dims[2,2] - box_dims[2,1]
    Lz = box_dims[3,2] - box_dims[3,1]
    L  = [Lx, Ly, Lz]
    lo = [box_dims[1,1], box_dims[2,1], box_dims[3,1]]

    molecules = Vector{Vector{Tuple{String,Float64,Float64,Float64}}}()
    n_outside = 0

    for mi in 1:nmols_total
        # ── Collect raw positions for this molecule ──
        raw = Matrix{Float64}(undef, apm, 3)
        for k in 1:apm
            p = all_positions[(mi-1)*apm + k]
            raw[k, :] = [Float64(p[1]), Float64(p[2]), Float64(p[3])]
        end

        # ── Wrap as a rigid unit ──
        # 1. Compute center of mass (unweighted, all UA sites ~same mass)
        com = vec(mean(raw, dims=1))

        # 2. Wrap COM into box [lo, lo+L)
        com_wrapped = copy(com)
        for d in 1:3
            com_wrapped[d] = lo[d] + mod(com[d] - lo[d], L[d])
        end

        # 3. Shift ALL atoms by the same displacement
        shift = com_wrapped - com
        wrapped = raw .+ shift'

        # 4. Verify internal geometry is preserved (bonds < 2 Å)
        ok = true
        for k in 1:(apm-1)
            bond_len = sqrt(sum((wrapped[k,:] - wrapped[k+1,:]).^2))
            if bond_len > 3.0
                @warn "Molecule $mi: bond $k-$(k+1) = $(round(bond_len,digits=2)) Å after wrapping"
                ok = false
            end
        end

        # 5. Check that COM is inside box (with tolerance)
        tol = cfg.box_tolerance
        inside = all(lo[d] - tol <= com_wrapped[d] <= lo[d] + L[d] + tol for d in 1:3)

        if !inside || !ok
            n_outside += 1
            continue
        end

        # Build molecule
        mol = Tuple{String,Float64,Float64,Float64}[]
        for k in 1:apm
            push!(mol, (cfg.eth_atom_names[k],
                        wrapped[k,1], wrapped[k,2], wrapped[k,3]))
        end
        push!(molecules, mol)
    end

    println("  Kept: $(length(molecules)) molecules inside box")
    n_outside > 0 && println("  Removed $n_outside molecules (outside box or broken)")

    # Write XYZ check
    if !isempty(molecules)
        ntot = sum(length(m) for m in molecules)
        open(cfg.eth_check_xyz, "w") do io
            println(io, ntot)
            println(io, "Ethanol ($(length(molecules)) molecules, COM-wrapped)")
            for mol in molecules
                for (name, x, y, z) in mol
                    elem = name == "H_eth" ? "H" : (name == "O_eth" ? "O" : "C")
                    @printf(io, "%s %.6f %.6f %.6f\n", elem, x, y, z)
                end
            end
        end
        println("  Wrote $(cfg.eth_check_xyz)")
    end

    return molecules
end


# ═══════════════════════════════════════════════════════════════════
# Step 4: Merge framework + ethanol
# ═══════════════════════════════════════════════════════════════════

function merge_framework_ethanol!(fw, molecules, cfg::ZeoliteConfig)
    println("\n═══ Step 4: Merging framework + ethanol ═══")
    nmols = length(molecules)
    # nmols == 0 && (println("  No ethanol."); return fw)

    # NEW (register masses before returning):
    if  nmols == 0
        println("  No ethanol (0 molecules).")
        for (name, t) in cfg.eth_types; fw.masses[t] = cfg.eth_masses[name]; end
        return fw
    end


    nfw = size(fw.coords, 1)
    neth_per = cfg.eth_atoms_per_mol
    neth = nmols * neth_per
    ntot = nfw + neth
    max_mol = maximum(fw.molecule_labels; init=0)

    # Extend arrays
    new_coords  = zeros(ntot, 3);   new_coords[1:nfw, :] = fw.coords
    new_ids     = zeros(Int, ntot); new_ids[1:nfw] = fw.atom_ids
    new_types   = zeros(Int, ntot); new_types[1:nfw] = fw.atom_labels
    new_q       = zeros(ntot);      new_q[1:nfw] = fw.atom_charges
    new_mols    = zeros(Int, ntot); new_mols[1:nfw] = fw.molecule_labels
    new_img     = zeros(Int, ntot, 3)
    !isempty(fw.image_flags) && size(fw.image_flags,1)==nfw &&
        (new_img[1:nfw,:] = fw.image_flags)

    for (mi, mol) in enumerate(molecules)
        mol_id = max_mol + mi
        for (k, (name, x, y, z)) in enumerate(mol)
            idx = nfw + (mi-1)*neth_per + k
            new_ids[idx]       = idx
            new_types[idx]     = cfg.eth_types[name]
            new_q[idx]         = cfg.eth_charges[name]
            new_mols[idx]      = mol_id
            new_coords[idx, :] = [x, y, z]
        end
    end

    fw.coords          = new_coords
    fw.atom_ids        = new_ids
    fw.atom_labels     = new_types
    fw.atom_charges    = new_q
    fw.molecule_labels = new_mols
    fw.image_flags     = new_img

    # ── Ethanol bonds: CH3-CH2(2), CH2-O(3), O-H(4) ──
    eth_bdefs =  cfg.eth_bond_defs   # [(2,1,2), (3,2,3), (4,3,4)]
    nfb = size(fw.bonds, 1)
    neb = nmols * length(eth_bdefs)    #3
    nb  = zeros(Int, nfb+neb, 2);  nb[1:nfb,:] = fw.bonds
    nbl = zeros(Int, nfb+neb);     nbl[1:nfb] = fw.bond_labels
    for mi in 1:nmols
        base = nfw + (mi-1)*neth_per
        for (bi,(bt,a1,a2)) in enumerate(eth_bdefs)
            i = nfb + (mi-1)*length(eth_bdefs) + bi               #3 + bi
            nbl[i] = bt; nb[i,1] = base+a1; nb[i,2] = base+a2
        end
    end
    fw.bonds = nb; fw.bond_labels = nbl
    fw.nbond_types = max(fw.nbond_types, maximum(bt for (bt,_,_) in eth_bdefs)) #; fw.nbond_types = max(fw.nbond_types, 4)

    # ── Ethanol angles: CH3-CH2-O(3), CH2-O-H(4) ──
    eth_adefs =  cfg.eth_angle_defs        # [(3,1,2,3), (4,2,3,4)]
    nfa = size(fw.angles, 1)
    nea = nmols *  length(eth_adefs)    # 2
    na  = zeros(Int, nfa+nea, 3); na[1:nfa,:] = fw.angles
    nal = zeros(Int, nfa+nea);    nal[1:nfa] = fw.angle_labels
    for mi in 1:nmols
        base = nfw + (mi-1)*neth_per
        for (ai,(at,a1,a2,a3)) in enumerate(eth_adefs)
            i = nfa + (mi-1)*length(eth_adefs) + ai     #2 + ai
            nal[i] = at; na[i,:] = [base+a1, base+a2, base+a3]
        end
    end
    #fw.angles = na; fw.angle_labels = nal; fw.nangle_types = max(fw.nangle_types, 4)
    fw.angles = na; fw.angle_labels = nal
    fw.nangle_types = max(fw.nangle_types, maximum(at for (at,_,_,_) in eth_adefs))


    # ── Ethanol dihedrals: CH3-CH2-O-H (type 2) ──
    eth_ddefs = cfg.eth_dihedral_defs
    nfd = size(fw.dihedrals, 1)
    ned = nmols * length(eth_ddefs)
    nd  = zeros(Int, nfd+ned, 4); nd[1:nfd,:] = fw.dihedrals
    ndl = zeros(Int, nfd+ned);    ndl[1:nfd] = fw.dihedral_labels
    for mi in 1:nmols
        base = nfw + (mi-1)*neth_per
        for (di,(dt,a1,a2,a3,a4)) in enumerate(eth_ddefs)
            i = nfd + (mi-1)*length(eth_ddefs) + di
            ndl[i] = dt; nd[i,:] = [base+a1, base+a2, base+a3, base+a4]
        end
    end
    fw.dihedrals = nd; fw.dihedral_labels = ndl
    fw.ndihedral_types = max(fw.ndihedral_types, maximum(dt for (dt,_,_,_,_) in eth_ddefs))

    
    #=
    nfd = size(fw.dihedrals, 1)
    ned = nmols
    nd  = zeros(Int, nfd+ned, 4); nd[1:nfd,:] = fw.dihedrals
    ndl = zeros(Int, nfd+ned);    ndl[1:nfd] = fw.dihedral_labels
    for mi in 1:nmols
        base = nfw + (mi-1)*neth_per
        i = nfd + mi
        ndl[i] = 2; nd[i,:] = [base+1, base+2, base+3, base+4]
    end
    fw.dihedrals = nd; fw.dihedral_labels = ndl; fw.ndihedral_types = max(fw.ndihedral_types, 2)
    =#

    # Update masses
    for (name, t) in cfg.eth_types; fw.masses[t] = cfg.eth_masses[name]; end

    qtot = sum(fw.atom_charges)
    println("  Total: $ntot atoms ($nfw fw + $neth ethanol, $nmols molecules)")
    @printf("  Charge: %.6f e\n", qtot)
    abs(qtot) > 0.01 && @warn "NOT charge-neutral: Q=$qtot"

    return fw
end

# ═══════════════════════════════════════════════════════════════════
# Step 5: Write LAMMPS data file with embedded class2 coefficients
# adding impropers
# ═══════════════════════════════════════════════════════════════════


# function that will generate the impropers

"""
    generate_impropers!(data; si_type=1, o_type=2)

Generate class2 improper terms for each Si atom.
For Si with 4 O neighbors (O1,O2,O3,O4), generates 4 impropers:
  O1-Si-O2-O3, O1-Si-O2-O4, O1-Si-O3-O4, O2-Si-O3-O4
All are improper type 1.

Must be called AFTER add_zeolite_topology! (needs bond list).
"""
function generate_impropers!(data; si_type::Int=1, o_type::Int=2)
    natoms = size(data.coords, 1)

    # Build neighbor list
    neighbors = Dict{Int, Vector{Int}}()
    for j in 1:natoms
        neighbors[data.atom_ids[j]] = Int[]
    end
    for k in 1:size(data.bonds, 1)
        a1, a2 = data.bonds[k, 1], data.bonds[k, 2]
        push!(neighbors[a1], a2)
        push!(neighbors[a2], a1)
    end

    id_to_idx = Dict(data.atom_ids[j] => j for j in 1:natoms)

    impropers = Vector{NTuple{4,Int}}()

    for j in 1:natoms
        data.atom_labels[j] == si_type || continue
        si_id = data.atom_ids[j]
        o_nbrs = filter(n -> data.atom_labels[id_to_idx[n]] == o_type, neighbors[si_id])
        length(o_nbrs) == 4 || continue

        # Generate C(4,3) = 4 impropers: Oa-Si-Ob-Oc
        for a in 1:4, b in (a+1):4, c in (b+1):4
            # LAMMPS class2 improper: i-j-k-l where j=center
            # Pick the remaining O as i, center=Si, k and l from the triple
            others = [o_nbrs[a], o_nbrs[b], o_nbrs[c]]
            push!(impropers, (others[1], si_id, others[2], others[3]))
        end
    end

    n_imp = length(impropers)
    if n_imp > 0
        data.nimproper_types = 1
        # Store as matrix + labels (matching bonds/angles/dihedrals convention)
        imp_mat = zeros(Int, n_imp, 4)
        imp_labels = ones(Int, n_imp)
        for k in 1:n_imp
            imp_mat[k, 1] = impropers[k][1]
            imp_mat[k, 2] = impropers[k][2]
            imp_mat[k, 3] = impropers[k][3]
            imp_mat[k, 4] = impropers[k][4]
        end
        data.impropers = imp_mat
        data.improper_labels = imp_labels
    end

    println("  Impropers: $n_imp ($(n_imp÷4) Si atoms × 4 each)")
    return data
end


function write_complete_data(fname::String, d, cfg::ZeoliteConfig)
    println("\n═══ Step 5: Writing $fname ═══")

    natoms = size(d.coords, 1)
    nbt = d.nbond_types; nat = d.nangle_types; ndt = d.ndihedral_types
    natypes = maximum(keys(d.masses))
    is_tri = any(d.tilt_factors .!= 0.0)

    has_eth = natypes > 2
    has_imp = hasproperty(d, :impropers) && size(d.impropers, 1) > 0
    n_imp = has_imp ? size(d.impropers, 1) : 0
    n_imp_types = has_imp ? d.nimproper_types : 0

    open(fname, "w") do io
        println(io, "LAMMPS data — loaded zeolite (build_loaded_zeolite.jl v1.1)\n")
        println(io, "$(natoms) atoms")
        println(io, "$(size(d.bonds,1)) bonds")
        println(io, "$(size(d.angles,1)) angles")
        println(io, "$(size(d.dihedrals,1)) dihedrals")
        println(io, "$(n_imp) impropers\n")
        println(io, "$natypes atom types")
        println(io, "$nbt bond types")
        println(io, "$nat angle types")
        println(io, "$ndt dihedral types")
        println(io, "$(n_imp_types) improper types\n")
        @printf(io, "%.10f %.10f xlo xhi\n", d.box_dimensions[1,1], d.box_dimensions[1,2])
        @printf(io, "%.10f %.10f ylo yhi\n", d.box_dimensions[2,1], d.box_dimensions[2,2])
        @printf(io, "%.10f %.10f zlo zhi\n", d.box_dimensions[3,1], d.box_dimensions[3,2])
        is_tri && @printf(io, "%.10f %.10f %.10f xy xz yz\n",
                          d.tilt_factors[1], d.tilt_factors[2], d.tilt_factors[3])

        mn = Dict(1=>"Si",2=>"O",3=>"CH3_eth",4=>"CH2_eth",5=>"O_eth",6=>"H_eth")
        println(io, "\nMasses\n")
        for t in sort(collect(keys(d.masses)))
            @printf(io, "  %d  %.6f  # %s\n", t, d.masses[t], get(mn,t,""))
        end

        # ── Bond Coeffs ──
        println(io, "\nBond Coeffs # class2\n")
        println(io, "  1  1.6104  459.0786  -672.4445  443.3651  # Si-O")
        if has_eth
            println(io, "  2  1.540  553.94  0.0  0.0  # CH3-CH2")
            println(io, "  3  1.430  553.94  0.0  0.0  # CH2-O_eth")
            println(io, "  4  0.945  553.94  0.0  0.0  # O_eth-H_eth")
        end

        # ── Angle Coeffs ──
        println(io, "\nAngle Coeffs # class2\n")
        println(io, "  1  150.0  20.7015  27.5506  10.9930  # Si-O-Si")
        println(io, "  2  113.0  81.9691  -36.5814  116.9558  # O-Si-O")
        if has_eth
            println(io, "  3  109.47  100.16  0.0  0.0  # CH3-CH2-O_eth")
            println(io, "  4  108.50  110.09  0.0  0.0  # CH2-O_eth-H_eth")
        end

        println(io, "\nBondBond Coeffs\n")
        println(io, "  1  151.8742  1.6104  1.6104")
        println(io, "  2  0.0  1.6104  1.6104")
        has_eth && println(io, "  3  0.0  1.540  1.430\n  4  0.0  1.430  0.945")

        println(io, "\nBondAngle Coeffs\n")
        println(io, "  1  9.2390  9.2390  1.6104  1.6104")
        println(io, "  2  78.1239  78.1239  1.6104  1.6104")
        has_eth && println(io, "  3  0.0  0.0  1.540  1.430\n  4  0.0  0.0  1.430  0.945")

        println(io, "\nDihedral Coeffs # class2\n")
        println(io, "  1  0.0306  0.0  -0.0105  0.0  0.0804  0.0  # O-Si-O-Si")
        has_eth && println(io, "  2  0.4169  0.0  0.0580  0.0  0.3734  0.0  # CH3-CH2-O-H")

        println(io, "\nMiddleBondTorsion Coeffs\n")
        println(io, "  1  0.0  0.0  0.0  1.6104")
        has_eth && println(io, "  2  0.0  0.0  0.0  1.430")

        println(io, "\nEndBondTorsion Coeffs\n")
        println(io, "  1  0.0  0.0  0.0  0.0  0.0  0.0  1.6104  1.6104")
        has_eth && println(io, "  2  0.0  0.0  0.0  0.0  0.0  0.0  1.540  0.945")

        println(io, "\nAngleTorsion Coeffs\n")
        println(io, "  1  0.0  0.0  0.0  0.0  0.0  0.0  113.0  150.0")
        has_eth && println(io, "  2  0.0  0.0  0.0  0.0  0.0  0.0  109.47  108.50")

        println(io, "\nAngleAngleTorsion Coeffs\n")
        println(io, "  1  -4.5150  113.0  150.0")
        has_eth && println(io, "  2  0.0  109.47  108.50")

        println(io, "\nBondBond13 Coeffs\n")
        println(io, "  1  0.0  1.6104  1.6104")
        has_eth && println(io, "  2  0.0  1.540  0.945")

        # ── Improper Coeffs (if impropers exist) ──
        if has_imp
            println(io, "\nImproper Coeffs # class2\n")
            println(io, "  1  0.0  0.0  # chi0=0, K0=0")

            println(io, "\nAngleAngle Coeffs\n")
            # M1, M2, M3, θ1, θ2, θ3
            # From Sholl: M1=-6.303, M2=0, M3=0, θ1=θ3=θ₀(O-Si-O)
            # Using original HSFF θ₀ = 112.02 (same convention as aat)
            println(io, "  1  -6.3030  0.0  0.0  112.0200  0.0  112.0200")
        end


        # ── Atoms ──
        println(io, "\nAtoms # full\n")
        for j in 1:natoms
            @printf(io, "  %d %d %d %.6f %.10f %.10f %.10f %d %d %d\n",
                    d.atom_ids[j], d.molecule_labels[j], d.atom_labels[j],
                    d.atom_charges[j], d.coords[j,1], d.coords[j,2], d.coords[j,3],
                    d.image_flags[j,1], d.image_flags[j,2], d.image_flags[j,3])
        end
### adcionar aqui <-


        # ── Topology sections ──
        for (secname, labels, atoms) in [
            ("Bonds", d.bond_labels, d.bonds),
            ("Angles", d.angle_labels, d.angles),
            ("Dihedrals", d.dihedral_labels, d.dihedrals)]
            size(atoms,1) == 0 && continue
            println(io, "\n$secname\n")
            for k in 1:size(atoms,1)
                print(io, "  $k $(labels[k])")
                for c in 1:size(atoms,2); print(io, " $(atoms[k,c])"); end
                println(io)
            end
        end
        
        # ── Impropers ──
        if has_imp
            println(io, "\nImpropers\n")
            for k in 1:n_imp
                print(io, "  $k $(d.improper_labels[k])")
                for c in 1:4; print(io, " $(d.impropers[k,c])"); end
                println(io)
            end
        end

    end
    println("  Done: $fname")
end

# ═══════════════════════════════════════════════════════════════════
# Step 6: Write LAMMPS input script
# ═══════════════════════════════════════════════════════════════════

function write_input_script(fname::String, data_file::String, cfg::ZeoliteConfig)
    println("\n═══ Step 6: Writing $fname ═══")
    open(fname, "w") do io
        print(io, """
# ============================================================
# LAMMPS: Loaded MFI zeolite — Hill-Sauer + TraPPE
# Generated by build_loaded_zeolite.jl v1.1
# Types: 1=Si 2=O 3=CH3 4=CH2 5=O_eth 6=H_eth
# ============================================================

units           real
boundary        p p p
atom_style      full

# Class2 cross-term coefficients are in the data file
bond_style      class2
angle_style     class2
dihedral_style  class2
improper_style  class2

read_data       $data_file

# ============================================================
# PAIR STYLE
# ============================================================
# fw-fw: table A/r^9 (Hill-Sauer)
# fw-ads + ads-ads: lj/cut (Bai x TraPPE)
# Coulomb: Ewald

pair_style      hybrid/overlay table spline 2000 lj/cut $(cfg.pair_cutoff) coul/long $(cfg.coul_cutoff)
pair_modify     tail no shift no
kspace_style    ewald 1.0E-06

pair_coeff      * * coul/long

# fw-fw: tabulated
pair_coeff  1 1  table $(cfg.table_file) Si_Si
pair_coeff  1 2  table $(cfg.table_file) Si_O
pair_coeff  2 2  table $(cfg.table_file) O_O

# fw-ads (Bai x TraPPE, Lorentz-Berthelot)
pair_coeff  1 3  lj/cut  0.09227  3.025   # Si-CH3
pair_coeff  1 4  lj/cut  0.06321  3.125   # Si-CH2
pair_coeff  1 5  lj/cut  0.08988  2.660   # Si-O_eth
pair_coeff  1 6  lj/cut  0.0      1.000   # Si-H_eth
pair_coeff  2 3  lj/cut  0.14320  3.525   # O-CH3
pair_coeff  2 4  lj/cut  0.09811  3.625   # O-CH2
pair_coeff  2 5  lj/cut  0.13949  3.160   # O-O_eth
pair_coeff  2 6  lj/cut  0.0      1.000   # O-H_eth

# ads-ads (TraPPE)
pair_coeff  3 3  lj/cut  0.19475  3.750   # CH3-CH3
pair_coeff  3 4  lj/cut  0.13341  3.850   # CH3-CH2
pair_coeff  3 5  lj/cut  0.18975  3.385   # CH3-O_eth
pair_coeff  3 6  lj/cut  0.0      1.000   # CH3-H_eth
pair_coeff  4 4  lj/cut  0.09141  3.950   # CH2-CH2
pair_coeff  4 5  lj/cut  0.12994  3.485   # CH2-O_eth
pair_coeff  4 6  lj/cut  0.0      1.000   # CH2-H_eth
pair_coeff  5 5  lj/cut  0.18481  3.020   # O_eth-O_eth
pair_coeff  5 6  lj/cut  0.0      1.000   # O_eth-H_eth
pair_coeff  6 6  lj/cut  0.0      1.000   # H_eth-H_eth

# ============================================================
special_bonds   lj 0.0 0.0 1.0 coul 0.0 0.0 1.0 angle yes dihedral yes
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes one 5000 page 150000

group           framework type 1 2
group           adsorbate type 3 4 5 6

variable        T equal 300.0
thermo          1000
thermo_style    custom step press temp etotal enthalpy vol &
                cella cellb cellc cellalpha cellbeta cellgamma

# ============================================================
# PHASE 1: Minimize (framework frozen)
# ============================================================
fix freeze framework setforce 0.0 0.0 0.0
velocity framework set 0.0 0.0 0.0
min_modify dmax 0.02
min_style cg
minimize 1e-6 1e-8 5000 50000
unfix freeze

# ============================================================
# PHASE 2: NVT ramp (10 -> 300 K)
# ============================================================
reset_timestep 0
fix freeze framework setforce 0.0 0.0 0.0
velocity adsorbate create 10.0 12345 dist gaussian mom yes rot yes
timestep 0.25

fix eq1 adsorbate nvt temp 10.0 50.0 100.0
run 400
unfix eq1

fix eq2 adsorbate nvt temp 50.0 150.0 100.0
run 800
unfix eq2

fix eq3 adsorbate nvt temp 150.0 300.0 100.0
run 1200
unfix eq3
unfix freeze

fix hold framework spring/self 10.0
fix eq4 all nvt temp \${T} \${T} 100.0
run 2000
unfix eq4
unfix hold

# ============================================================
# PHASE 3: NPT production
# ============================================================
reset_timestep 0

change_box all triclinic
fix prod all npt temp \${T} \${T} 100.0 tri 1.0 1.0 5000.0
fix mom all momentum 1000 linear 1 1 1 angular

dump d1 all custom 1000 traj_loaded.lammpstrj id type x y z
dump_modify d1 sort id
run 10000

write_data loaded_npt_final.lmp
print "=== Complete ==="
""")
    end
    println("  Done: $fname")
end


#=
#
# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

function main(args=ARGS)
    println("╔══════════════════════════════════════════════╗")
    println("║  build_loaded_zeolite.jl v1.1               ║")
    println("╚══════════════════════════════════════════════╝\n")

    cfg = ZeoliteConfig(
        ovito_data    = get(args, 1, "test.data"),
        raspa_restart = get(args, 2, ""),
        output_data   = get(args, 3, "loaded_zeolite.lmp"),
    )

    fw = read_and_remap_framework(cfg)
    add_framework_topology!(fw, cfg)
    mols = read_raspa3_ethanol(cfg, fw.box_dimensions)
    merge_framework_ethanol!(fw, mols, cfg)
    write_complete_data(cfg.output_data, fw, cfg)
    write_input_script(cfg.output_input, cfg.output_data, cfg)

    println("\n═══ Pipeline complete ═══")
    println("  $(cfg.output_data) + $(cfg.output_input)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
=#
