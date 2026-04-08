#!/usr/bin/env julia
"""
add_zeolite_topology.jl  —  v3

Add full bonded topology (bonds, angles, dihedrals) to a LAMMPS data file
exported from Ovito for zeolite simulations.

Integrates with the LammpsDataReader module (read_lammps_data.jl).

─────────────────────────────────────────────────────────────────────
  TOPOLOGY TYPE MAP  (Si and Al fully distinguished)
─────────────────────────────────────────────────────────────────────
  Bonds:
    1 = Si─O
    2 = Al─O
    3 = O─H   (Brønsted)

  Angles:
    1 = Si─O─Si
    2 = O─Si─O
    3 = Si─O─Al
    4 = O─Al─O
    5 = Si─O─H
    6 = Al─O─H

  Dihedrals:
    1 = O─Si─O─Si
    2 = O─Si─O─Al
    3 = O─Al─O─Si
    4 = O─Al─O─Al
    5 = O─Si─O─H
    6 = O─Al─O─H

  Impropers: (none)
─────────────────────────────────────────────────────────────────────
  For pure silica (no Al/H), only types 1/2/1 appear in bonds/
  angles/dihedrals — fully backward-compatible with v2.
─────────────────────────────────────────────────────────────────────

Atom type conventions (Ovito/CIF):
    1 = O, 2 = Si, 3 = Al (optional), 4 = H_acid (optional)

Usage standalone:
    julia add_zeolite_topology.jl input.data output.data [options]
    julia add_zeolite_topology.jl input.data output.data --al-type 3 --h-type 4
    

Usage from module:
    include("read_lammps_data.jl")
    include("add_zeolite_topology.jl")
    data = LammpsDataReader.read_lammps_data("test.data")
    add_zeolite_topology!(data; si_type=2, o_type=1, al_type=3, h_type=4)
    LammpsDataReader.write_lammps_data("output.data", data)
"""


# ============================================================================
# Try to load the module if available
# ============================================================================
#if !@isdefined(LammpsDataReader)
#    local_module = joinpath(@__DIR__, "read_lammps_data.jl")
#    if isfile(local_module)
#        include(local_module)
#    end
#end

# ============================================================================
# Minimum-image distance (orthogonal + triclinic)
# ============================================================================

function min_image_distance(r1, r2, box_dims, tilt)
    lx = box_dims[1,2] - box_dims[1,1]
    ly = box_dims[2,2] - box_dims[2,1]
    lz = box_dims[3,2] - box_dims[3,1]
    xy, xz, yz = tilt

    dx = r1[1] - r2[1]
    dy = r1[2] - r2[2]
    dz = r1[3] - r2[3]

    while dz >  lz/2; dz -= lz; dy -= yz; dx -= xz; end
    while dz < -lz/2; dz += lz; dy += yz; dx += xz; end
    while dy >  ly/2; dy -= ly; dx -= xy; end
    while dy < -ly/2; dy += ly; dx += xy; end
    while dx >  lx/2; dx -= lx; end
    while dx < -lx/2; dx += lx; end

    return sqrt(dx^2 + dy^2 + dz^2)
end

# ============================================================================
# Type classification helpers
# ============================================================================

function _bond_type(t_atom_type, si_type, al_type)
    t_atom_type == si_type && return 1   # Si-O
    t_atom_type == al_type && return 2   # Al-O
    error("Unknown T-atom type: $t_atom_type")
end

function _angle_type_TOT(t1_type, t2_type, si_type, al_type)
    a, b = minmax(t1_type, t2_type)
    a == si_type && b == si_type && return 1   # Si-O-Si
    a == si_type && b == al_type && return 3   # Si-O-Al
    a == al_type && b == al_type && return 3   # Al-O-Al (rare, grouped with Si-O-Al)
    error("Unknown T-O-T: $t1_type, $t2_type")
end

function _angle_type_OTO(center_type, si_type, al_type)
    center_type == si_type && return 2   # O-Si-O
    center_type == al_type && return 4   # O-Al-O
    error("Unknown O-T-O center: $center_type")
end

function _angle_type_TOH(t_type, si_type, al_type)
    t_type == si_type && return 5   # Si-O-H
    t_type == al_type && return 6   # Al-O-H
    error("Unknown T-O-H: $t_type")
end

function _dihedral_type_OTOT(t_inner, t_outer, si_type, al_type)
    t_inner == si_type && t_outer == si_type && return 1   # O-Si-O-Si
    t_inner == si_type && t_outer == al_type && return 2   # O-Si-O-Al
    t_inner == al_type && t_outer == si_type && return 3   # O-Al-O-Si
    t_inner == al_type && t_outer == al_type && return 4   # O-Al-O-Al
    error("Unknown dihedral T types: $t_inner, $t_outer")
end

function _dihedral_type_OTOH(t_inner, si_type, al_type)
    t_inner == si_type && return 5   # O-Si-O-H
    t_inner == al_type && return 6   # O-Al-O-H
    error("Unknown Bronsted dihedral T: $t_inner")
end


# ============================================================================
# Main topology builder
# ============================================================================

"""
    add_zeolite_topology!(data; kwargs...)

Add bonds, angles, and dihedrals to a LammpsData object in-place.
Si and Al are fully distinguished in all topology types.

# Keyword arguments
- `si_type=2`:       atom type for Si
- `o_type=1`:        atom type for framework O
- `al_type=nothing`: atom type for Al (nothing or 0 to skip)
- `h_type=nothing`:  atom type for Bronsted H (nothing or 0 to skip)
- `si_o_cutoff=1.85`: Si-O bond cutoff (A)
- `al_o_cutoff=2.00`: Al-O bond cutoff (A)
- `o_h_cutoff=1.05`:  O-H bond cutoff (A)
- `verbose=true`:     print topology summary
"""
function add_zeolite_topology!(data;
        si_type::Int    = 2,
        o_type::Int     = 1,
        al_type         = nothing,
        h_type          = nothing,
        si_o_cutoff     = 1.85,
        al_o_cutoff     = 2.00,
        o_h_cutoff      = 1.05,
        verbose::Bool   = true)

    natoms = size(data.coords, 1)
    box    = data.box_dimensions
    tilt   = data.tilt_factors

    has_al = al_type !== nothing && al_type > 0
    has_h  = h_type  !== nothing && h_type  > 0
    _al = has_al ? al_type : -1
    _h  = has_h  ? h_type  : -1

    t_types = Set{Int}([si_type])
    has_al && push!(t_types, _al)

    type_of = data.atom_labels

    t_indices = [j for j in 1:natoms if type_of[j] in t_types]
    o_indices = [j for j in 1:natoms if type_of[j] == o_type]
    h_indices = has_h ? [j for j in 1:natoms if type_of[j] == _h] : Int[]

    verbose && println("  Framework: $(length(t_indices)) T-atoms ",
        "($(count(j->type_of[j]==si_type, t_indices)) Si",
        has_al ? ", $(count(j->type_of[j]==_al, t_indices)) Al" : "",
        "), $(length(o_indices)) O",
        has_h ? ", $(length(h_indices)) H" : "")

    # ── STEP 1: Bonds ─────────────────────────────────────────────
    neighbors = [Tuple{Int,Int}[] for _ in 1:natoms]
    bond_list = Tuple{Int,Int,Int}[]

    for ti in t_indices
        cutoff = type_of[ti] == si_type ? si_o_cutoff : al_o_cutoff
        ri = @view data.coords[ti, :]
        bt = _bond_type(type_of[ti], si_type, _al)
        for oi in o_indices
            ro = @view data.coords[oi, :]
            r = min_image_distance(ri, ro, box, tilt)
            if r < cutoff
                push!(bond_list, (bt, ti, oi))
                push!(neighbors[ti], (oi, bt))
                push!(neighbors[oi], (ti, bt))
            end
        end
    end

    for hi in h_indices
        rh = @view data.coords[hi, :]
        best_oi = 0; best_r = o_h_cutoff
        for oi in o_indices
            ro = @view data.coords[oi, :]
            r = min_image_distance(rh, ro, box, tilt)
            if r < best_r; best_r = r; best_oi = oi; end
        end
        if best_oi > 0
            push!(bond_list, (3, best_oi, hi))
            push!(neighbors[best_oi], (hi, 3))
            push!(neighbors[hi], (best_oi, 3))
        else
            verbose && @warn "Bronsted H idx=$hi: no O within $o_h_cutoff A"
        end
    end

    # ── Coordination report ──
    if verbose
        si_idxs = [j for j in 1:natoms if type_of[j] == si_type]
        if !isempty(si_idxs)
            c = [count(nb -> type_of[nb[1]] == o_type, neighbors[j]) for j in si_idxs]
            @printf("  Si coordination: min=%d max=%d avg=%.2f (expect 4)\n",
                    minimum(c), maximum(c), sum(c)/length(c))
        end
        if has_al
            al_idxs = [j for j in 1:natoms if type_of[j] == _al]
            if !isempty(al_idxs)
                c = [count(nb -> type_of[nb[1]] == o_type, neighbors[j]) for j in al_idxs]
                @printf("  Al coordination: min=%d max=%d avg=%.2f (expect 4)\n",
                        minimum(c), maximum(c), sum(c)/length(c))
            end
        end
        c = [count(nb -> type_of[nb[1]] in t_types, neighbors[j]) for j in o_indices]
        if !isempty(c)
            @printf("  O  coordination: min=%d max=%d avg=%.2f (expect 2)\n",
                    minimum(c), maximum(c), sum(c)/length(c))
            n1 = count(==(1), c)
            n1 > 0 && @warn "$n1 terminal O atoms (bonded to only 1 T-atom)"
        end

        for hi in h_indices
            o_nbrs = [nb[1] for nb in neighbors[hi] if type_of[nb[1]] == o_type]
            length(o_nbrs) != 1 && (@warn "H idx=$hi bonded to $(length(o_nbrs)) O"; continue)
            oi = o_nbrs[1]
            t_nbrs = [type_of[nb[1]] for nb in neighbors[oi] if type_of[nb[1]] in t_types]
            if has_al
                if any(==(al_type), t_nbrs) && any(==(si_type), t_nbrs)
                    println("  H idx=$hi bridges Al-O-Si (correct Bronsted site)")
                elseif !any(==(al_type), t_nbrs)
                    @warn "H idx=$hi: bridging O has no Al neighbor (T-types=$t_nbrs)"
                end
            end
        end
    end

    # ── STEP 2: Angles ────────────────────────────────────────────
    angle_list = Tuple{Int,Int,Int,Int}[]

    for center in 1:natoms
        nbrs = neighbors[center]
        length(nbrs) < 2 && continue
        ct = type_of[center]

        for i in 1:length(nbrs)
            for j in (i+1):length(nbrs)
                n1, n2 = nbrs[i][1], nbrs[j][1]
                t1, t2 = type_of[n1], type_of[n2]

                if ct == o_type
                    if t1 in t_types && t2 in t_types
                        at = _angle_type_TOT(t1, t2, si_type, _al)
                        push!(angle_list, (at, n1, center, n2))
                    elseif t1 in t_types && t2 == _h
                        at = _angle_type_TOH(t1, si_type, _al)
                        push!(angle_list, (at, n1, center, n2))
                    elseif t1 == _h && t2 in t_types
                        at = _angle_type_TOH(t2, si_type, _al)
                        push!(angle_list, (at, n2, center, n1))
                    end
                elseif ct in t_types
                    if t1 == o_type && t2 == o_type
                        at = _angle_type_OTO(ct, si_type, _al)
                        push!(angle_list, (at, n1, center, n2))
                    end
                end
            end
        end
    end

    # ── STEP 3: Dihedrals ─────────────────────────────────────────
    dihedral_list = Tuple{Int,Int,Int,Int,Int}[]
    seen = Set{NTuple{4,Int}}()

    for center_t in 1:natoms
        type_of[center_t] in t_types || continue
        ct = type_of[center_t]

        o_nbrs = [nb[1] for nb in neighbors[center_t] if type_of[nb[1]] == o_type]

        for oi in o_nbrs
            far_t = [nb[1] for nb in neighbors[oi]
                     if type_of[nb[1]] in t_types && nb[1] != center_t]
            far_h = has_h ? [nb[1] for nb in neighbors[oi]
                            if type_of[nb[1]] == _h] : Int[]

            for oj in o_nbrs
                oj == oi && continue

                for tf in far_t
                    key = (oj, center_t, oi, tf)
                    rkey = (tf, oi, center_t, oj)
                    if !(key in seen) && !(rkey in seen)
                        push!(seen, key)
                        dt = _dihedral_type_OTOT(ct, type_of[tf], si_type, _al)
                        push!(dihedral_list, (dt, oj, center_t, oi, tf))
                    end
                end

                for hf in far_h
                    key = (oj, center_t, oi, hf)
                    rkey = (hf, oi, center_t, oj)
                    if !(key in seen) && !(rkey in seen)
                        push!(seen, key)
                        dt = _dihedral_type_OTOH(ct, si_type, _al)
                        push!(dihedral_list, (dt, oj, center_t, oi, hf))
                    end
                end
            end
        end
    end

    # ── STEP 4: Store ─────────────────────────────────────────────
    nbonds     = length(bond_list)
    nangles    = length(angle_list)
    ndihedrals = length(dihedral_list)

    nbond_types     = nbonds > 0     ? maximum(b[1] for b in bond_list)     : 0
    nangle_types    = nangles > 0    ? maximum(a[1] for a in angle_list)    : 0
    ndihedral_types = ndihedrals > 0 ? maximum(d[1] for d in dihedral_list) : 0

    data.bonds       = zeros(Int, nbonds, 2)
    data.bond_labels = zeros(Int, nbonds)
    for (k, (bt, i1, i2)) in enumerate(bond_list)
        data.bond_labels[k] = bt
        data.bonds[k, 1] = data.atom_ids[i1]
        data.bonds[k, 2] = data.atom_ids[i2]
    end
    data.nbond_types = nbond_types

    data.angles       = zeros(Int, nangles, 3)
    data.angle_labels = zeros(Int, nangles)
    for (k, (at, i1, ic, i2)) in enumerate(angle_list)
        data.angle_labels[k] = at
        data.angles[k, 1] = data.atom_ids[i1]
        data.angles[k, 2] = data.atom_ids[ic]
        data.angles[k, 3] = data.atom_ids[i2]
    end
    data.nangle_types = nangle_types

    data.dihedrals       = zeros(Int, ndihedrals, 4)
    data.dihedral_labels = zeros(Int, ndihedrals)
    for (k, (dt, i1, i2, i3, i4)) in enumerate(dihedral_list)
        data.dihedral_labels[k] = dt
        data.dihedrals[k, 1] = data.atom_ids[i1]
        data.dihedrals[k, 2] = data.atom_ids[i2]
        data.dihedrals[k, 3] = data.atom_ids[i3]
        data.dihedrals[k, 4] = data.atom_ids[i4]
    end
    data.ndihedral_types = ndihedral_types

    data.impropers       = zeros(Int, 0, 4)
    data.improper_labels = zeros(Int, 0)
    data.nimproper_types = 0

    if verbose
        bond_names     = Dict(1=>"Si-O", 2=>"Al-O", 3=>"O-H")
        angle_names    = Dict(1=>"Si-O-Si", 2=>"O-Si-O", 3=>"Si-O-Al",
                              4=>"O-Al-O",  5=>"Si-O-H", 6=>"Al-O-H")
        dihedral_names = Dict(1=>"O-Si-O-Si", 2=>"O-Si-O-Al",
                              3=>"O-Al-O-Si", 4=>"O-Al-O-Al",
                              5=>"O-Si-O-H",  6=>"O-Al-O-H")

        println("\n  === Topology summary ===")
        println("  Bonds:     $nbonds  ($nbond_types types)")
        for bt in sort(unique(data.bond_labels))
            n = count(==(bt), data.bond_labels)
            println("    type $bt  $(get(bond_names, bt, "?")): $n")
        end
        println("  Angles:    $nangles  ($nangle_types types)")
        for at in sort(unique(data.angle_labels))
            n = count(==(at), data.angle_labels)
            println("    type $at  $(get(angle_names, at, "?")): $n")
        end
        println("  Dihedrals: $ndihedrals  ($ndihedral_types types)")
        for dt in sort(unique(data.dihedral_labels))
            n = count(==(dt), data.dihedral_labels)
            println("    type $dt  $(get(dihedral_names, dt, "?")): $n")
        end
        println("  Impropers: 0")

        println("\n  === Coefficient map for LAMMPS input ===")
        println("  bond_coeff  1  ...  # Si-O")
        has_al && println("  bond_coeff  2  ...  # Al-O")
        has_h  && println("  bond_coeff  3  ...  # O-H")
        println("  angle_coeff 1  ...  # Si-O-Si")
        println("  angle_coeff 2  ...  # O-Si-O")
        if has_al
            println("  angle_coeff 3  ...  # Si-O-Al")
            println("  angle_coeff 4  ...  # O-Al-O")
        end
        if has_h
            println("  angle_coeff 5  ...  # Si-O-H")
            has_al && println("  angle_coeff 6  ...  # Al-O-H")
        end
        println("  dihedral_coeff 1 ... # O-Si-O-Si")
        if has_al
            println("  dihedral_coeff 2 ... # O-Si-O-Al")
            println("  dihedral_coeff 3 ... # O-Al-O-Si")
            println("  dihedral_coeff 4 ... # O-Al-O-Al")
        end
        if has_h
            println("  dihedral_coeff 5 ... # O-Si-O-H")
            has_al && println("  dihedral_coeff 6 ... # O-Al-O-H")
        end
    end

    return data
end


# ============================================================================
# Standalone CLI
# ============================================================================

function main(args=ARGS)
    if length(args) < 2
        println("""
        Usage: julia add_zeolite_topology.jl input.data output.data [options]

        Options:
          --si-type N        Si atom type (default: 2)
          --o-type N         O atom type  (default: 1)
          --al-type N        Al atom type (0 or omit to skip)
          --h-type N         Bronsted H type (0 or omit to skip)
          --si-o-cutoff R    Si-O cutoff in A (default: 1.85)
          --al-o-cutoff R    Al-O cutoff in A (default: 2.00)
          --o-h-cutoff R     O-H cutoff in A  (default: 1.05)

        Type map (Si/Al fully distinguished):
          Bonds:     1=Si-O  2=Al-O  3=O-H
          Angles:    1=Si-O-Si  2=O-Si-O  3=Si-O-Al  4=O-Al-O
                     5=Si-O-H   6=Al-O-H
          Dihedrals: 1=O-Si-O-Si  2=O-Si-O-Al  3=O-Al-O-Si
                     4=O-Al-O-Al  5=O-Si-O-H   6=O-Al-O-H
        """)
        return
    end

    input_file  = args[1]
    output_file = args[2]

    si_type = 2; o_type = 1; al_type = nothing; h_type = nothing
    si_o_cutoff = 1.85; al_o_cutoff = 2.00; o_h_cutoff = 1.05

    i = 3
    while i <= length(args)
        a = args[i]
        if a == "--si-type" && i < length(args)
            si_type = parse(Int, args[i+1]); i += 2
        elseif a == "--o-type" && i < length(args)
            o_type = parse(Int, args[i+1]); i += 2
        elseif a == "--al-type" && i < length(args)
            v = parse(Int, args[i+1])
            al_type = v > 0 ? v : nothing; i += 2
        elseif a == "--h-type" && i < length(args)
            v = parse(Int, args[i+1])
            h_type = v > 0 ? v : nothing; i += 2
        elseif a == "--si-o-cutoff" && i < length(args)
            si_o_cutoff = parse(Float64, args[i+1]); i += 2
        elseif a == "--al-o-cutoff" && i < length(args)
            al_o_cutoff = parse(Float64, args[i+1]); i += 2
        elseif a == "--o-h-cutoff" && i < length(args)
            o_h_cutoff = parse(Float64, args[i+1]); i += 2
        else
            @warn "Unknown argument: $a"; i += 1
        end
    end

    println("Reading $input_file...")
    data = LammpsDataReader.read_lammps_data(input_file; verbose=false)
    println("  $(size(data.coords,1)) atoms")

    println("Building zeolite topology...")
    add_zeolite_topology!(data;
        si_type=si_type, o_type=o_type,
        al_type=al_type, h_type=h_type,
        si_o_cutoff=si_o_cutoff, al_o_cutoff=al_o_cutoff,
        o_h_cutoff=o_h_cutoff)

    println("\nWriting $output_file...")
    LammpsDataReader.write_lammps_data(output_file, data;
        comment="LAMMPS data - zeolite topology (add_zeolite_topology.jl v3)")
    println("Done!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
