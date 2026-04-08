#module LammpsDataReader

#using Printf
#using LinearAlgebra
#using Statistics

#export LammpsData, read_lammps_data, write_lammps_data
#export RaspaForceField, read_raspa3_forcefield, map_raspa3_forcefield!
#export make_supercell, box_to_matrix
#export compute_rdf, compute_msd, read_lammps_dump, parse_raspa3_isotherm

# ════════════════════════════════════════════════════════════════════
#  Core Data Structure  (features 4, 5, 7 integrated)
# ════════════════════════════════════════════════════════════════════

"""
    LammpsData

Stores all data from a LAMMPS data file.

Supports orthogonal and triclinic boxes (`tilt_factors`),
image flags (`image_flags`), and atom styles `:full`, `:molecular`,
`:charge`, `:atomic` (`atom_style`).
"""
mutable struct LammpsData
    # core
    masses::Dict{Int,Float64}
    coords::Matrix{Float64}          # Natoms × 3
    molecule_labels::Vector{Int}
    atom_charges::Vector{Float64}
    atom_labels::Vector{Int}
    atom_ids::Vector{Int}
    box_dimensions::Matrix{Float64}  # 3×2  [xlo xhi; ylo yhi; zlo zhi]
    # (4) triclinic
    tilt_factors::Vector{Float64}    # [xy, xz, yz]  — all zero if orthogonal
    # (5) image flags
    image_flags::Matrix{Int}         # Natoms × 3  (ix iy iz)
    # topology
    bonds::Matrix{Int}
    bond_labels::Vector{Int}
    nbond_types::Int
    angles::Matrix{Int}
    angle_labels::Vector{Int}
    nangle_types::Int
    dihedrals::Matrix{Int}
    dihedral_labels::Vector{Int}
    ndihedral_types::Int
    impropers::Matrix{Int}
    improper_labels::Vector{Int}
    nimproper_types::Int
    # coefficients
    pair_coeffs::Dict{Int,Vector{Float64}}
    bond_coeffs::Dict{Int,Vector{Float64}}
    angle_coeffs::Dict{Int,Vector{Float64}}
    dihedral_coeffs::Dict{Int,Vector{Float64}}
    improper_coeffs::Dict{Int,Vector{Float64}}
    # (7) atom style tag
    atom_style::Symbol
end

function LammpsData()
    LammpsData(
        Dict{Int,Float64}(),
        zeros(0, 3), zeros(Int, 0), zeros(0), zeros(Int, 0), zeros(Int, 0),
        zeros(3, 2),
        [0.0, 0.0, 0.0],        # tilt
        zeros(Int, 0, 3),        # image flags
        zeros(Int, 0, 2), zeros(Int, 0), 0,
        zeros(Int, 0, 3), zeros(Int, 0), 0,
        zeros(Int, 0, 4), zeros(Int, 0), 0,
        zeros(Int, 0, 4), zeros(Int, 0), 0,
        Dict{Int,Vector{Float64}}(),
        Dict{Int,Vector{Float64}}(),
        Dict{Int,Vector{Float64}}(),
        Dict{Int,Vector{Float64}}(),
        Dict{Int,Vector{Float64}}(),
        :full,
    )
end

"""
    box_to_matrix(data::LammpsData) -> Matrix{Float64}

Return the 3×3 cell matrix (rows = lattice vectors) accounting for tilt.

    | lx   0   0 |
    | xy  ly   0 |
    | xz  yz  lz |
"""
function box_to_matrix(data::LammpsData)
    lx = data.box_dimensions[1,2] - data.box_dimensions[1,1]
    ly = data.box_dimensions[2,2] - data.box_dimensions[2,1]
    lz = data.box_dimensions[3,2] - data.box_dimensions[3,1]
    xy, xz, yz = data.tilt_factors
    return [lx  0.0 0.0;
            xy  ly  0.0;
            xz  yz  lz]
end


# ════════════════════════════════════════════════════════════════════
#  Internal parsing helpers
# ════════════════════════════════════════════════════════════════════

mutable struct _ParseState
    natoms::Int; nbonds::Int; nangles::Int; ndihedrals::Int; nimpropers::Int
    natom_types::Int; nbond_types::Int; nangle_types::Int
    ndihedral_types::Int; nimproper_types::Int
    xlo::Float64; xhi::Float64
    ylo::Float64; yhi::Float64
    zlo::Float64; zhi::Float64
    xy::Float64; xz::Float64; yz::Float64
    is_triclinic::Bool
    detected_style::Symbol
end
_ParseState() = _ParseState(
    0,0,0,0,0, 0,0,0,0,0,
    0.0,0.0, 0.0,0.0, 0.0,0.0,
    0.0, 0.0, 0.0,
    false, :full,
)

_parse_line(lines, i) = (split(strip(lines[i])), i + 1)


# ════════════════════════════════════════════════════════════════════
#  (7) Atom-style-aware atom reader
# ════════════════════════════════════════════════════════════════════

# Minimum data columns (excluding trailing image flags):
#   :full       → id mol type charge x y z           (7)
#   :molecular  → id mol type x y z                  (6)
#   :charge     → id type charge x y z               (6)
#   :atomic     → id type x y z                      (5)
const _STYLE_MIN_COLS = Dict(:full=>7, :molecular=>6, :charge=>6, :atomic=>5)

"""
Auto-detect atom style from column count of the first atom line.
"""
function _detect_atom_style(tok, natom_types::Int)
    n = length(tok)
    base = n >= 10 ? n - 3 : n   # strip possible image flags
    base >= 7 && return :full
    base == 6 && return :molecular   # heuristic default for 6-col
    base == 5 && return :atomic
    return :full
end

function _read_atoms_styled!(data, st, lines, i)
    n = st.natoms
    first_tok = split(strip(lines[i]))
    style = _detect_atom_style(first_tok, st.natom_types)
    st.detected_style = style

    data.atom_ids        = zeros(Int, n)
    data.coords          = zeros(n, 3)
    data.molecule_labels = zeros(Int, n)
    data.atom_labels     = zeros(Int, n)
    data.atom_charges    = zeros(n)
    data.image_flags     = zeros(Int, n, 3)

    min_cols = _STYLE_MIN_COLS[style]

    for j in 1:n
        tok, i = _parse_line(lines, i)
        ntok = length(tok)
        data.atom_ids[j] = parse(Int, tok[1])

        if style == :full
            data.molecule_labels[j] = parse(Int, tok[2])
            data.atom_labels[j]     = parse(Int, tok[3])
            data.atom_charges[j]    = parse(Float64, tok[4])
            data.coords[j, :]      .= parse.(Float64, tok[5:7])
            ntok >= 10 && (data.image_flags[j,:] .= parse.(Int, tok[8:10]))

        elseif style == :molecular
            data.molecule_labels[j] = parse(Int, tok[2])
            data.atom_labels[j]     = parse(Int, tok[3])
            data.coords[j, :]      .= parse.(Float64, tok[4:6])
            ntok >= 9 && (data.image_flags[j,:] .= parse.(Int, tok[7:9]))

        elseif style == :charge
            data.atom_labels[j]  = parse(Int, tok[2])
            data.atom_charges[j] = parse(Float64, tok[3])
            data.coords[j, :]   .= parse.(Float64, tok[4:6])
            ntok >= 9 && (data.image_flags[j,:] .= parse.(Int, tok[7:9]))

        elseif style == :atomic
            data.atom_labels[j] = parse(Int, tok[2])
            data.coords[j, :]  .= parse.(Float64, tok[3:5])
            ntok >= 8 && (data.image_flags[j,:] .= parse.(Int, tok[6:8]))
        end
    end
    return i
end


# ════════════════════════════════════════════════════════════════════
#  Main reader
# ════════════════════════════════════════════════════════════════════

"""
    read_lammps_data(filename; verbose=true, atom_style=:auto)

Read a LAMMPS data file.  Detects atom style from column count when
`atom_style=:auto`.  Handles triclinic tilt factors and image flags.
"""
function read_lammps_data(filename::AbstractString;
                          verbose::Bool=true,
                          atom_style::Symbol=:auto)
    data = LammpsData()
    st = _ParseState()
    lines = readlines(filename)
    i = 1
    while i <= length(lines)
        tokens = split(strip(lines[i]))
        i = _analyse!(data, st, tokens, lines, i, verbose)
    end
    data.box_dimensions  = [st.xlo st.xhi; st.ylo st.yhi; st.zlo st.zhi]
    st.is_triclinic && (data.tilt_factors = [st.xy, st.xz, st.yz])
    data.nbond_types     = st.nbond_types
    data.nangle_types    = st.nangle_types
    data.ndihedral_types = st.ndihedral_types
    data.nimproper_types = st.nimproper_types
    data.atom_style      = atom_style == :auto ? st.detected_style : atom_style
    verbose && _validate(data)
    return data
end

function _analyse!(data, st, tokens, lines, i, verbose)
    isempty(tokens) && return i + 1
    first = string(tokens[1])
    startswith(first, "#") && return i + 1

    header_numbers = Dict("atoms"=>:natoms,"bonds"=>:nbonds,"angles"=>:nangles,
                          "dihedrals"=>:ndihedrals,"impropers"=>:nimpropers)
    header_types   = Dict("atom"=>:natom_types,"bond"=>:nbond_types,"angle"=>:nangle_types,
                          "dihedral"=>:ndihedral_types,"improper"=>:nimproper_types)
    box_keys = [["xlo","xhi"],["ylo","yhi"],["zlo","zhi"]]
    ncoeff   = Dict("Pair"=>:natom_types,"Bond"=>:nbond_types,"Angle"=>:nangle_types,
                    "Dihedral"=>:ndihedral_types,"Improper"=>:nimproper_types)

    if tryparse(Float64, first) !== nothing
        n = length(tokens)
        if n == 2
            f = get(header_numbers, string(tokens[2]), nothing)
            f !== nothing ? setproperty!(st, f, parse(Int, first)) :
                (verbose && @warn "Unknown line $i")
        elseif n == 3
            f = get(header_types, string(tokens[2]), nothing)
            f !== nothing ? setproperty!(st, f, parse(Int, first)) :
                (verbose && @warn "Unknown line $i")
        elseif n == 4
            pair = [string(tokens[3]), string(tokens[4])]
            if pair in box_keys
                setproperty!(st, Symbol(pair[1]), parse(Float64, string(tokens[1])))
                setproperty!(st, Symbol(pair[2]), parse(Float64, string(tokens[2])))
            elseif verbose
                @warn "Unknown line $i"
            end
        # ── (4) triclinic tilt: val val val xy xz yz ──
        elseif n == 6 && string(tokens[4])=="xy" && string(tokens[5])=="xz" && string(tokens[6])=="yz"
            st.xy = parse(Float64, string(tokens[1]))
            st.xz = parse(Float64, string(tokens[2]))
            st.yz = parse(Float64, string(tokens[3]))
            st.is_triclinic = true
        end
        return i + 1

    elseif length(tokens) == 2 && haskey(ncoeff, first)
        verbose && println("reading $first Coeffs")
        N = getproperty(st, ncoeff[first])
        name = lowercase(first) * "_" * lowercase(string(tokens[2]))
        i += 2
        return _read_coeffs!(data, name, N, lines, i)

    else
        sections = Dict(
            "Masses"     => (d,s,ls,idx)->_read_masses!(d,s,ls,idx),
            "Atoms"      => (d,s,ls,idx)->_read_atoms_styled!(d,s,ls,idx),
            "Velocities" => (d,s,ls,idx)->(_skip_lines(s.natoms,ls,idx)),
            "Bonds"      => (d,s,ls,idx)->_read_topo!(d,:bond_labels,:bonds,s.nbonds,2,ls,idx),
            "Angles"     => (d,s,ls,idx)->_read_topo!(d,:angle_labels,:angles,s.nangles,3,ls,idx),
            "Dihedrals"  => (d,s,ls,idx)->_read_topo!(d,:dihedral_labels,:dihedrals,s.ndihedrals,4,ls,idx),
            "Impropers"  => (d,s,ls,idx)->_read_topo!(d,:improper_labels,:impropers,s.nimpropers,4,ls,idx),
        )
        func = get(sections, first, nothing)
        if func !== nothing
            verbose && println("reading $first")
            i += 2
            return func(data, st, lines, i)
        else
            verbose && @warn "Unknown line $i: $(join(tokens,' '))"
            return i + 1
        end
    end
end

# ── section readers ─────────────────────────────────────────────────

function _read_masses!(data, st, lines, i)
    data.masses = Dict{Int,Float64}()
    for _ in 1:st.natom_types
        tok, i = _parse_line(lines, i)
        data.masses[parse(Int, tok[1])] = parse(Float64, tok[2])
    end
    return i
end

function _skip_lines(n, lines, i)
    for _ in 1:n; _, i = _parse_line(lines, i); end
    return i
end

function _read_topo!(data, label_field::Symbol, data_field::Symbol, n, ncols, lines, i)
    labels = zeros(Int, n)
    atoms  = zeros(Int, n, ncols)
    for _ in 1:n
        tok, i = _parse_line(lines, i)
        idx = parse(Int, tok[1])
        labels[idx] = parse(Int, tok[2])
        atoms[idx, :] .= parse.(Int, tok[3:2+ncols])
    end
    setfield!(data, label_field, labels)
    setfield!(data, data_field, atoms)
    return i
end

function _read_coeffs!(data, name, n, lines, i)
    coeffs = Dict{Int,Vector{Float64}}()
    for _ in 1:n
        tok, i = _parse_line(lines, i)
        coeffs[parse(Int, tok[1])] = parse.(Float64, tok[2:end])
    end
    sym = Symbol(name)
    hasfield(LammpsData, sym) && setfield!(data, sym, coeffs)
    return i
end

function _validate(data::LammpsData)
    skip = (:masses,:pair_coeffs,:bond_coeffs,:angle_coeffs,:dihedral_coeffs,:improper_coeffs,:atom_style,:tilt_factors)
    for field in fieldnames(LammpsData)
        field in skip && continue
        val = getfield(data, field)
        (val isa AbstractArray && isempty(val)) && @warn "Undefined or empty: $field"
    end
end


# ════════════════════════════════════════════════════════════════════
#  (3) Force-field mapper:  RASPA3 JSON → LAMMPS
# ════════════════════════════════════════════════════════════════════

"""
    RaspaForceField

Parsed contents of RASPA3 `force_field.json` and `pseudo_atoms.json`.
"""
struct RaspaForceField
    lj_params::Dict{String, Tuple{Float64,Float64}}          # name → (ε/K, σ/Å)
    pseudo_atoms::Dict{String, @NamedTuple{mass::Float64, charge::Float64, element::String}}
end

"""
    read_raspa3_forcefield(ff_json, pseudo_json) -> RaspaForceField

Read RASPA3 force-field definition files.

Expected `force_field.json` layout (either array or dict under key
`"LennardJones"` / `"lennard_jones"`):

```json
{ "LennardJones": [
    { "Name": "CH4_sp3", "epsilon": 148.0, "sigma": 3.73 }, ...
  ]
}
```

Expected `pseudo_atoms.json`:

```json
{ "PseudoAtoms": [
    { "Name": "CH4_sp3", "mass": 16.04, "charge": 0.0, "element": "C" }, ...
  ]
}
```

Also accepts top-level arrays in either file.
"""
function read_raspa3_forcefield(ff_json::AbstractString, pseudo_json::AbstractString)
    ff_raw = _load_json(ff_json)
    pa_raw = _load_json(pseudo_json)

    # ── LJ parameters ──
    lj = Dict{String, Tuple{Float64,Float64}}()
    lj_entries = _get_any_key(ff_raw, ["LennardJones","lennard_jones","Lennard-Jones"])
    if lj_entries isa AbstractVector
        for e in lj_entries
            lj[e["Name"]] = (Float64(e["epsilon"]), Float64(e["sigma"]))
        end
    elseif lj_entries isa AbstractDict
        for (name, p) in lj_entries
            lj[name] = (Float64(p["epsilon"]), Float64(p["sigma"]))
        end
    end

    # ── Pseudo atoms ──
    pa = Dict{String, @NamedTuple{mass::Float64, charge::Float64, element::String}}()
    pa_list = pa_raw isa AbstractVector ? pa_raw :
              _get_any_key(pa_raw, ["PseudoAtoms","pseudo_atoms","Pseudo_Atoms"])
    if pa_list !== nothing
        for e in pa_list
            pa[e["Name"]] = (
                mass    = Float64(get(e, "mass", 0.0)),
                charge  = Float64(get(e, "charge", 0.0)),
                element = string(get(e, "element", "X")),
            )
        end
    end
    return RaspaForceField(lj, pa)
end

# Minimal JSON reader (avoids hard dependency on JSON.jl).
# Falls back to shell `python3 -c 'import json; ...'` if JSON.jl
# is not installed; but you really should just `Pkg.add("JSON")`.
function _load_json(path::AbstractString)
    try
        # try JSON.jl first
        JSON_mod = Base.require(Base.PkgId(Base.UUID("682c06a0-de6a-54ab-a142-c8b1cf79cde6"), "JSON"))
        return JSON_mod.parsefile(path)
    catch
        # fallback: read and eval (only safe for trusted files!)
        raw = read(path, String)
        # very naive JSON → Julia dict via replace
        raw = replace(raw, "true"=>"true", "false"=>"false", "null"=>"nothing")
        @warn "JSON.jl not found — using naive parser.  `Pkg.add(\"JSON\")` is recommended."
        return Meta.eval(Meta.parse(raw))
    end
end

function _get_any_key(d, keys)
    d isa AbstractDict || return d
    for k in keys
        haskey(d, k) && return d[k]
    end
    return nothing
end

"""
    map_raspa3_forcefield!(data, raspa_ff, type_map; energy_units=:kcalmol)

Apply RASPA3 LJ parameters to `data.pair_coeffs` and masses.

`type_map` maps LAMMPS atom-type IDs → RASPA pseudo-atom names,
e.g. `Dict(1 => "Ow", 2 => "Hw")`.

RASPA3 ε is in Kelvin; set `energy_units=:kcalmol` (default) to
convert via kB = 0.001987 kcal/(mol·K), or `:kelvin` to keep as-is.
"""
function map_raspa3_forcefield!(data::LammpsData,
                                 raspa_ff::RaspaForceField,
                                 type_map::Dict{Int,String};
                                 energy_units::Symbol=:kcalmol)
    kB = 0.0019872041  # kcal/(mol·K)
    conv = energy_units == :kcalmol ? kB : 1.0

    data.pair_coeffs = Dict{Int,Vector{Float64}}()
    for (ltype, rname) in type_map
        if !haskey(raspa_ff.lj_params, rname)
            @warn "No LJ parameters for pseudo-atom '$rname'"
            continue
        end
        ε_K, σ = raspa_ff.lj_params[rname]
        data.pair_coeffs[ltype] = [ε_K * conv, σ]

        if haskey(raspa_ff.pseudo_atoms, rname)
            data.masses[ltype] = raspa_ff.pseudo_atoms[rname].mass
        end
    end
    return data
end


# ════════════════════════════════════════════════════════════════════
#  (6) Supercell builder
# ════════════════════════════════════════════════════════════════════

"""
    make_supercell(data::LammpsData, nx, ny, nz) -> LammpsData

Replicate the unit cell `nx × ny × nz` times.  All atom IDs,
topology indices, and molecule labels are correctly remapped.
Works for both orthogonal and triclinic cells.
"""
function make_supercell(data::LammpsData, nx::Int, ny::Int, nz::Int)
    cell = box_to_matrix(data)
    a, b, c = cell[1,:], cell[2,:], cell[3,:]

    natoms_old = size(data.coords, 1)
    nrep = nx * ny * nz
    natoms_new = natoms_old * nrep

    sc = LammpsData()
    sc.atom_style      = data.atom_style
    sc.masses          = copy(data.masses)
    sc.pair_coeffs     = deepcopy(data.pair_coeffs)
    sc.bond_coeffs     = deepcopy(data.bond_coeffs)
    sc.angle_coeffs    = deepcopy(data.angle_coeffs)
    sc.dihedral_coeffs = deepcopy(data.dihedral_coeffs)
    sc.improper_coeffs = deepcopy(data.improper_coeffs)
    sc.nbond_types     = data.nbond_types
    sc.nangle_types    = data.nangle_types
    sc.ndihedral_types = data.ndihedral_types
    sc.nimproper_types = data.nimproper_types

    # expanded box
    bd = data.box_dimensions
    sc.box_dimensions = [
        bd[1,1]  bd[1,1]+nx*(bd[1,2]-bd[1,1]);
        bd[2,1]  bd[2,1]+ny*(bd[2,2]-bd[2,1]);
        bd[3,1]  bd[3,1]+nz*(bd[3,2]-bd[3,1])
    ]
    sc.tilt_factors = data.tilt_factors .* [nx, nx, ny]  # xy→nx, xz→nx, yz→ny

    # allocate atoms
    sc.coords          = zeros(natoms_new, 3)
    sc.atom_ids        = zeros(Int, natoms_new)
    sc.atom_labels     = zeros(Int, natoms_new)
    sc.atom_charges    = zeros(natoms_new)
    sc.molecule_labels = zeros(Int, natoms_new)
    sc.image_flags     = zeros(Int, natoms_new, 3)

    max_mol = maximum(data.molecule_labels; init=0)
    idx = 0; rep = 0
    for ix in 0:nx-1, iy in 0:ny-1, iz in 0:nz-1
        shift = ix .* a .+ iy .* b .+ iz .* c
        for j in 1:natoms_old
            idx += 1
            sc.coords[idx, :]       = data.coords[j, :] .+ shift
            sc.atom_ids[idx]        = idx
            sc.atom_labels[idx]     = data.atom_labels[j]
            sc.atom_charges[idx]    = data.atom_charges[j]
            sc.molecule_labels[idx] = data.molecule_labels[j] + rep * max_mol
        end
        rep += 1
    end

    # replicate topology
    function _rep_topo(old_data, old_labels, ncols)
        nold = size(old_data, 1)
        nold == 0 && return (zeros(Int, 0, ncols), zeros(Int, 0))
        new_data   = zeros(Int, nold * nrep, ncols)
        new_labels = zeros(Int, nold * nrep)
        for r in 0:nrep-1
            ao = r * natoms_old
            to = r * nold
            for k in 1:nold
                new_labels[to+k]    = old_labels[k]
                new_data[to+k, :]  .= old_data[k, :] .+ ao
            end
        end
        return new_data, new_labels
    end

    sc.bonds,     sc.bond_labels     = _rep_topo(data.bonds,     data.bond_labels,     2)
    sc.angles,    sc.angle_labels    = _rep_topo(data.angles,    data.angle_labels,    3)
    sc.dihedrals, sc.dihedral_labels = _rep_topo(data.dihedrals, data.dihedral_labels, 4)
    sc.impropers, sc.improper_labels = _rep_topo(data.impropers, data.improper_labels, 4)

    return sc
end


# ════════════════════════════════════════════════════════════════════
#  (8) Post-processing helpers
# ════════════════════════════════════════════════════════════════════

# ── 8a  Radial Distribution Function ───────────────────────────────

"""
    compute_rdf(coords, box_dims, type_a, type_b, atom_labels;
                nbins=200, rmax=nothing) -> (r, g_r)

Pair RDF g(r) between atom types `type_a` and `type_b` under the
minimum-image convention (orthogonal box).

Returns bin centres `r` and normalised `g_r`.
"""
function compute_rdf(coords::Matrix{Float64},
                     box_dims::Matrix{Float64},
                     type_a::Int, type_b::Int,
                     atom_labels::Vector{Int};
                     nbins::Int=200,
                     rmax::Union{Nothing,Float64}=nothing)
    L = box_dims[:,2] .- box_dims[:,1]
    halfL = minimum(L) / 2
    rmax = something(rmax, halfL)
    dr = rmax / nbins

    idxA = findall(==(type_a), atom_labels)
    idxB = findall(==(type_b), atom_labels)
    same = type_a == type_b

    hist = zeros(nbins)
    for ia in idxA
        ra = @view coords[ia, :]
        for ib in idxB
            same && ib <= ia && continue
            dx = ra .- @view(coords[ib, :])
            dx .-= L .* round.(dx ./ L)
            r = sqrt(sum(dx .^ 2))
            if r < rmax
                bin = clamp(ceil(Int, r / dr), 1, nbins)
                hist[bin] += 1.0
            end
        end
    end

    vol = prod(L)
    ρB = length(idxB) / vol
    nA = length(idxA)
    r_centers = collect(range(dr/2, step=dr, length=nbins))
    g_r = similar(r_centers)
    for k in 1:nbins
        shell = (4π/3) * ((k*dr)^3 - ((k-1)*dr)^3)
        ideal = nA * ρB * shell
        g_r[k] = same ? hist[k] * 2 / ideal : hist[k] / ideal
    end
    return r_centers, g_r
end


# ── 8b  Mean Squared Displacement ──────────────────────────────────

"""
    compute_msd(trajectory; skip=1) -> (frame_indices, msd)

MSD from a vector of `Natoms × 3` coordinate matrices (unwrapped).
"""
function compute_msd(trajectory::Vector{Matrix{Float64}}; skip::Int=1)
    ref = trajectory[1]
    frames = 1:skip:length(trajectory)
    msd = zeros(length(frames))
    for (k, f) in enumerate(frames)
        disp = trajectory[f] .- ref
        msd[k] = mean(sum(disp .^ 2, dims=2))
    end
    return collect(frames), msd
end

"""
    read_lammps_dump(filename; skip_every=1) -> Vector{Matrix{Float64}}

Minimal reader for LAMMPS custom dump (`id xu yu zu` columns).
Returns one coordinate matrix per retained frame.
"""
function read_lammps_dump(filename::AbstractString; skip_every::Int=1)
    trajectory = Matrix{Float64}[]
    lines = readlines(filename)
    i = 1; frame = 0
    while i <= length(lines)
        if startswith(lines[i], "ITEM: NUMBER OF ATOMS")
            i += 1
            natoms = parse(Int, strip(lines[i])); i += 1
            while i <= length(lines) && !startswith(lines[i], "ITEM: ATOMS")
                i += 1
            end
            i += 1  # skip header
            frame += 1
            coords = zeros(natoms, 3)
            for _ in 1:natoms
                tok = split(strip(lines[i]))
                aid = parse(Int, tok[1])
                coords[aid, :] .= parse.(Float64, tok[2:4])
                i += 1
            end
            frame % skip_every == 0 && push!(trajectory, coords)
        else
            i += 1
        end
    end
    return trajectory
end


# ── 8c  RASPA3 isotherm parser ─────────────────────────────────────

"""
    parse_raspa3_isotherm(output_dir; component="methane")
        -> (pressures, loadings, errors)

Walk RASPA3 output sub-directories (one per pressure point) and
extract the adsorption isotherm.  Looks for lines matching:

    Average loading absolute [mol/kg]:  <value> +/- <error>
    Average loading absolute [molecules/uc]:  <value> +/- <error>

Returns vectors sorted by ascending pressure.
"""
function parse_raspa3_isotherm(output_dir::AbstractString;
                                component::String="methane")
    pressures = Float64[]; loadings = Float64[]; errors = Float64[]

    for entry in sort(readdir(output_dir))
        path = joinpath(output_dir, entry)
        isdir(path) || continue

        for fname in readdir(path)
            (endswith(fname, ".data") || fname == "output.txt") || continue
            content = read(joinpath(path, fname), String)

            # pressure
            mp = match(r"(?i)pressure\s*[:=]\s*([\d.eE+-]+)", content)
            p = mp !== nothing ? parse(Float64, mp.captures[1]) : tryparse(Float64, entry)
            p === nothing && continue

            # loading (try mol/kg first, then molecules/uc)
            ml = match(r"Average loading absolute \[mol/kg\]\s*:\s*([\d.eE+-]+)\s*\+/-\s*([\d.eE+-]+)", content)
            if ml === nothing
                ml = match(r"Average loading absolute \[molecules/uc\]\s*:\s*([\d.eE+-]+)\s*\+/-\s*([\d.eE+-]+)", content)
            end
            ml === nothing && continue

            push!(pressures, p)
            push!(loadings, parse(Float64, ml.captures[1]))
            push!(errors,   parse(Float64, ml.captures[2]))
        end
    end

    perm = sortperm(pressures)
    return pressures[perm], loadings[perm], errors[perm]
end


# ════════════════════════════════════════════════════════════════════
#  Writer  (features 4, 5, 7 aware)
# ════════════════════════════════════════════════════════════════════

"""
    write_lammps_data(filename, data; comment="...")

Write a LAMMPS data file respecting the stored atom style,
triclinic tilt factors, and image flags.
"""
function write_lammps_data(filename::AbstractString, data::LammpsData;
                           comment::String="LAMMPS data via LammpsDataReader.jl")
    natoms     = size(data.coords, 1)
    nbonds     = size(data.bonds, 1)
    nangles    = size(data.angles, 1)
    ndihedrals = size(data.dihedrals, 1)
    nimpropers = size(data.impropers, 1)
    natom_types = isempty(data.masses) ? maximum(data.atom_labels; init=0) :
                                         maximum(keys(data.masses))
    has_img = !isempty(data.image_flags) && any(data.image_flags .!= 0)
    is_tri  = any(data.tilt_factors .!= 0.0)

    open(filename, "w") do io
        println(io, comment, "\n")
        println(io, "$natoms atoms")
        nbonds > 0     && println(io, "$nbonds bonds")
        nangles > 0    && println(io, "$nangles angles")
        ndihedrals > 0 && println(io, "$ndihedrals dihedrals")
        nimpropers > 0 && println(io, "$nimpropers impropers")
        println(io)
        println(io, "$natom_types atom types")
        data.nbond_types > 0     && println(io, "$(data.nbond_types) bond types")
        data.nangle_types > 0    && println(io, "$(data.nangle_types) angle types")
        data.ndihedral_types > 0 && println(io, "$(data.ndihedral_types) dihedral types")
        data.nimproper_types > 0 && println(io, "$(data.nimproper_types) improper types")
        println(io)
        @printf(io, "%.10f %.10f xlo xhi\n", data.box_dimensions[1,1], data.box_dimensions[1,2])
        @printf(io, "%.10f %.10f ylo yhi\n", data.box_dimensions[2,1], data.box_dimensions[2,2])
        @printf(io, "%.10f %.10f zlo zhi\n", data.box_dimensions[3,1], data.box_dimensions[3,2])
        is_tri && @printf(io, "%.10f %.10f %.10f xy xz yz\n",
                          data.tilt_factors[1], data.tilt_factors[2], data.tilt_factors[3])

        # Masses
        if !isempty(data.masses)
            println(io, "\nMasses\n")
            for t in sort(collect(keys(data.masses)))
                @printf(io, "%d %.6f\n", t, data.masses[t])
            end
        end

        # Pair Coeffs
        if !isempty(data.pair_coeffs)
            println(io, "\nPair Coeffs\n")
            for t in sort(collect(keys(data.pair_coeffs)))
                print(io, "$t")
                for v in data.pair_coeffs[t]; @printf(io, " %.8f", v); end
                println(io)
            end
        end

        # Atoms
        style = data.atom_style
        println(io, "\nAtoms # $style\n")
        for j in 1:natoms
            id = data.atom_ids[j]
            if style == :full
                @printf(io, "%d %d %d %.8f %.10f %.10f %.10f",
                        id, data.molecule_labels[j], data.atom_labels[j],
                        data.atom_charges[j],
                        data.coords[j,1], data.coords[j,2], data.coords[j,3])
            elseif style == :molecular
                @printf(io, "%d %d %d %.10f %.10f %.10f",
                        id, data.molecule_labels[j], data.atom_labels[j],
                        data.coords[j,1], data.coords[j,2], data.coords[j,3])
            elseif style == :charge
                @printf(io, "%d %d %.8f %.10f %.10f %.10f",
                        id, data.atom_labels[j], data.atom_charges[j],
                        data.coords[j,1], data.coords[j,2], data.coords[j,3])
            elseif style == :atomic
                @printf(io, "%d %d %.10f %.10f %.10f",
                        id, data.atom_labels[j],
                        data.coords[j,1], data.coords[j,2], data.coords[j,3])
            end
            has_img && @printf(io, " %d %d %d",
                               data.image_flags[j,1], data.image_flags[j,2], data.image_flags[j,3])
            println(io)
        end

        # Topology
        function _write_sec(io, name, labels, atoms)
            size(atoms,1) == 0 && return
            println(io, "\n$name\n")
            for k in 1:size(atoms,1)
                print(io, "$k $(labels[k])")
                for c in 1:size(atoms,2); print(io, " $(atoms[k,c])"); end
                println(io)
            end
        end
        _write_sec(io, "Bonds",     data.bond_labels,     data.bonds)
        _write_sec(io, "Angles",    data.angle_labels,    data.angles)
        _write_sec(io, "Dihedrals", data.dihedral_labels, data.dihedrals)
        _write_sec(io, "Impropers", data.improper_labels, data.impropers)
    end
    return nothing
end

#end # module
