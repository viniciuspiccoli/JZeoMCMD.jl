# ════════════════════════════════════════════════════════════════
# analysis.jl — Strain tracking, LAMMPS log parsing, channel occupancy
# ════════════════════════════════════════════════════════════════

"""
    parse_lammps_log(filename) -> Dict{String, Vector{Float64}}

Parse a LAMMPS thermo log file and return columns as vectors.
Handles the standard thermo output format with column headers.
"""
function parse_lammps_log(filename::String)
    lines = readlines(filename)
    data = Dict{String, Vector{Float64}}()
    headers = String[]
    in_thermo = false

    for line in lines
        s = strip(line)

        # Detect header line (starts with "Step")
        if startswith(s, "Step")
            headers = split(s)
            for h in headers
                data[h] = Float64[]
            end
            in_thermo = true
            continue
        end

        # Detect end of thermo block
        if in_thermo && (startswith(s, "Loop") || startswith(s, "WARNING") ||
                         isempty(s) || startswith(s, "="))
            in_thermo = false
            continue
        end

        if in_thermo && !isempty(headers)
            tokens = split(s)
            length(tokens) == length(headers) || continue
            for (i, h) in enumerate(headers)
                val = tryparse(Float64, tokens[i])
                val !== nothing && push!(data[h], val)
            end
        end
    end

    return data
end

"""
    extract_cell_params(log_data) -> NamedTuple

Extract final cell parameters from parsed LAMMPS log.
Returns (a, b, c, alpha, beta, gamma, volume).
"""
function extract_cell_params(log_data::Dict)
    get_last(key) = haskey(log_data, key) && !isempty(log_data[key]) ?
                    log_data[key][end] : NaN

    return (
        a     = get_last("Cella"),
        b     = get_last("Cellb"),
        c     = get_last("Cellc"),
        alpha = get_last("CellAlpha"),
        beta  = get_last("CellBeta"),
        gamma = get_last("CellGamma"),
        volume = get_last("Volume"),
    )
end

"""
    compute_strain(cell_current, cell_reference) -> NamedTuple

Compute linear strain ε = (d - d₀)/d₀ for each cell parameter.
Returns (εa, εb, εc, εV) — volumetric strain is (V-V₀)/V₀.
"""
function compute_strain(current::NamedTuple, reference::NamedTuple)
    strain(x, x0) = isnan(x) || isnan(x0) || x0 == 0 ? NaN : (x - x0) / x0
    return (
        εa = strain(current.a, reference.a),
        εb = strain(current.b, reference.b),
        εc = strain(current.c, reference.c),
        εV = strain(current.volume, reference.volume),
    )
end

"""
    parse_raspa3_output(output_dir) -> Dict

Parse RASPA3 output for average loading, energy, etc.
Searches for the output data file in Output/System_0/.
"""
function parse_raspa3_output(output_dir::String)
    result = Dict{String, Any}(
        "n_ads" => NaN,
        "n_ads_err" => NaN,
        "energy" => NaN,
        "volume" => NaN,
    )

    # Search for output files
    out_dir = joinpath(output_dir, "Output", "System_0")
    !isdir(out_dir) && return result

    for fname in readdir(out_dir)
        endswith(fname, ".data") || continue
        content = read(joinpath(out_dir, fname), String)

        # Average loading [molecules/uc]
        m = match(r"Average loading absolute \[molecules/uc\]\s*:\s*([\d.eE+-]+)\s*\+/-\s*([\d.eE+-]+)", content)
        if m !== nothing
            result["n_ads"] = parse(Float64, m.captures[1])
            result["n_ads_err"] = parse(Float64, m.captures[2])
        end

        # Try mol/kg too
        if isnan(result["n_ads"])
            m = match(r"Average loading absolute \[mol/kg\]\s*:\s*([\d.eE+-]+)", content)
            m !== nothing && (result["n_ads"] = parse(Float64, m.captures[1]))
        end

        # Average volume
        m = match(r"Average Volume:\s*([\d.eE+-]+)", content)
        m !== nothing && (result["volume"] = parse(Float64, m.captures[1]))
    end

    return result
end

# ════════════════════════════════════════════════════════════════
# MFI Channel Occupancy Analysis
# ════════════════════════════════════════════════════════════════

"""
    analyze_channel_occupancy(coords, atom_types, ads_types, box_dims)

Classify adsorbate molecules into MFI channel types based on
their center-of-mass position within the unit cell.

MFI has three distinct adsorption sites:
  - Straight channels (along b): y ≈ 0.25, 0.75 in fractional coords
  - Sinusoidal channels (along a): z ≈ 0.0, 0.5
  - Intersections: where straight meets sinusoidal

Returns Dict with counts and molecule lists for each channel type.
"""
function analyze_channel_occupancy(coords::Matrix{Float64},
                                    atom_types::Vector{Int},
                                    ads_types::Set{Int},
                                    box_dims::Matrix{Float64};
                                    atoms_per_mol::Int = 4)
    # Identify adsorbate atom indices
    ads_indices = findall(t -> t in ads_types, atom_types)
    n_ads = length(ads_indices)
    n_mols = n_ads ÷ atoms_per_mol

    L = [box_dims[d,2] - box_dims[d,1] for d in 1:3]
    lo = [box_dims[d,1] for d in 1:3]

    # MFI unit cell dimensions (for 1×1×1)
    # Straight channel runs along b, centered at x/a ≈ 0.0, 0.5
    # Sinusoidal runs along a, oscillating in z
    # For 2×2×2 supercell, there are 8 unit cells

    straight = Int[]
    sinusoidal = Int[]
    intersection = Int[]

    for mi in 1:n_mols
        # COM of molecule
        start_idx = ads_indices[(mi-1)*atoms_per_mol + 1]
        com = zeros(3)
        for k in 0:(atoms_per_mol-1)
            idx = ads_indices[(mi-1)*atoms_per_mol + 1 + k]
            com .+= coords[idx, :]
        end
        com ./= atoms_per_mol

        # Convert to fractional coordinates within a single unit cell
        frac = [(com[d] - lo[d]) / L[d] for d in 1:3]
        # Map to [0,1) within one unit cell (for 2x2x2, multiply by 2 and mod 1)
        # Assuming 2x2x2 supercell:
        frac_uc = [mod(f * 2, 1.0) for f in frac]

        fx, fy, fz = frac_uc

        # Classification based on MFI topology:
        # Intersection: near (0.0±0.1, 0.25±0.1) or (0.5±0.1, 0.25±0.1)
        # in fractional unit cell coordinates
        near_x_channel = (fx < 0.15 || fx > 0.85 || abs(fx - 0.5) < 0.15)
        near_y_quarter = (abs(fy - 0.25) < 0.12 || abs(fy - 0.75) < 0.12)

        if near_x_channel && near_y_quarter
            push!(intersection, mi)
        elseif near_y_quarter
            push!(straight, mi)
        else
            push!(sinusoidal, mi)
        end
    end

    return Dict(
        "total" => n_mols,
        "straight" => length(straight),
        "sinusoidal" => length(sinusoidal),
        "intersection" => length(intersection),
        "straight_mols" => straight,
        "sinusoidal_mols" => sinusoidal,
        "intersection_mols" => intersection,
    )
end

"""
    write_cycle_summary(filename, cycle_data; append=true)

Append one line to a CSV tracking file.
`cycle_data` is a Dict with keys: iteration, pressure, n_ads, a, b, c, V, εV, ...
"""
function write_cycle_summary(filename::String, d::Dict; append::Bool=true)
    header = "iteration,pressure,n_ads,n_straight,n_sinusoidal,n_intersection," *
             "a,b,c,alpha,beta,gamma,volume,strain_a,strain_b,strain_c,strain_V,timestamp"

    mode = append && isfile(filename) ? "a" : "w"
    open(filename, mode) do io
        mode == "w" && println(io, header)
        @printf(io, "%d,%.1f,%.2f,%d,%d,%d,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.2f,%.6e,%.6e,%.6e,%.6e,%s\n",
                get(d, "iteration", 0),
                get(d, "pressure", 0.0),
                get(d, "n_ads", NaN),
                get(d, "n_straight", 0),
                get(d, "n_sinusoidal", 0),
                get(d, "n_intersection", 0),
                get(d, "a", NaN), get(d, "b", NaN), get(d, "c", NaN),
                get(d, "alpha", NaN), get(d, "beta", NaN), get(d, "gamma", NaN),
                get(d, "volume", NaN),
                get(d, "strain_a", NaN), get(d, "strain_b", NaN),
                get(d, "strain_c", NaN), get(d, "strain_V", NaN),
                Dates.format(now(), "yyyy-mm-dd_HH:MM:SS"))
    end
end
