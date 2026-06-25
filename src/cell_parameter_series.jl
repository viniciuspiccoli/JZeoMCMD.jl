# Robust parser for cell_params.dat-style production time series.
#
# The parser is independent of the statistical calculations. It converts a
# whitespace-delimited file into a canonical, typed representation that later
# steps can pass to block averaging and drift analysis.

const _CELL_PARAMETER_CORE_COLUMNS = (
    :a,
    :b,
    :c,
    :alpha,
    :beta,
    :gamma,
    :volume,
)

const _CELL_PARAMETER_OPTIONAL_COLUMNS = (
    :temperature,
    :pressure,
    :potential_energy,
    :total_energy,
    :enthalpy,
    :density,
)

const _CELL_PARAMETER_ALIASES = Dict{String,Symbol}(
    # independent variables
    "step" => :step,
    "timestep" => :step,
    "time_step" => :step,
    "time" => :time,

    # cell lengths
    "a" => :a,
    "cella" => :a,
    "cell_a" => :a,
    "lx" => :a,
    "xlength" => :a,
    "x_length" => :a,

    "b" => :b,
    "cellb" => :b,
    "cell_b" => :b,
    "ly" => :b,
    "ylength" => :b,
    "y_length" => :b,

    "c" => :c,
    "cellc" => :c,
    "cell_c" => :c,
    "lz" => :c,
    "zlength" => :c,
    "z_length" => :c,

    # cell angles
    "alpha" => :alpha,
    "cellalpha" => :alpha,
    "cell_alpha" => :alpha,

    "beta" => :beta,
    "cellbeta" => :beta,
    "cell_beta" => :beta,

    "gamma" => :gamma,
    "cellgamma" => :gamma,
    "cell_gamma" => :gamma,

    # thermodynamic observables
    "volume" => :volume,
    "vol" => :volume,
    "cellvolume" => :volume,
    "cell_volume" => :volume,

    "temp" => :temperature,
    "temperature" => :temperature,

    "press" => :pressure,
    "pressure" => :pressure,

    "pe" => :potential_energy,
    "poteng" => :potential_energy,
    "potentialenergy" => :potential_energy,
    "potential_energy" => :potential_energy,

    "etotal" => :total_energy,
    "toteng" => :total_energy,
    "totalenergy" => :total_energy,
    "total_energy" => :total_energy,

    "enthalpy" => :enthalpy,
    "density" => :density,

    # triclinic box terms, retained for future analyses
    "xy" => :xy,
    "xz" => :xz,
    "yz" => :yz,
)

"""
    CellParameterParseError

Error raised when a cell-parameter time-series file cannot be interpreted
safely. `line == 0` denotes a file-level error rather than a particular row.
"""
struct CellParameterParseError <: Exception
    path::String
    line::Int
    message::String
end

function Base.showerror(io::IO, err::CellParameterParseError)
    if err.line > 0
        print(io, "failed to parse ", err.path, " at line ", err.line, ": ", err.message)
    else
        print(io, "failed to parse ", err.path, ": ", err.message)
    end
end

"""
    CellParameterSeries

Canonical representation of a whitespace-delimited production time series.
Rows are samples and columns are observables. Column names are normalized to
symbols such as `:a`, `:volume`, `:temperature`, and `:step`.

`header_source` is one of:

- `:explicit`: names were supplied through the `columns` keyword;
- `:file`: names were read from a commented or plain-text header;
- `:inferred`: a supported legacy layout was inferred from the data.
"""
struct CellParameterSeries
    source::String
    columns::Vector{Symbol}
    values::Matrix{Float64}
    line_numbers::Vector{Int}
    header_source::Symbol

    function CellParameterSeries(
        source::AbstractString,
        columns::AbstractVector{Symbol},
        values::AbstractMatrix{<:Real},
        line_numbers::AbstractVector{<:Integer},
        header_source::Symbol,
    )
        normalized_source = abspath(String(source))
        normalized_columns = Symbol.(columns)
        length(unique(normalized_columns)) == length(normalized_columns) ||
            throw(ArgumentError("CellParameterSeries columns must be unique"))

        matrix = Matrix{Float64}(values)
        size(matrix, 2) == length(normalized_columns) ||
            throw(ArgumentError(
                "number of matrix columns does not match the column-name count",
            ))
        size(matrix, 1) == length(line_numbers) ||
            throw(ArgumentError(
                "number of matrix rows does not match the line-number count",
            ))
        header_source in (:explicit, :file, :inferred) ||
            throw(ArgumentError(
                "header_source must be :explicit, :file, or :inferred",
            ))

        return new(
            normalized_source,
            normalized_columns,
            matrix,
            Int.(line_numbers),
            header_source,
        )
    end
end

Base.length(series::CellParameterSeries) = size(series.values, 1)
Base.size(series::CellParameterSeries) = size(series.values)
Base.size(series::CellParameterSeries, dimension::Integer) =
    size(series.values, dimension)
Base.isempty(series::CellParameterSeries) = length(series) == 0

function Base.show(io::IO, series::CellParameterSeries)
    print(
        io,
        "CellParameterSeries(",
        length(series),
        " samples, ",
        length(series.columns),
        " columns, header_source=:",
        series.header_source,
        ")",
    )
end

"""Return the canonical observable names in file order."""
cell_parameter_names(series::CellParameterSeries) = copy(series.columns)

"""Return whether a canonical observable exists in the series."""
has_cell_parameter(series::CellParameterSeries, name::Union{Symbol,AbstractString}) =
    _canonical_cell_parameter_name(name) in series.columns

"""
    cell_parameter_column(series, name)

Return a view of one observable column. The returned object shares storage with
`series.values` and should be copied before mutation.
"""
function cell_parameter_column(
    series::CellParameterSeries,
    name::Union{Symbol,AbstractString},
)
    canonical = _canonical_cell_parameter_name(name)
    index = findfirst(==(canonical), series.columns)
    isnothing(index) && throw(KeyError(canonical))
    return @view series.values[:, index]
end

Base.getindex(
    series::CellParameterSeries,
    name::Union{Symbol,AbstractString},
) = cell_parameter_column(series, name)

function _normalize_header_token(token::AbstractString)::String
    normalized = lowercase(strip(String(token)))
    normalized = replace(normalized, r"^[\$\{\[]+" => "")
    normalized = replace(normalized, r"[\}\],;:]+$" => "")
    normalized = replace(normalized, r"^(v|c|f)_" => "")
    normalized = replace(normalized, r"[^a-z0-9]+" => "_")
    normalized = strip(normalized, ['_'])
    return normalized
end

function _canonical_cell_parameter_name(
    name::Union{Symbol,AbstractString},
)::Symbol
    normalized = _normalize_header_token(String(name))
    isempty(normalized) && throw(ArgumentError("column name cannot be empty"))
    return get(_CELL_PARAMETER_ALIASES, normalized, Symbol(normalized))
end

function _canonicalize_columns(
    names::AbstractVector{<:Union{Symbol,AbstractString}},
    path::AbstractString,
    line::Integer,
)::Vector{Symbol}
    canonical = Symbol[_canonical_cell_parameter_name(name) for name in names]

    duplicates = Symbol[]
    for name in unique(canonical)
        count(==(name), canonical) > 1 && push!(duplicates, name)
    end

    isempty(duplicates) || throw(CellParameterParseError(
        String(path),
        Int(line),
        "duplicate columns after normalization: $(join(string.(duplicates), ", "))",
    ))

    return canonical
end

function _try_parse_row(tokens::AbstractVector{<:AbstractString})
    values = Float64[]
    sizehint!(values, length(tokens))

    for token in tokens
        value = tryparse(Float64, token)
        isnothing(value) && return nothing
        push!(values, value)
    end

    return values
end

function _is_step_like(values::AbstractVector{<:Real})::Bool
    length(values) >= 2 || return false
    all(isfinite, values) || return false
    all(value -> value >= 0.0, values) || return false
    all(value -> isapprox(value, round(value); atol=1.0e-8, rtol=1.0e-10), values) ||
        return false

    differences = diff(values)
    all(difference -> difference >= 0.0, differences) || return false
    any(difference -> difference > 0.0, differences) || return false
    return true
end

function _infer_legacy_columns(
    values::Matrix{Float64},
    path::AbstractString,
)::Vector{Symbol}
    nrows, ncolumns = size(values)
    nrows >= 1 || throw(CellParameterParseError(
        String(path),
        0,
        "cannot infer columns from an empty file",
    ))

    if ncolumns == length(_CELL_PARAMETER_CORE_COLUMNS)
        return collect(_CELL_PARAMETER_CORE_COLUMNS)
    end

    has_step = _is_step_like(@view values[:, 1])
    maximum_optional = length(_CELL_PARAMETER_OPTIONAL_COLUMNS)

    if has_step
        optional_count = ncolumns - 1 - length(_CELL_PARAMETER_CORE_COLUMNS)
        if 0 <= optional_count <= maximum_optional
            return vcat(
                [:step],
                collect(_CELL_PARAMETER_CORE_COLUMNS),
                collect(_CELL_PARAMETER_OPTIONAL_COLUMNS[1:optional_count]),
            )
        end
    else
        optional_count = ncolumns - length(_CELL_PARAMETER_CORE_COLUMNS)
        if 1 <= optional_count <= maximum_optional
            return vcat(
                collect(_CELL_PARAMETER_CORE_COLUMNS),
                collect(_CELL_PARAMETER_OPTIONAL_COLUMNS[1:optional_count]),
            )
        end
    end

    throw(CellParameterParseError(
        String(path),
        0,
        "headerless layout with $(ncolumns) columns is ambiguous; " *
        "supply `columns=[...]` explicitly",
    ))
end

function _header_matches(
    tokens::AbstractVector{<:AbstractString},
    columns::AbstractVector{Symbol},
)::Bool
    length(tokens) == length(columns) || return false
    try
        return _canonicalize_columns(tokens, "<header>", 0) == columns
    catch
        return false
    end
end

function _validate_required_columns(
    columns::AbstractVector{Symbol},
    required_columns,
    path::AbstractString,
)
    required = Symbol[
        _canonical_cell_parameter_name(name) for name in required_columns
    ]
    missing = setdiff(required, columns)
    isempty(missing) || throw(CellParameterParseError(
        String(path),
        0,
        "missing required columns: $(join(string.(missing), ", "))",
    ))
    return nothing
end

function _validate_physical_cell_values(
    values::Matrix{Float64},
    columns::Vector{Symbol},
    line_numbers::Vector{Int},
    path::AbstractString,
)
    for name in (:a, :b, :c, :volume)
        index = findfirst(==(name), columns)
        isnothing(index) && continue
        for row in axes(values, 1)
            value = values[row, index]
            isfinite(value) || continue
            value > 0.0 || throw(CellParameterParseError(
                String(path),
                line_numbers[row],
                "$(name) must be positive, found $(value)",
            ))
        end
    end

    for name in (:alpha, :beta, :gamma)
        index = findfirst(==(name), columns)
        isnothing(index) && continue
        for row in axes(values, 1)
            value = values[row, index]
            isfinite(value) || continue
            0.0 < value < 180.0 || throw(CellParameterParseError(
                String(path),
                line_numbers[row],
                "$(name) must be between 0 and 180 degrees, found $(value)",
            ))
        end
    end

    step_index = findfirst(==(:step), columns)
    if !isnothing(step_index)
        steps = @view values[:, step_index]
        for row in eachindex(steps)
            value = steps[row]
            isfinite(value) || continue
            value >= 0.0 || throw(CellParameterParseError(
                String(path),
                line_numbers[row],
                "step must be non-negative, found $(value)",
            ))
        end
        for row in 2:length(steps)
            previous = steps[row - 1]
            current = steps[row]
            if isfinite(previous) && isfinite(current) && current < previous
                throw(CellParameterParseError(
                    String(path),
                    line_numbers[row],
                    "step values must be nondecreasing ($(current) < $(previous))",
                ))
            end
        end
    end

    return nothing
end

"""
    read_cell_parameter_series(path; kwargs...) -> CellParameterSeries

Parse a whitespace-delimited `cell_params.dat`-style file.

Supported input forms:

1. a commented header, for example `# Step Cella Cellb ...`;
2. a plain-text header, for example `Step Cella Cellb ...`;
3. a headerless file with explicit `columns=[...]`;
4. a supported legacy headerless layout inferred from seven core cell columns,
   optionally preceded by a monotonic integer step column and followed by
   temperature, pressure, and energy columns.

Canonical required cell columns default to `a`, `b`, `c`, `alpha`, `beta`,
`gamma`, and `volume`.

Keyword arguments:

- `columns=nothing`: explicit names for a headerless file;
- `required_columns=_CELL_PARAMETER_CORE_COLUMNS`: columns that must exist;
- `allow_nonfinite=false`: permit `NaN`/`Inf` values for downstream failure
  diagnostics;
- `validate_physical=true`: check positive cell lengths/volume, valid angles,
  and nondecreasing steps;
- `allow_repeated_headers=true`: skip repeated headers in appended files.
"""
function read_cell_parameter_series(
    path::AbstractString;
    columns::Union{Nothing,AbstractVector}=nothing,
    required_columns=_CELL_PARAMETER_CORE_COLUMNS,
    allow_nonfinite::Bool=false,
    validate_physical::Bool=true,
    allow_repeated_headers::Bool=true,
)::CellParameterSeries
    normalized_path = abspath(String(path))
    isfile(normalized_path) || throw(CellParameterParseError(
        normalized_path,
        0,
        "file does not exist or is not a regular file",
    ))

    explicit_columns = if isnothing(columns)
        nothing
    else
        isempty(columns) && throw(ArgumentError("columns cannot be empty"))
        _canonicalize_columns(columns, normalized_path, 0)
    end

    rows = Vector{Vector{Float64}}()
    row_lines = Int[]
    file_header = nothing
    file_header_line = 0
    comment_header_candidates = Vector{Tuple{Int,Vector{String}}}()
    expected_columns = isnothing(explicit_columns) ? nothing : length(explicit_columns)

    open(normalized_path, "r") do io
        for (line_number, raw_line) in enumerate(eachline(io))
            stripped = strip(raw_line)
            isempty(stripped) && continue

            if startswith(stripped, "#") || startswith(stripped, ";")
                content = strip(chop(stripped; head=1, tail=0))
                isempty(content) && continue
                tokens = split(content)
                parsed = _try_parse_row(tokens)
                isnothing(parsed) && push!(
                    comment_header_candidates,
                    (line_number, String.(tokens)),
                )
                continue
            end

            tokens = split(stripped)
            parsed = _try_parse_row(tokens)

            if isnothing(parsed)
                if isempty(rows) && isnothing(explicit_columns) && isnothing(file_header)
                    file_header = String.(tokens)
                    file_header_line = line_number
                    expected_columns = length(tokens)
                    continue
                end

                active_columns = if !isnothing(explicit_columns)
                    explicit_columns
                elseif !isnothing(file_header)
                    _canonicalize_columns(
                        file_header,
                        normalized_path,
                        file_header_line,
                    )
                else
                    nothing
                end

                if allow_repeated_headers &&
                   !isnothing(active_columns) &&
                   _header_matches(tokens, active_columns)
                    continue
                end

                throw(CellParameterParseError(
                    normalized_path,
                    line_number,
                    "row contains non-numeric values and is not a recognized header",
                ))
            end

            if isnothing(expected_columns)
                expected_columns = length(parsed)
            elseif length(parsed) != expected_columns
                throw(CellParameterParseError(
                    normalized_path,
                    line_number,
                    "expected $(expected_columns) columns, found $(length(parsed))",
                ))
            end

            if !allow_nonfinite
                nonfinite_index = findfirst(value -> !isfinite(value), parsed)
                if !isnothing(nonfinite_index)
                    throw(CellParameterParseError(
                        normalized_path,
                        line_number,
                        "column $(nonfinite_index) contains a non-finite value",
                    ))
                end
            end

            push!(rows, parsed)
            push!(row_lines, line_number)
        end
    end

    isempty(rows) && throw(CellParameterParseError(
        normalized_path,
        0,
        "no numeric data rows were found",
    ))

    nrows = length(rows)
    ncolumns = length(rows[1])
    matrix = Matrix{Float64}(undef, nrows, ncolumns)
    for row in 1:nrows
        matrix[row, :] = rows[row]
    end

    local resolved_columns::Vector{Symbol}
    local header_source::Symbol

    if !isnothing(explicit_columns)
        length(explicit_columns) == ncolumns || throw(CellParameterParseError(
            normalized_path,
            0,
            "explicit column count $(length(explicit_columns)) does not match " *
            "the data width $(ncolumns)",
        ))
        resolved_columns = explicit_columns
        header_source = :explicit
    elseif !isnothing(file_header)
        length(file_header) == ncolumns || throw(CellParameterParseError(
            normalized_path,
            file_header_line,
            "header has $(length(file_header)) columns but data have $(ncolumns)",
        ))
        resolved_columns = _canonicalize_columns(
            file_header,
            normalized_path,
            file_header_line,
        )
        header_source = :file
    else
        matching_candidate = nothing
        for candidate in Iterators.reverse(comment_header_candidates)
            if length(candidate[2]) == ncolumns
                matching_candidate = candidate
                break
            end
        end

        if !isnothing(matching_candidate)
            resolved_columns = _canonicalize_columns(
                matching_candidate[2],
                normalized_path,
                matching_candidate[1],
            )
            header_source = :file
        else
            resolved_columns = _infer_legacy_columns(matrix, normalized_path)
            header_source = :inferred
        end
    end

    _validate_required_columns(
        resolved_columns,
        required_columns,
        normalized_path,
    )

    validate_physical && _validate_physical_cell_values(
        matrix,
        resolved_columns,
        row_lines,
        normalized_path,
    )

    return CellParameterSeries(
        normalized_path,
        resolved_columns,
        matrix,
        row_lines,
        header_source,
    )
end
