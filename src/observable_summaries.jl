# Simple statistical summaries for production time series.
#
# This layer computes means, raw standard deviations, and naive standard
# errors only. Block averaging, autocorrelation corrections, drift analysis,
# and between-cycle convergence remain separate later steps.

const _DEFAULT_OBSERVABLE_UNITS = Dict{Symbol,String}(
    :step => "step",
    :a => "Å",
    :b => "Å",
    :c => "Å",
    :alpha => "deg",
    :beta => "deg",
    :gamma => "deg",
    :volume => "Å^3",
    :temperature => "K",
    :density => "g/cm^3",
)

"""
    default_observable_unit(name)

Return the package default display unit for a canonical observable name.
Unknown observables, pressure, energies, and time return an empty string
because their units depend on the external simulation configuration.
"""
function default_observable_unit(
    name::Union{Symbol,AbstractString},
)::String
    canonical = _canonical_cell_parameter_name(name)
    return get(_DEFAULT_OBSERVABLE_UNITS, canonical, "")
end

function _summary_unit(
    name::Union{Symbol,AbstractString},
    unit::Union{Nothing,AbstractString},
)::String
    return isnothing(unit) ? default_observable_unit(name) : String(unit)
end

"""
    summarize_observable(values; name, stage=:md, unit=nothing,
                         minimum_samples=2, corrected=true)

Calculate the arithmetic mean, raw standard deviation, and naive standard
error of one production time series.

The result status is:

- `:summary_only` when at least `minimum_samples` finite values are available;
- `:insufficient_samples` when the series is finite but too short;
- `:invalid` when one or more values are non-finite.

This function deliberately does not mark an observable as `:valid`, because
block-average and drift diagnostics have not yet been evaluated.
"""
function summarize_observable(
    values::AbstractVector{<:Real};
    name::Union{Symbol,AbstractString},
    stage::Union{Symbol,AbstractString}=:md,
    unit::Union{Nothing,AbstractString}=nothing,
    minimum_samples::Integer=2,
    corrected::Bool=true,
)::ObservableStatistics
    minimum_samples >= 2 ||
        throw(ArgumentError("minimum_samples must be at least 2"))

    data = Float64.(collect(values))
    n_samples = length(data)
    resolved_unit = _summary_unit(name, unit)

    if n_samples == 0
        return ObservableStatistics(
            name=name,
            stage=stage,
            unit=resolved_unit,
            n_samples=0,
            status=:insufficient_samples,
            warnings=["time series contains no samples"],
        )
    end

    nonfinite_count = count(value -> !isfinite(value), data)
    if nonfinite_count > 0
        noun = nonfinite_count == 1 ? "value" : "values"
        return ObservableStatistics(
            name=name,
            stage=stage,
            unit=resolved_unit,
            n_samples=n_samples,
            status=:invalid,
            warnings=[
                "time series contains $(nonfinite_count) non-finite $(noun)",
            ],
        )
    end

    sample_mean = mean(data)
    sample_std = n_samples >= 2 ?
        std(data; corrected=corrected) :
        NaN
    naive_sem = n_samples >= 2 ?
        sample_std / sqrt(n_samples) :
        NaN

    if n_samples < minimum_samples
        return ObservableStatistics(
            name=name,
            stage=stage,
            unit=resolved_unit,
            n_samples=n_samples,
            mean=sample_mean,
            standard_deviation=sample_std,
            naive_sem=naive_sem,
            status=:insufficient_samples,
            warnings=[
                "time series contains $(n_samples) samples; " *
                "at least $(minimum_samples) are required",
            ],
        )
    end

    return ObservableStatistics(
        name=name,
        stage=stage,
        unit=resolved_unit,
        n_samples=n_samples,
        mean=sample_mean,
        standard_deviation=sample_std,
        naive_sem=naive_sem,
        status=:summary_only,
    )
end

function _normalize_summary_units(
    units::AbstractDict,
)::Dict{Symbol,String}
    normalized = Dict{Symbol,String}()
    for (raw_name, raw_unit) in units
        raw_name isa Union{Symbol,AbstractString} || throw(ArgumentError(
            "unit-map keys must be symbols or strings",
        ))
        raw_unit isa AbstractString || throw(ArgumentError(
            "unit-map values must be strings",
        ))
        canonical = _canonical_cell_parameter_name(raw_name)
        normalized[canonical] = String(raw_unit)
    end
    return normalized
end

function _selected_cell_parameter_columns(
    series::CellParameterSeries,
    columns,
    include_independent::Bool,
)::Vector{Symbol}
    if isnothing(columns)
        selected = copy(series.columns)
        if !include_independent
            filter!(name -> !(name in (:step, :time)), selected)
        end
        return selected
    end

    columns isa AbstractVector || throw(ArgumentError(
        "columns must be `nothing` or a vector of names",
    ))

    selected = Symbol[]
    for name in columns
        name isa Union{Symbol,AbstractString} || throw(ArgumentError(
            "column names must be symbols or strings",
        ))
        push!(selected, _canonical_cell_parameter_name(name))
    end

    length(unique(selected)) == length(selected) || throw(ArgumentError(
        "selected columns must be unique after canonicalization",
    ))

    return selected
end

"""
    summarize_cell_parameter_series(series; columns=nothing, units=Dict(),
                                    minimum_samples=2, corrected=true,
                                    include_independent=false)

Convert selected columns from a `CellParameterSeries` into
`ObservableStatistics` objects.

By default, all parsed observables except `:step` and `:time` are summarized.
Supplying `columns` selects an explicit subset and permits independent columns
such as `:step`. The returned dictionary is keyed by canonical observable name.

`units` can override display units, for example
`Dict(:pressure => "atm", :potential_energy => "kcal/mol")`.
"""
function summarize_cell_parameter_series(
    series::CellParameterSeries;
    columns=nothing,
    units::AbstractDict=Dict{Symbol,String}(),
    minimum_samples::Integer=2,
    corrected::Bool=true,
    include_independent::Bool=false,
)::Dict{Symbol,ObservableStatistics}
    selected = _selected_cell_parameter_columns(
        series,
        columns,
        include_independent,
    )

    isempty(selected) && throw(ArgumentError(
        "no cell-parameter columns were selected for summarization",
    ))

    for name in selected
        has_cell_parameter(series, name) || throw(KeyError(name))
    end

    unit_map = _normalize_summary_units(units)
    summaries = Dict{Symbol,ObservableStatistics}()

    for name in selected
        resolved_unit = get(unit_map, name, default_observable_unit(name))
        summaries[name] = summarize_observable(
            cell_parameter_column(series, name);
            name=name,
            stage=:md,
            unit=resolved_unit,
            minimum_samples=minimum_samples,
            corrected=corrected,
        )
    end

    return summaries
end

"""
    summarize_md_cycle(cycle, series; kwargs...)

Create a `CycleStatistics` object containing summary-only MD observables from a
parsed `cell_params.dat` series.

The cycle is intentionally returned with `valid == false`: mean, standard
deviation, and naive SEM alone are not sufficient to establish sampling
validity. Later steps will add block-average and drift diagnostics.
"""
function summarize_md_cycle(
    cycle::Integer,
    series::CellParameterSeries;
    kwargs...,
)::CycleStatistics
    summaries = summarize_cell_parameter_series(series; kwargs...)
    warnings = String[
        "MD observables contain summary-only statistics; " *
        "block averaging and drift analysis have not been evaluated",
    ]

    for name in sort!(collect(keys(summaries)); by=String)
        statistics = summaries[name]
        if statistics.status in (:invalid, :insufficient_samples)
            for warning in statistics.warnings
                push!(warnings, "$(name): $(warning)")
            end
        end
    end

    return CycleStatistics(
        cycle=cycle,
        md=summaries,
        valid=false,
        warnings=warnings,
    )
end
