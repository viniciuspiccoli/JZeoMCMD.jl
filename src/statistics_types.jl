# Typed result containers for the statistical and convergence layers.
#
# This file intentionally contains no trajectory parsing, block averaging, or
# convergence calculations. It defines the stable data model that those later
# steps will populate.

const _STATISTICAL_STAGES = (:gcmc, :md, :derived)
const _OBSERVABLE_STATUSES = (
    :not_evaluated,
    :summary_only,
    :valid,
    :insufficient_samples,
    :nonstationary,
    :invalid,
)

_nonempty_symbol(value::Symbol, field::AbstractString) = begin
    isempty(String(value)) && throw(ArgumentError("$(field) cannot be empty"))
    value
end

_nonempty_symbol(value::AbstractString, field::AbstractString) = begin
    text = strip(String(value))
    isempty(text) && throw(ArgumentError("$(field) cannot be empty"))
    Symbol(text)
end

function _finite_float(value::Real, field::AbstractString)::Float64
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("$(field) must be finite"))
    return result
end

function _nonnegative_float(value::Real, field::AbstractString)::Float64
    result = _finite_float(value, field)
    result >= 0.0 || throw(ArgumentError("$(field) must be non-negative"))
    return result
end

function _finite_or_nan(value::Real, field::AbstractString)::Float64
    result = Float64(value)
    (isfinite(result) || isnan(result)) ||
        throw(ArgumentError("$(field) must be finite or NaN"))
    return result
end

function _nonnegative_or_nan(value::Real, field::AbstractString)::Float64
    result = _finite_or_nan(value, field)
    (isnan(result) || result >= 0.0) ||
        throw(ArgumentError("$(field) must be non-negative or NaN"))
    return result
end

"""
    BlockAveragingStatistics

Diagnostics obtained from a block-average or autocorrelation analysis of one
production time series.

This type is only constructed when block analysis succeeds. If block analysis
has not yet been performed or failed, `ObservableStatistics.block` should be
`nothing` instead.
"""
struct BlockAveragingStatistics
    corrected_sem::Float64
    autocorrelation_time::Float64
    effective_sample_size::Float64
    block_size::Int
    number_of_blocks::Int
    plateau_detected::Bool
end

function BlockAveragingStatistics(;
    corrected_sem::Real,
    autocorrelation_time::Real,
    effective_sample_size::Real,
    block_size::Integer,
    number_of_blocks::Integer,
    plateau_detected::Bool,
)
    block_size >= 1 || throw(ArgumentError("block_size must be at least 1"))
    number_of_blocks >= 1 ||
        throw(ArgumentError("number_of_blocks must be at least 1"))

    return BlockAveragingStatistics(
        _nonnegative_float(corrected_sem, "corrected_sem"),
        _nonnegative_float(autocorrelation_time, "autocorrelation_time"),
        _nonnegative_float(effective_sample_size, "effective_sample_size"),
        Int(block_size),
        Int(number_of_blocks),
        plateau_detected,
    )
end

"""
    DriftStatistics

Stationarity diagnostics for one production time series.

`slope` and `slope_standard_error` describe a linear drift estimate.
`drift_score` is a non-negative normalized measure whose interpretation and
threshold are defined by the later analysis layer.
"""
struct DriftStatistics
    slope::Float64
    slope_standard_error::Float64
    early_mean::Float64
    late_mean::Float64
    drift_score::Float64
    stationary::Bool
end

function DriftStatistics(;
    slope::Real,
    slope_standard_error::Real,
    early_mean::Real,
    late_mean::Real,
    drift_score::Real,
    stationary::Bool,
)
    return DriftStatistics(
        _finite_float(slope, "slope"),
        _nonnegative_float(slope_standard_error, "slope_standard_error"),
        _finite_float(early_mean, "early_mean"),
        _finite_float(late_mean, "late_mean"),
        _nonnegative_float(drift_score, "drift_score"),
        stationary,
    )
end

"""
    ObservableStatistics

Statistical summary of one observable from one simulation stage.

Allowed stages are `:gcmc`, `:md`, and `:derived`. Allowed statuses are:

- `:not_evaluated`: no time-series analysis has been attempted;
- `:summary_only`: mean, standard deviation, and naive SEM are available;
- `:valid`: block and drift diagnostics are available and passed;
- `:insufficient_samples`: the series is too short or has too few effective
  samples;
- `:nonstationary`: significant within-cycle drift remains;
- `:invalid`: parsing or statistical evaluation failed.

The nested block and drift results are optional so that simple summaries can be
introduced before the full statistical analysis is implemented.
"""
struct ObservableStatistics
    name::Symbol
    stage::Symbol
    unit::String
    n_samples::Int
    mean::Float64
    standard_deviation::Float64
    naive_sem::Float64
    block::Union{Nothing,BlockAveragingStatistics}
    drift::Union{Nothing,DriftStatistics}
    status::Symbol
    warnings::Vector{String}
end

function ObservableStatistics(;
    name::Union{Symbol,AbstractString},
    stage::Union{Symbol,AbstractString},
    unit::AbstractString="",
    n_samples::Integer=0,
    mean::Real=NaN,
    standard_deviation::Real=NaN,
    naive_sem::Real=NaN,
    block::Union{Nothing,BlockAveragingStatistics}=nothing,
    drift::Union{Nothing,DriftStatistics}=nothing,
    status::Union{Symbol,AbstractString}=:not_evaluated,
    warnings::AbstractVector{<:AbstractString}=String[],
)
    observable_name = _nonempty_symbol(name, "name")
    observable_stage = _nonempty_symbol(stage, "stage")
    observable_stage in _STATISTICAL_STAGES ||
        throw(ArgumentError(
            "stage must be one of $(join(_STATISTICAL_STAGES, ", "))",
        ))

    observable_status = _nonempty_symbol(status, "status")
    observable_status in _OBSERVABLE_STATUSES ||
        throw(ArgumentError(
            "status must be one of $(join(_OBSERVABLE_STATUSES, ", "))",
        ))

    n_samples >= 0 || throw(ArgumentError("n_samples must be non-negative"))
    n = Int(n_samples)

    value_mean = _finite_or_nan(mean, "mean")
    value_std = _nonnegative_or_nan(standard_deviation, "standard_deviation")
    value_sem = _nonnegative_or_nan(naive_sem, "naive_sem")

    if observable_status == :not_evaluated
        n == 0 || throw(ArgumentError(
            "a :not_evaluated observable must have n_samples == 0",
        ))
        isnothing(block) || throw(ArgumentError(
            "a :not_evaluated observable cannot contain block diagnostics",
        ))
        isnothing(drift) || throw(ArgumentError(
            "a :not_evaluated observable cannot contain drift diagnostics",
        ))
    elseif observable_status in (:summary_only, :valid, :nonstationary)
        n >= 2 || throw(ArgumentError(
            "status $(observable_status) requires at least two samples",
        ))
        all(isfinite, (value_mean, value_std, value_sem)) ||
            throw(ArgumentError(
                "status $(observable_status) requires finite summary values",
            ))
    end

    if !isnothing(block) && n > 0
        block.block_size <= n || throw(ArgumentError(
            "block_size cannot exceed n_samples",
        ))
        block.block_size * block.number_of_blocks <= n ||
            throw(ArgumentError(
                "block_size * number_of_blocks cannot exceed n_samples",
            ))
    end

    if observable_status == :valid
        isnothing(block) && throw(ArgumentError(
            "status :valid requires block diagnostics",
        ))
        isnothing(drift) && throw(ArgumentError(
            "status :valid requires drift diagnostics",
        ))
        block.plateau_detected || throw(ArgumentError(
            "status :valid requires a detected block-error plateau",
        ))
        drift.stationary || throw(ArgumentError(
            "status :valid requires a stationary time series",
        ))
    elseif observable_status == :nonstationary
        isnothing(drift) && throw(ArgumentError(
            "status :nonstationary requires drift diagnostics",
        ))
        !drift.stationary || throw(ArgumentError(
            "status :nonstationary requires stationary == false",
        ))
    end

    return ObservableStatistics(
        observable_name,
        observable_stage,
        String(unit),
        n,
        value_mean,
        value_std,
        value_sem,
        block,
        drift,
        observable_status,
        String.(warnings),
    )
end

"""Return the best currently available standard error for an observable."""
statistical_uncertainty(statistics::ObservableStatistics) =
    isnothing(statistics.block) ? statistics.naive_sem : statistics.block.corrected_sem

"""Return `true` only when all required sampling diagnostics passed."""
sampling_valid(statistics::ObservableStatistics) = statistics.status == :valid

"""Return whether block-average diagnostics are available."""
has_block_analysis(statistics::ObservableStatistics) = !isnothing(statistics.block)

"""Return whether drift diagnostics are available."""
has_drift_analysis(statistics::ObservableStatistics) = !isnothing(statistics.drift)

function Base.show(io::IO, statistics::ObservableStatistics)
    print(io, "ObservableStatistics(", statistics.stage, ":", statistics.name,
          ", n=", statistics.n_samples, ", mean=", statistics.mean,
          ", uncertainty=", statistical_uncertainty(statistics),
          ", status=", statistics.status, ")")
end

function _normalize_observable_dictionary(
    observations::AbstractDict,
    expected_stage::Symbol,
)::Dict{Symbol,ObservableStatistics}
    normalized = Dict{Symbol,ObservableStatistics}()
    for (raw_name, statistics) in observations
        statistics isa ObservableStatistics || throw(ArgumentError(
            "all cycle observables must be ObservableStatistics objects",
        ))
        name = raw_name isa Symbol ? raw_name : Symbol(String(raw_name))
        name == statistics.name || throw(ArgumentError(
            "dictionary key $(name) does not match observable name " *
            "$(statistics.name)",
        ))
        statistics.stage == expected_stage || throw(ArgumentError(
            "observable $(name) has stage $(statistics.stage), expected " *
            "$(expected_stage)",
        ))
        normalized[name] = statistics
    end
    return normalized
end

"""
    CycleStatistics

All statistical results produced for one coupled GCMC/NPT-MD cycle.

The stage-specific dictionaries are keyed by observable name. `valid` is a
cycle-level sampling-quality flag; it is deliberately separate from iterative
between-cycle convergence.
"""
struct CycleStatistics
    cycle::Int
    gcmc::Dict{Symbol,ObservableStatistics}
    md::Dict{Symbol,ObservableStatistics}
    derived::Dict{Symbol,ObservableStatistics}
    valid::Bool
    warnings::Vector{String}
end

function CycleStatistics(;
    cycle::Integer,
    gcmc::AbstractDict=Dict{Symbol,ObservableStatistics}(),
    md::AbstractDict=Dict{Symbol,ObservableStatistics}(),
    derived::AbstractDict=Dict{Symbol,ObservableStatistics}(),
    valid::Bool=false,
    warnings::AbstractVector{<:AbstractString}=String[],
)
    cycle >= 1 || throw(ArgumentError("cycle must be at least 1"))

    gcmc_statistics = _normalize_observable_dictionary(gcmc, :gcmc)
    md_statistics = _normalize_observable_dictionary(md, :md)
    derived_statistics = _normalize_observable_dictionary(derived, :derived)

    if valid
        observations = vcat(
            collect(values(gcmc_statistics)),
            collect(values(md_statistics)),
            collect(values(derived_statistics)),
        )
        isempty(observations) && throw(ArgumentError(
            "a valid cycle must contain at least one observable",
        ))
        all(sampling_valid, observations) || throw(ArgumentError(
            "a valid cycle cannot contain observables that failed sampling checks",
        ))
    end

    return CycleStatistics(
        Int(cycle),
        gcmc_statistics,
        md_statistics,
        derived_statistics,
        valid,
        String.(warnings),
    )
end

"""Return the requested observable from a cycle, or `nothing` if absent."""
function get_observable(statistics::CycleStatistics,
                        stage::Union{Symbol,AbstractString},
                        name::Union{Symbol,AbstractString})
    stage_symbol = stage isa Symbol ? stage : Symbol(String(stage))
    name_symbol = name isa Symbol ? name : Symbol(String(name))

    stage_symbol == :gcmc && return get(statistics.gcmc, name_symbol, nothing)
    stage_symbol == :md && return get(statistics.md, name_symbol, nothing)
    stage_symbol == :derived && return get(statistics.derived, name_symbol, nothing)
    throw(ArgumentError("stage must be one of $(join(_STATISTICAL_STAGES, ", "))"))
end

"""Return all observables contained in a cycle."""
function all_observables(statistics::CycleStatistics)
    return vcat(
        collect(values(statistics.gcmc)),
        collect(values(statistics.md)),
        collect(values(statistics.derived)),
    )
end

"""Return the cycle-level sampling-quality flag."""
cycle_statistics_valid(statistics::CycleStatistics) = statistics.valid

function Base.show(io::IO, statistics::CycleStatistics)
    print(io, "CycleStatistics(cycle=", statistics.cycle,
          ", observables=", length(all_observables(statistics)),
          ", valid=", statistics.valid, ")")
end

"""
    ConvergenceDecision

Result of a future between-cycle convergence evaluation.

This type does not implement a convergence criterion. It records the overall
result, the number of consecutive passing cycles, per-observable pass/fail
flags, and explanatory reasons.
"""
struct ConvergenceDecision
    cycle::Int
    converged::Bool
    consecutive_count::Int
    required_consecutive::Int
    observable_results::Dict{Symbol,Bool}
    reasons::Vector{String}
end

function ConvergenceDecision(;
    cycle::Integer,
    converged::Bool,
    consecutive_count::Integer,
    required_consecutive::Integer,
    observable_results::AbstractDict=Dict{Symbol,Bool}(),
    reasons::AbstractVector{<:AbstractString}=String[],
)
    cycle >= 1 || throw(ArgumentError("cycle must be at least 1"))
    consecutive_count >= 0 ||
        throw(ArgumentError("consecutive_count must be non-negative"))
    required_consecutive >= 1 ||
        throw(ArgumentError("required_consecutive must be at least 1"))

    normalized_results = Dict{Symbol,Bool}()
    for (raw_name, result) in observable_results
        name = raw_name isa Symbol ? raw_name : Symbol(String(raw_name))
        result isa Bool || throw(ArgumentError(
            "observable convergence results must be Bool values",
        ))
        normalized_results[name] = result
    end

    if converged
        consecutive_count >= required_consecutive || throw(ArgumentError(
            "a converged decision requires enough consecutive passing cycles",
        ))
        isempty(normalized_results) && throw(ArgumentError(
            "a converged decision requires at least one observable result",
        ))
        all(values(normalized_results)) || throw(ArgumentError(
            "a converged decision cannot contain failed observables",
        ))
    end

    return ConvergenceDecision(
        Int(cycle),
        converged,
        Int(consecutive_count),
        Int(required_consecutive),
        normalized_results,
        String.(reasons),
    )
end

"""Return the sorted names of observables that did not pass convergence."""
function failed_observables(decision::ConvergenceDecision)
    failed = Symbol[
        name for (name, passed) in decision.observable_results if !passed
    ]
    sort!(failed; by=String)
    return failed
end

function Base.show(io::IO, decision::ConvergenceDecision)
    print(io, "ConvergenceDecision(cycle=", decision.cycle,
          ", converged=", decision.converged,
          ", consecutive=", decision.consecutive_count, "/",
          decision.required_consecutive, ")")
end
