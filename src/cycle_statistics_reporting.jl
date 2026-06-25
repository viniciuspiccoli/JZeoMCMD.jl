# Per-cycle statistical report generation.
#
# This layer integrates the simple production summaries into the simulation
# workflow. It deliberately leaves the legacy between-cycle convergence
# decision unchanged. Block averaging, autocorrelation corrections, and drift
# diagnostics are added in later steps.

import Dates
import JSON

"""
    CycleStatisticsReport

Files and typed statistics generated for one simulation cycle.

`csv_path` contains one row per observable. `json_path` contains the complete
cycle-level record and is intended for programmatic post-processing.
"""
struct CycleStatisticsReport
    statistics::CycleStatistics
    directory::String
    csv_path::String
    json_path::String
    source_path::String
end

"""
    CycleStatisticsReportError

Raised when per-cycle statistics cannot be generated in strict mode.
"""
struct CycleStatisticsReportError <: Exception
    message::String
end

Base.showerror(io::IO, error::CycleStatisticsReportError) =
    print(io, error.message)

_json_statistic(value::Real) =
    isfinite(Float64(value)) ? Float64(value) : nothing

function _block_statistics_dictionary(
    block::Union{Nothing,BlockAveragingStatistics},
)
    isnothing(block) && return nothing
    return Dict{String,Any}(
        "corrected_sem" => _json_statistic(block.corrected_sem),
        "autocorrelation_time" =>
            _json_statistic(block.autocorrelation_time),
        "effective_sample_size" =>
            _json_statistic(block.effective_sample_size),
        "block_size" => block.block_size,
        "number_of_blocks" => block.number_of_blocks,
        "plateau_detected" => block.plateau_detected,
    )
end

function _drift_statistics_dictionary(
    drift::Union{Nothing,DriftStatistics},
)
    isnothing(drift) && return nothing
    return Dict{String,Any}(
        "slope" => _json_statistic(drift.slope),
        "slope_standard_error" =>
            _json_statistic(drift.slope_standard_error),
        "early_mean" => _json_statistic(drift.early_mean),
        "late_mean" => _json_statistic(drift.late_mean),
        "drift_score" => _json_statistic(drift.drift_score),
        "stationary" => drift.stationary,
    )
end

"""
    observable_statistics_dictionary(statistics)

Convert one `ObservableStatistics` object into a JSON-compatible dictionary.
Non-finite floating-point values are represented as `nothing` and therefore
written as JSON `null`.
"""
function observable_statistics_dictionary(
    statistics::ObservableStatistics,
)::Dict{String,Any}
    return Dict{String,Any}(
        "name" => String(statistics.name),
        "stage" => String(statistics.stage),
        "unit" => statistics.unit,
        "n_samples" => statistics.n_samples,
        "mean" => _json_statistic(statistics.mean),
        "standard_deviation" =>
            _json_statistic(statistics.standard_deviation),
        "naive_sem" => _json_statistic(statistics.naive_sem),
        "block" => _block_statistics_dictionary(statistics.block),
        "drift" => _drift_statistics_dictionary(statistics.drift),
        "status" => String(statistics.status),
        "warnings" => copy(statistics.warnings),
    )
end

function _stage_statistics_dictionary(
    observations::Dict{Symbol,ObservableStatistics},
)::Dict{String,Any}
    result = Dict{String,Any}()
    for name in sort!(collect(keys(observations)); by=String)
        result[String(name)] =
            observable_statistics_dictionary(observations[name])
    end
    return result
end

"""
    cycle_statistics_dictionary(statistics; source_path="")

Convert `CycleStatistics` into a JSON-compatible dictionary.
"""
function cycle_statistics_dictionary(
    statistics::CycleStatistics;
    source_path::AbstractString="",
)::Dict{String,Any}
    return Dict{String,Any}(
        "schema_version" => 1,
        "generated_at" => string(Dates.now()),
        "cycle" => statistics.cycle,
        "valid" => statistics.valid,
        "source_path" => String(source_path),
        "warnings" => copy(statistics.warnings),
        "gcmc" => _stage_statistics_dictionary(statistics.gcmc),
        "md" => _stage_statistics_dictionary(statistics.md),
        "derived" => _stage_statistics_dictionary(statistics.derived),
    )
end

function _csv_field(value)::String
    text = string(value)
    if occursin(',', text) || occursin('"', text) ||
       occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

_csv_optional(value::Real) =
    isfinite(Float64(value)) ? string(Float64(value)) : ""

function _statistics_csv_row(
    cycle::Integer,
    statistics::ObservableStatistics,
)::Vector{String}
    block = statistics.block
    drift = statistics.drift

    return String[
        string(cycle),
        String(statistics.stage),
        String(statistics.name),
        statistics.unit,
        string(statistics.n_samples),
        _csv_optional(statistics.mean),
        _csv_optional(statistics.standard_deviation),
        _csv_optional(statistics.naive_sem),
        isnothing(block) ? "" : _csv_optional(block.corrected_sem),
        isnothing(block) ? "" :
            _csv_optional(block.autocorrelation_time),
        isnothing(block) ? "" :
            _csv_optional(block.effective_sample_size),
        isnothing(block) ? "" : string(block.block_size),
        isnothing(block) ? "" : string(block.number_of_blocks),
        isnothing(block) ? "" : string(block.plateau_detected),
        isnothing(drift) ? "" : _csv_optional(drift.slope),
        isnothing(drift) ? "" :
            _csv_optional(drift.slope_standard_error),
        isnothing(drift) ? "" : _csv_optional(drift.early_mean),
        isnothing(drift) ? "" : _csv_optional(drift.late_mean),
        isnothing(drift) ? "" : _csv_optional(drift.drift_score),
        isnothing(drift) ? "" : string(drift.stationary),
        String(statistics.status),
        join(statistics.warnings, " | "),
    ]
end

const _CYCLE_STATISTICS_CSV_HEADER = String[
    "cycle",
    "stage",
    "observable",
    "unit",
    "n_samples",
    "mean",
    "standard_deviation",
    "naive_sem",
    "corrected_sem",
    "autocorrelation_time",
    "effective_sample_size",
    "block_size",
    "number_of_blocks",
    "plateau_detected",
    "slope",
    "slope_standard_error",
    "early_mean",
    "late_mean",
    "drift_score",
    "stationary",
    "status",
    "warnings",
]

function _ordered_cycle_observables(
    statistics::CycleStatistics,
)::Vector{ObservableStatistics}
    stage_order = Dict(:gcmc => 1, :md => 2, :derived => 3)
    observations = all_observables(statistics)
    sort!(
        observations;
        by=observation -> (
            get(stage_order, observation.stage, typemax(Int)),
            String(observation.name),
        ),
    )
    return observations
end

"""
    write_cycle_statistics_csv(path, statistics; overwrite=true)

Write a tidy CSV file with one row per observable. Fields reserved for block
and drift diagnostics are included now and remain empty until those analyses
are implemented.
"""
function write_cycle_statistics_csv(
    path::AbstractString,
    statistics::CycleStatistics;
    overwrite::Bool=true,
)::String
    destination = normpath(abspath(expanduser(String(path))))
    if ispath(destination) && !overwrite
        throw(CycleStatisticsReportError(
            "statistics CSV already exists: $(destination)",
        ))
    end

    mkpath(dirname(destination))
    open(destination, "w") do io
        println(io, join(_csv_field.(_CYCLE_STATISTICS_CSV_HEADER), ","))
        for observation in _ordered_cycle_observables(statistics)
            row = _statistics_csv_row(statistics.cycle, observation)
            println(io, join(_csv_field.(row), ","))
        end
    end
    return destination
end

"""
    write_cycle_statistics_json(path, statistics; source_path="", overwrite=true)

Write the complete typed cycle summary as standards-compliant JSON.
"""
function write_cycle_statistics_json(
    path::AbstractString,
    statistics::CycleStatistics;
    source_path::AbstractString="",
    overwrite::Bool=true,
)::String
    destination = normpath(abspath(expanduser(String(path))))
    if ispath(destination) && !overwrite
        throw(CycleStatisticsReportError(
            "statistics JSON already exists: $(destination)",
        ))
    end

    document = cycle_statistics_dictionary(
        statistics;
        source_path=source_path,
    )

    mkpath(dirname(destination))
    open(destination, "w") do io
        JSON.print(io, document, 2)
        println(io)
    end
    return destination
end

"""
    write_cycle_statistics_report(statistics, directory; kwargs...)

Write `md_statistics.csv` and `cycle_statistics.json` in `directory`.
"""
function write_cycle_statistics_report(
    statistics::CycleStatistics,
    directory::AbstractString;
    source_path::AbstractString="",
    csv_filename::AbstractString="md_statistics.csv",
    json_filename::AbstractString="cycle_statistics.json",
    overwrite::Bool=true,
)::CycleStatisticsReport
    report_directory =
        normpath(abspath(expanduser(String(directory))))
    mkpath(report_directory)

    csv_path = write_cycle_statistics_csv(
        joinpath(report_directory, String(csv_filename)),
        statistics;
        overwrite=overwrite,
    )
    json_path = write_cycle_statistics_json(
        joinpath(report_directory, String(json_filename)),
        statistics;
        source_path=source_path,
        overwrite=overwrite,
    )

    return CycleStatisticsReport(
        statistics,
        report_directory,
        csv_path,
        json_path,
        String(source_path),
    )
end

function _summary_display_name(name::Symbol)::String
    names = Dict(
        :a => "a",
        :b => "b",
        :c => "c",
        :alpha => "alpha",
        :beta => "beta",
        :gamma => "gamma",
        :volume => "volume",
        :temperature => "temperature",
        :pressure => "pressure",
        :density => "density",
        :potential_energy => "potential energy",
        :total_energy => "total energy",
        :enthalpy => "enthalpy",
    )
    return get(names, name, replace(String(name), "_" => " "))
end

function _print_observable_summary(
    io::IO,
    statistics::ObservableStatistics,
)
    label = _summary_display_name(statistics.name)
    if isfinite(statistics.mean) && isfinite(statistics.naive_sem)
        unit_suffix = isempty(statistics.unit) ? "" : " $(statistics.unit)"
        @printf(
            io,
            "    %-20s %14.6g ± %-12.4g%s  n=%d  [%s]\n",
            label,
            statistics.mean,
            statistics.naive_sem,
            unit_suffix,
            statistics.n_samples,
            String(statistics.status),
        )
    else
        @printf(
            io,
            "    %-20s %-30s n=%d  [%s]\n",
            label,
            "statistics unavailable",
            statistics.n_samples,
            String(statistics.status),
        )
    end
end

"""
    print_md_cycle_statistics([io], statistics)

Print a concise summary of the MD production statistics. This report explicitly
labels uncertainties as naive SEM values until block averaging is available.
"""
function print_md_cycle_statistics(
    io::IO,
    statistics::CycleStatistics,
)
    println(io, "  [STATS] MD production summary (naive SEM)")
    preferred = Symbol[
        :a, :b, :c, :alpha, :beta, :gamma,
        :volume, :temperature, :pressure,
    ]
    printed = Set{Symbol}()

    for name in preferred
        observation = get(statistics.md, name, nothing)
        isnothing(observation) && continue
        _print_observable_summary(io, observation)
        push!(printed, name)
    end

    for name in sort!(collect(keys(statistics.md)); by=String)
        name in printed && continue
        _print_observable_summary(io, statistics.md[name])
    end

    if isempty(statistics.md)
        println(io, "    no MD observables were available")
    end

    for warning in statistics.warnings
        println(io, "    WARNING: ", warning)
    end
    return nothing
end

print_md_cycle_statistics(
    statistics::CycleStatistics,
) = print_md_cycle_statistics(stdout, statistics)

function _failed_md_cycle_statistics(
    cycle::Integer,
    message::AbstractString,
)::CycleStatistics
    return CycleStatistics(
        cycle=cycle,
        valid=false,
        warnings=[String(message)],
    )
end

"""
    analyze_md_cycle_statistics(cycle, lammps_directory; kwargs...)
        -> CycleStatisticsReport

Parse `cell_params.dat`, calculate simple MD production summaries, write the
per-cycle CSV/JSON report, and print a concise console summary.

The default report directory is the sibling `statistics/` directory of the
cycle's `lammps/` directory.

When `strict=false` (default), a missing or malformed time-series file produces
a report containing warnings but does not abort the simulation. With
`strict=true`, the underlying failure is rethrown as
`CycleStatisticsReportError`.
"""
function analyze_md_cycle_statistics(
    cycle::Integer,
    lammps_directory::AbstractString;
    filename::AbstractString="cell_params.dat",
    report_directory::Union{Nothing,AbstractString}=nothing,
    columns=nothing,
    units::AbstractDict=Dict{Symbol,String}(),
    minimum_samples::Integer=2,
    corrected::Bool=true,
    strict::Bool=false,
    overwrite::Bool=true,
    print_summary::Bool=true,
)::CycleStatisticsReport
    cycle >= 1 || throw(ArgumentError("cycle must be at least 1"))

    lammps_dir =
        normpath(abspath(expanduser(String(lammps_directory))))
    source_path = joinpath(lammps_dir, String(filename))
    statistics_dir = isnothing(report_directory) ?
        joinpath(dirname(lammps_dir), "statistics") :
        normpath(abspath(expanduser(String(report_directory))))

    statistics = if !isfile(source_path)
        message = "MD time-series file not found: $(source_path)"
        strict && throw(CycleStatisticsReportError(message))
        _failed_md_cycle_statistics(cycle, message)
    else
        try
            series = read_cell_parameter_series(source_path)
            summarize_md_cycle(
                cycle,
                series;
                columns=columns,
                units=units,
                minimum_samples=minimum_samples,
                corrected=corrected,
            )
        catch error
            message = "failed to analyze $(source_path): " *
                      sprint(showerror, error)
            strict && throw(CycleStatisticsReportError(message))
            _failed_md_cycle_statistics(cycle, message)
        end
    end

    report = write_cycle_statistics_report(
        statistics,
        statistics_dir;
        source_path=source_path,
        overwrite=overwrite,
    )
    print_summary && print_md_cycle_statistics(statistics)
    return report
end
