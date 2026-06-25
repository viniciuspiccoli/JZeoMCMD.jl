function valid_block_statistics()
    return BlockAveragingStatistics(
        corrected_sem=0.12,
        autocorrelation_time=3.5,
        effective_sample_size=40.0,
        block_size=10,
        number_of_blocks=10,
        plateau_detected=true,
    )
end

function stationary_drift_statistics()
    return DriftStatistics(
        slope=1.0e-6,
        slope_standard_error=2.0e-6,
        early_mean=100.0,
        late_mean=100.01,
        drift_score=0.25,
        stationary=true,
    )
end

function valid_observable(name::Symbol, stage::Symbol; unit::String="")
    return ObservableStatistics(
        name=name,
        stage=stage,
        unit=unit,
        n_samples=100,
        mean=100.0,
        standard_deviation=1.2,
        naive_sem=0.12,
        block=valid_block_statistics(),
        drift=stationary_drift_statistics(),
        status=:valid,
    )
end

@testset "Statistical diagnostic result types" begin
    block = valid_block_statistics()
    @test block.corrected_sem == 0.12
    @test block.autocorrelation_time == 3.5
    @test block.effective_sample_size == 40.0
    @test block.block_size == 10
    @test block.number_of_blocks == 10
    @test block.plateau_detected

    drift = stationary_drift_statistics()
    @test drift.stationary
    @test drift.drift_score == 0.25

    @test_throws ArgumentError BlockAveragingStatistics(
        corrected_sem=-1.0,
        autocorrelation_time=1.0,
        effective_sample_size=10.0,
        block_size=2,
        number_of_blocks=5,
        plateau_detected=true,
    )
    @test_throws ArgumentError BlockAveragingStatistics(
        corrected_sem=0.1,
        autocorrelation_time=1.0,
        effective_sample_size=10.0,
        block_size=0,
        number_of_blocks=5,
        plateau_detected=true,
    )
    @test_throws ArgumentError DriftStatistics(
        slope=0.0,
        slope_standard_error=-0.1,
        early_mean=1.0,
        late_mean=1.0,
        drift_score=0.0,
        stationary=true,
    )
end

@testset "Observable statistics" begin
    summary = ObservableStatistics(
        name=:volume,
        stage=:md,
        unit="A^3",
        n_samples=100,
        mean=45_000.0,
        standard_deviation=50.0,
        naive_sem=5.0,
        status=:summary_only,
    )

    @test summary.name == :volume
    @test summary.stage == :md
    @test summary.unit == "A^3"
    @test summary.status == :summary_only
    @test !sampling_valid(summary)
    @test !has_block_analysis(summary)
    @test !has_drift_analysis(summary)
    @test statistical_uncertainty(summary) == 5.0

    complete = valid_observable(:volume, :md; unit="A^3")
    @test sampling_valid(complete)
    @test has_block_analysis(complete)
    @test has_drift_analysis(complete)
    @test statistical_uncertainty(complete) == 0.12
    @test occursin("md:volume", sprint(show, complete))

    pending = ObservableStatistics(name="loading", stage="gcmc")
    @test pending.name == :loading
    @test pending.stage == :gcmc
    @test pending.status == :not_evaluated
    @test pending.n_samples == 0
    @test isnan(pending.mean)

    original_warnings = ["low effective sample size"]
    warned = ObservableStatistics(
        name=:loading,
        stage=:gcmc,
        n_samples=10,
        mean=1.0,
        standard_deviation=0.2,
        naive_sem=0.05,
        status=:summary_only,
        warnings=original_warnings,
    )
    push!(original_warnings, "later mutation")
    @test warned.warnings == ["low effective sample size"]

    @test_throws ArgumentError ObservableStatistics(
        name=:volume,
        stage=:unknown,
    )
    @test_throws ArgumentError ObservableStatistics(
        name=:volume,
        stage=:md,
        status=:unknown,
    )
    @test_throws ArgumentError ObservableStatistics(
        name=:volume,
        stage=:md,
        n_samples=1,
        mean=1.0,
        standard_deviation=0.0,
        naive_sem=0.0,
        status=:summary_only,
    )
    @test_throws ArgumentError ObservableStatistics(
        name=:volume,
        stage=:md,
        n_samples=100,
        mean=1.0,
        standard_deviation=0.1,
        naive_sem=0.01,
        status=:valid,
    )

    nonstationary = DriftStatistics(
        slope=0.1,
        slope_standard_error=0.01,
        early_mean=1.0,
        late_mean=2.0,
        drift_score=8.0,
        stationary=false,
    )
    drifting = ObservableStatistics(
        name=:volume,
        stage=:md,
        n_samples=100,
        mean=1.5,
        standard_deviation=0.3,
        naive_sem=0.03,
        drift=nonstationary,
        status=:nonstationary,
    )
    @test drifting.status == :nonstationary
    @test !drifting.drift.stationary
end

@testset "Cycle statistics" begin
    loading = valid_observable(:loading, :gcmc; unit="molecules/uc")
    volume = valid_observable(:volume, :md; unit="A^3")
    strain = valid_observable(:volumetric_strain, :derived)

    cycle = CycleStatistics(
        cycle=4,
        gcmc=Dict(:loading => loading),
        md=Dict("volume" => volume),
        derived=Dict(:volumetric_strain => strain),
        valid=true,
        warnings=["example cycle warning"],
    )

    @test cycle.cycle == 4
    @test cycle_statistics_valid(cycle)
    @test get_observable(cycle, :gcmc, :loading) === loading
    @test get_observable(cycle, "md", "volume") === volume
    @test get_observable(cycle, :md, :missing) === nothing
    @test length(all_observables(cycle)) == 3
    @test occursin("cycle=4", sprint(show, cycle))

    @test_throws ArgumentError CycleStatistics(cycle=0)
    @test_throws ArgumentError CycleStatistics(
        cycle=1,
        md=Dict(:wrong_key => volume),
    )
    @test_throws ArgumentError CycleStatistics(
        cycle=1,
        gcmc=Dict(:volume => volume),
    )
    @test_throws ArgumentError CycleStatistics(cycle=1, valid=true)
    @test_throws ArgumentError CycleStatistics(
        cycle=1,
        md=Dict(:volume => ObservableStatistics(
            name=:volume,
            stage=:md,
            n_samples=10,
            mean=1.0,
            standard_deviation=0.1,
            naive_sem=0.03,
            status=:summary_only,
        )),
        valid=true,
    )
    @test_throws ArgumentError get_observable(cycle, :unknown, :volume)
end

@testset "Convergence decision data model" begin
    decision = ConvergenceDecision(
        cycle=8,
        converged=false,
        consecutive_count=2,
        required_consecutive=3,
        observable_results=Dict(
            :loading => true,
            :volume => true,
            :cell_b => false,
        ),
        reasons=["cell_b remains outside tolerance"],
    )

    @test !decision.converged
    @test failed_observables(decision) == [:cell_b]
    @test occursin("consecutive=2/3", sprint(show, decision))

    converged = ConvergenceDecision(
        cycle=9,
        converged=true,
        consecutive_count=3,
        required_consecutive=3,
        observable_results=Dict(:loading => true, :volume => true),
    )
    @test converged.converged
    @test isempty(failed_observables(converged))

    @test_throws ArgumentError ConvergenceDecision(
        cycle=1,
        converged=true,
        consecutive_count=2,
        required_consecutive=3,
        observable_results=Dict(:volume => true),
    )
    @test_throws ArgumentError ConvergenceDecision(
        cycle=1,
        converged=true,
        consecutive_count=3,
        required_consecutive=3,
        observable_results=Dict(:volume => false),
    )
    @test_throws ArgumentError ConvergenceDecision(
        cycle=1,
        converged=true,
        consecutive_count=3,
        required_consecutive=3,
    )
end
