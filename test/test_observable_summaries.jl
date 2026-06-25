@testset "Observable summary calculations" begin
    statistics = summarize_observable(
        [1.0, 2.0, 3.0, 4.0];
        name=:volume,
        stage=:md,
        unit="A^3",
    )

    @test statistics.status == :summary_only
    @test statistics.n_samples == 4
    @test statistics.mean == 2.5
    @test isapprox(
        statistics.standard_deviation,
        sqrt(5.0 / 3.0);
        atol=1.0e-12,
    )
    @test isapprox(
        statistics.naive_sem,
        sqrt(5.0 / 3.0) / 2.0;
        atol=1.0e-12,
    )
    @test statistical_uncertainty(statistics) == statistics.naive_sem
    @test !sampling_valid(statistics)

    constant = summarize_observable(
        fill(7.5, 10);
        name=:temperature,
        stage=:md,
    )
    @test constant.mean == 7.5
    @test constant.standard_deviation == 0.0
    @test constant.naive_sem == 0.0
    @test constant.unit == "K"

    population = summarize_observable(
        [1.0, 2.0, 3.0, 4.0];
        name=:custom,
        stage=:derived,
        corrected=false,
    )
    @test isapprox(
        population.standard_deviation,
        sqrt(5.0 / 4.0);
        atol=1.0e-12,
    )

    short = summarize_observable(
        [42.0];
        name=:loading,
        stage=:gcmc,
    )
    @test short.status == :insufficient_samples
    @test short.n_samples == 1
    @test short.mean == 42.0
    @test isnan(short.standard_deviation)
    @test isnan(short.naive_sem)
    @test !isempty(short.warnings)

    empty_summary = summarize_observable(
        Float64[];
        name=:loading,
        stage=:gcmc,
    )
    @test empty_summary.status == :insufficient_samples
    @test empty_summary.n_samples == 0
    @test isnan(empty_summary.mean)

    invalid = summarize_observable(
        [1.0, Inf, 2.0, NaN];
        name=:pressure,
        stage=:md,
    )
    @test invalid.status == :invalid
    @test invalid.n_samples == 4
    @test any(warning -> occursin("2 non-finite", warning), invalid.warnings)

    @test default_observable_unit(:a) == "Å"
    @test default_observable_unit(:volume) == "Å^3"
    @test default_observable_unit(:pressure) == ""
    @test default_observable_unit("CellAlpha") == "deg"

    @test_throws ArgumentError summarize_observable(
        [1.0, 2.0];
        name=:volume,
        minimum_samples=1,
    )
end

@testset "Cell-parameter series summaries" begin
    mktempdir() do root
        path = joinpath(root, "cell_params.dat")
        write(
            path,
            """
            # Step Cella Cellb Cellc CellAlpha CellBeta CellGamma Volume Temp Press PotEng
            0   40.0 41.0 27.0 90.0 90.0 90.0 44280.0 372.8 1.2 -1000.0
            100 40.1 41.1 27.1 89.9 90.1 90.0 44655.4 373.2 0.8 -1001.0
            200 40.2 41.2 27.2 89.8 90.2 90.0 45034.0 373.0 1.1 -999.5
            """,
        )

        series = read_cell_parameter_series(path)
        summaries = summarize_cell_parameter_series(
            series;
            units=Dict(
                :pressure => "atm",
                "potential_energy" => "kcal/mol",
            ),
        )

        @test !haskey(summaries, :step)
        @test length(summaries) == 10
        @test summaries[:a].unit == "Å"
        @test summaries[:volume].unit == "Å^3"
        @test summaries[:pressure].unit == "atm"
        @test summaries[:potential_energy].unit == "kcal/mol"
        @test summaries[:volume].status == :summary_only
        @test summaries[:volume].n_samples == 3
        @test isapprox(
            summaries[:temperature].mean,
            373.0;
            atol=1.0e-12,
        )

        selected = summarize_cell_parameter_series(
            series;
            columns=["vol", :temp],
            units=Dict(:volume => "custom-volume"),
        )
        @test Set(keys(selected)) == Set([:volume, :temperature])
        @test selected[:volume].unit == "custom-volume"

        with_step = summarize_cell_parameter_series(
            series;
            include_independent=true,
        )
        @test haskey(with_step, :step)
        @test with_step[:step].unit == "step"

        explicit_step = summarize_cell_parameter_series(
            series;
            columns=[:step],
        )
        @test Set(keys(explicit_step)) == Set([:step])

        cycle = summarize_md_cycle(
            3,
            series;
            columns=[:volume, :temperature],
        )
        @test cycle.cycle == 3
        @test !cycle.valid
        @test Set(keys(cycle.md)) == Set([:volume, :temperature])
        @test occursin("summary-only", cycle.warnings[1])

        @test_throws KeyError summarize_cell_parameter_series(
            series;
            columns=[:missing],
        )
        @test_throws ArgumentError summarize_cell_parameter_series(
            series;
            columns=[:volume, "vol"],
        )
        @test_throws ArgumentError summarize_cell_parameter_series(
            series;
            units=Dict(:volume => 123),
        )
        @test_throws ArgumentError summarize_md_cycle(
            0,
            series;
            columns=[:volume],
        )
    end
end
