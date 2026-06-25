@testset "Cell-parameter time-series parser" begin
    mktempdir() do root
        path = joinpath(root, "cell_params.dat")
        write(path, """
# Production cell parameters
# Step Cella Cellb Cellc CellAlpha CellBeta CellGamma Volume Temp Press PotEng TotEng Enthalpy
0 40.0 41.0 27.0 90.0 90.0 90.0 44280.0 372.8 1.2 -1000.0 -990.0 -989.0
100 40.1 41.1 27.1 89.9 90.1 90.0 44655.381 373.2 0.8 -1001.0 -991.0 -990.0
""")

        series = read_cell_parameter_series(path)

        @test length(series) == 2
        @test size(series) == (2, 13)
        @test series.header_source == :file
        @test series.source == abspath(path)
        @test series.line_numbers == [3, 4]
        @test cell_parameter_names(series) == [
            :step,
            :a,
            :b,
            :c,
            :alpha,
            :beta,
            :gamma,
            :volume,
            :temperature,
            :pressure,
            :potential_energy,
            :total_energy,
            :enthalpy,
        ]
        @test has_cell_parameter(series, :volume)
        @test has_cell_parameter(series, "CellAlpha")
        @test !has_cell_parameter(series, :density)
        @test collect(series[:a]) == [40.0, 40.1]
        @test collect(cell_parameter_column(series, "Temp")) == [372.8, 373.2]
        @test_throws KeyError series[:density]
        @test occursin("2 samples", sprint(show, series))
    end
end

@testset "Plain and repeated cell-parameter headers" begin
    mktempdir() do root
        path = joinpath(root, "appended_cell_params.dat")
        write(path, """
Step lx ly lz alpha beta gamma vol
0 10.0 11.0 12.0 90.0 90.0 90.0 1320.0
100 10.1 11.1 12.1 90.0 90.0 90.0 1356.531
Step lx ly lz alpha beta gamma vol
200 10.2 11.2 12.2 90.0 90.0 90.0 1393.728
""")

        series = read_cell_parameter_series(path)
        @test length(series) == 3
        @test series.header_source == :file
        @test collect(series[:step]) == [0.0, 100.0, 200.0]
        @test collect(series[:volume]) == [1320.0, 1356.531, 1393.728]

        @test_throws CellParameterParseError read_cell_parameter_series(
            path;
            allow_repeated_headers=false,
        )
    end
end

@testset "Explicit and inferred cell-parameter layouts" begin
    mktempdir() do root
        explicit_path = joinpath(root, "explicit.dat")
        write(explicit_path, """
0.0 20.0 21.0 22.0 90.0 91.0 92.0 9200.0 300.0
1.0 20.1 21.1 22.1 90.0 91.0 92.0 9374.031 301.0
""")

        explicit = read_cell_parameter_series(
            explicit_path;
            columns=[
                :time,
                :a,
                :b,
                :c,
                :alpha,
                :beta,
                :gamma,
                :volume,
                :temperature,
            ],
        )
        @test explicit.header_source == :explicit
        @test collect(explicit[:time]) == [0.0, 1.0]
        @test collect(explicit[:temperature]) == [300.0, 301.0]

        core_path = joinpath(root, "core_only.dat")
        write(core_path, """
10.0 11.0 12.0 90.0 90.0 90.0 1320.0
10.1 11.1 12.1 90.0 90.0 90.0 1356.531
""")
        core = read_cell_parameter_series(core_path)
        @test core.header_source == :inferred
        @test cell_parameter_names(core) == collect(
            JZeoMCMD._CELL_PARAMETER_CORE_COLUMNS,
        )

        step_path = joinpath(root, "step_and_core.dat")
        write(step_path, """
0 10.0 11.0 12.0 90.0 90.0 90.0 1320.0
100 10.1 11.1 12.1 90.0 90.0 90.0 1356.531
""")
        step_series = read_cell_parameter_series(step_path)
        @test step_series.header_source == :inferred
        @test cell_parameter_names(step_series) == [
            :step,
            :a,
            :b,
            :c,
            :alpha,
            :beta,
            :gamma,
            :volume,
        ]

        temperature_path = joinpath(root, "core_and_temperature.dat")
        write(temperature_path, """
10.01 11.0 12.0 90.0 90.0 90.0 1321.32 299.9
10.02 11.1 12.1 90.0 90.0 90.0 1345.7862 300.1
""")
        temperature_series = read_cell_parameter_series(temperature_path)
        @test cell_parameter_names(temperature_series)[end] == :temperature
        @test collect(temperature_series[:temperature]) == [299.9, 300.1]
    end
end

@testset "Cell-parameter parser diagnostics" begin
    mktempdir() do root
        missing = joinpath(root, "missing.dat")
        @test_throws CellParameterParseError read_cell_parameter_series(missing)

        empty_path = joinpath(root, "empty.dat")
        write(empty_path, "# no data\n")
        @test_throws CellParameterParseError read_cell_parameter_series(empty_path)

        malformed = joinpath(root, "malformed.dat")
        write(malformed, """
# a b c alpha beta gamma volume
10 11 12 90 90 90 1320
10 11 wrong 90 90 90 1320
""")
        error = try
            read_cell_parameter_series(malformed)
            nothing
        catch err
            err
        end
        @test error isa CellParameterParseError
        @test error.line == 3
        @test occursin("non-numeric", error.message)

        inconsistent = joinpath(root, "inconsistent.dat")
        write(inconsistent, """
# a b c alpha beta gamma volume
10 11 12 90 90 90 1320
10 11 12 90 90 90
""")
        @test_throws CellParameterParseError read_cell_parameter_series(inconsistent)

        duplicate = joinpath(root, "duplicate.dat")
        write(duplicate, """
# Cella lx c alpha beta gamma volume
10 10 12 90 90 90 1200
""")
        @test_throws CellParameterParseError read_cell_parameter_series(duplicate)

        missing_required = joinpath(root, "missing_required.dat")
        write(missing_required, """
# a b c alpha beta gamma temperature
10 11 12 90 90 90 300
""")
        @test_throws CellParameterParseError read_cell_parameter_series(
            missing_required,
        )

        nonfinite = joinpath(root, "nonfinite.dat")
        write(nonfinite, """
# a b c alpha beta gamma volume temperature
10 11 12 90 90 90 1320 NaN
10.1 11.1 12.1 90 90 90 1356.531 300
""")
        @test_throws CellParameterParseError read_cell_parameter_series(nonfinite)
        accepted = read_cell_parameter_series(nonfinite; allow_nonfinite=true)
        @test isnan(accepted[:temperature][1])

        invalid_angle = joinpath(root, "invalid_angle.dat")
        write(invalid_angle, """
# a b c alpha beta gamma volume
10 11 12 0 90 90 1320
""")
        @test_throws CellParameterParseError read_cell_parameter_series(invalid_angle)
        unvalidated = read_cell_parameter_series(
            invalid_angle;
            validate_physical=false,
        )
        @test unvalidated[:alpha][1] == 0.0

        decreasing_step = joinpath(root, "decreasing_step.dat")
        write(decreasing_step, """
# step a b c alpha beta gamma volume
100 10 11 12 90 90 90 1320
50 10 11 12 90 90 90 1320
""")
        @test_throws CellParameterParseError read_cell_parameter_series(
            decreasing_step,
        )

        ambiguous = joinpath(root, "ambiguous.dat")
        write(ambiguous, "0 10 11 12 90 90 90 1320\n")
        @test_throws CellParameterParseError read_cell_parameter_series(ambiguous)
    end
end

@testset "Cell-parameter parser preserves additional observables" begin
    mktempdir() do root
        path = joinpath(root, "extra_columns.dat")
        write(path, """
# Step Cella Cellb Cellc CellAlpha CellBeta CellGamma Volume xy xz yz custom_metric
0 10 11 12 90 90 90 1320 0.1 0.2 0.3 7.5
100 10.1 11.1 12.1 90 90 90 1356.531 0.1 0.2 0.3 7.7
""")

        series = read_cell_parameter_series(path)
        @test has_cell_parameter(series, :xy)
        @test has_cell_parameter(series, :custom_metric)
        @test collect(series[:custom_metric]) == [7.5, 7.7]
    end
end
