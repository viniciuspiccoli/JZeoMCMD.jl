@testset "Analysis utilities" begin
    @testset "strain" begin
        reference = (
            a=10.0,
            b=20.0,
            c=30.0,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            volume=6000.0,
        )
        current = (
            a=10.1,
            b=19.8,
            c=30.3,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            volume=6060.0,
        )

        strain = compute_strain(current, reference)
        @test strain.εa ≈ 0.01
        @test strain.εb ≈ -0.01
        @test strain.εc ≈ 0.01
        @test strain.εV ≈ 0.01
    end

    @testset "undefined strain is NaN" begin
        reference = (
            a=0.0,
            b=20.0,
            c=30.0,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            volume=6000.0,
        )
        current = (
            a=10.0,
            b=NaN,
            c=30.0,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            volume=6000.0,
        )

        strain = compute_strain(current, reference)
        @test isnan(strain.εa)
        @test isnan(strain.εb)
        @test strain.εc == 0.0
        @test strain.εV == 0.0
    end

    @testset "single thermo block" begin
        log_text = """
LAMMPS test log
Step Temp Press Cella Cellb Cellc CellAlpha CellBeta CellGamma Volume
0 300.0 1.0 10.0 11.0 12.0 90.0 90.0 90.0 1320.0
100 301.0 2.0 10.1 11.1 12.1 90.0 90.0 90.0 1356.531
Loop time of 1.0 on 1 procs for 100 steps
"""

        mktempdir() do dir
            path = joinpath(dir, "log.lammps")
            write(path, log_text)

            thermo = parse_lammps_log(path)
            @test thermo["Step"] == [0.0, 100.0]
            @test thermo["Temp"] == [300.0, 301.0]

            cell = extract_cell_params(thermo)
            @test cell.a ≈ 10.1
            @test cell.b ≈ 11.1
            @test cell.c ≈ 12.1
            @test cell.alpha ≈ 90.0
            @test cell.beta ≈ 90.0
            @test cell.gamma ≈ 90.0
            @test cell.volume ≈ 1356.531
        end
    end
end
