@testset "CIF writing" begin
    data = synthetic_lammps_data(triclinic=true)

    mktempdir() do dir
        path = joinpath(dir, "framework.cif")
        write_cif(
            path,
            data;
            framework_types=[1, 2],
            type_elements=Dict(1 => "Si", 2 => "O"),
            comment="synthetic framework",
        )

        @test isfile(path)
        content = read(path, String)

        @test occursin("_cell_length_a", content)
        @test occursin("_cell_angle_alpha", content)
        @test occursin("_cell_volume", content)
        @test occursin("1320.00", content)
        @test occursin("Si1", content)
        @test occursin("O1", content)
        @test !occursin("C1", content)
    end
end
