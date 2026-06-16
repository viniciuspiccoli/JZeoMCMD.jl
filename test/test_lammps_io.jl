@testset "LAMMPS data I/O" begin
    @testset "orthorhombic cell matrix" begin
        data = synthetic_lammps_data()
        H = box_to_matrix(data)

        @test H == [
            10.0  0.0  0.0
             0.0 11.0  0.0
             0.0  0.0 12.0
        ]
        @test det(H) ≈ 1320.0
    end

    @testset "triclinic cell matrix" begin
        data = synthetic_lammps_data(triclinic=true)
        H = box_to_matrix(data)

        @test H == [
            10.0   0.0   0.0
             1.25 11.0   0.0
            -0.50  0.75 12.0
        ]
        @test det(H) ≈ 1320.0
    end

    @testset "full-style round trip" begin
        original = synthetic_lammps_data(triclinic=true)

        mktempdir() do dir
            path = joinpath(dir, "roundtrip.data")
            write_lammps_data(path, original; comment="JZeoMCMD test fixture")
            recovered = read_lammps_data(path; verbose=false)

            @test recovered.atom_style == :full
            @test recovered.masses == original.masses
            @test recovered.atom_ids == original.atom_ids
            @test recovered.molecule_labels == original.molecule_labels
            @test recovered.atom_labels == original.atom_labels
            @test recovered.atom_charges ≈ original.atom_charges
            @test recovered.coords ≈ original.coords
            @test recovered.box_dimensions ≈ original.box_dimensions
            @test recovered.tilt_factors ≈ original.tilt_factors
            @test recovered.image_flags == original.image_flags

            @test recovered.bonds == original.bonds
            @test recovered.bond_labels == original.bond_labels
            @test recovered.nbond_types == original.nbond_types
        end
    end
end
