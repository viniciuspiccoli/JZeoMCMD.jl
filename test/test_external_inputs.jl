@testset "External input-file specification" begin
    inputs = ExternalInputFiles(
        initial_cif="framework.cif",
        initial_data="framework.data",
        raspa_simulation_initial="simulation_initial.json",
        raspa_simulation_iterative="simulation_iterative.json",
        raspa_force_field="force_field.json",
        raspa_molecule_files=["ethanol.json", "water.json"],
        raspa_auxiliary_files=["charges.json"],
        lammps_input="run_npt.in",
        lammps_force_field_files=["framework.ff", "guest.ff"],
        lammps_auxiliary_files=["framework.table"],
    )

    @test inputs.initial_cif == "framework.cif"
    @test inputs.initial_data == "framework.data"
    @test inputs.raspa_force_field == "force_field.json"
    @test inputs.raspa_molecule_files == ["ethanol.json", "water.json"]
    @test inputs.lammps_force_field_files == ["framework.ff", "guest.ff"]

    @test raspa_simulation_for_cycle(inputs, 1) == "simulation_initial.json"
    @test raspa_simulation_for_cycle(inputs, 2) == "simulation_iterative.json"
    @test raspa_simulation_for_cycle(inputs, 50) == "simulation_iterative.json"

    fallback = ExternalInputFiles(
        initial_cif="framework.cif",
        initial_data="framework.data",
        raspa_simulation_initial="simulation.json",
        raspa_force_field="force_field.json",
        lammps_input="run_npt.in",
    )

    @test isnothing(fallback.raspa_simulation_iterative)
    @test isempty(fallback.raspa_molecule_files)
    @test isempty(fallback.raspa_auxiliary_files)
    @test isempty(fallback.lammps_force_field_files)
    @test isempty(fallback.lammps_auxiliary_files)
    @test raspa_simulation_for_cycle(fallback, 1) == "simulation.json"
    @test raspa_simulation_for_cycle(fallback, 2) == "simulation.json"

    @test_throws ArgumentError raspa_simulation_for_cycle(inputs, 0)
    @test_throws ArgumentError raspa_simulation_for_cycle(inputs, -1)
end
