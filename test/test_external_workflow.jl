@testset "External workflow cycle preparation" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write(
            joinpath(root, "simulation.json.template"),
            """{
              \"NumberOfInitializationCycles\": 1,
              \"NumberOfEquilibrationCycles\": 2,
              \"NumberOfCycles\": 3,
              \"PrintEvery\": 4,
              \"Systems\": [{
                \"Name\": \"__FRAMEWORK__\",
                \"ExternalTemperature\": __TEMPERATURE__,
                \"ExternalPressure\": __PRESSURE__
              }]
            }""",
        )

        wp = WorkflowParams(
            base_dir=joinpath(root, "work"),
            temperature=373.0,
            pressure=25000.0,
            raspa_n_init=111,
            raspa_n_equil=222,
            raspa_n_prod=333,
            raspa_print_every=44,
        )

        prepared = prepare_external_workflow_cycle(
            wp,
            inputs,
            1;
            input_base_dir=root,
        )

        @test prepared isa PreparedExternalWorkflowCycle
        @test prepared.staged.cycle == 1
        @test prepared.raspa_runtime_input ==
              joinpath(root, "work", "cycle_01", "raspa", "simulation.json")
        @test isfile(prepared.raspa_runtime_input)
        @test isfile(prepared.staged.raspa_simulation)
        @test prepared.staged.raspa_simulation != prepared.raspa_runtime_input

        document = JSON.parsefile(prepared.raspa_runtime_input)
        system = only(document["Systems"])
        @test system["Name"] == "framework"
        @test system["ExternalTemperature"] == 373.0
        @test system["ExternalPressure"] == 25000.0
        @test document["NumberOfInitializationCycles"] == 111
        @test document["NumberOfEquilibrationCycles"] == 222
        @test document["NumberOfCycles"] == 333
        @test document["PrintEvery"] == 44

        # Preparation must not mutate the caller's input specification.
        @test inputs.initial_cif == "framework.cif"
        @test inputs.initial_data == "framework.data"
    end
end

@testset "External workflow patches plain JSON" begin
    mktempdir() do root
        inputs0 = make_valid_external_inputs(root)
        plain = joinpath(root, "simulation_initial.json")
        write(
            plain,
            """{
              \"NumberOfInitializationCycles\": 10,
              \"NumberOfEquilibrationCycles\": 10,
              \"NumberOfCycles\": 10,
              \"PrintEvery\": 10,
              \"Systems\": [{
                \"Name\": \"old_framework\",
                \"ExternalTemperature\": 100.0,
                \"ExternalPressure\": 1.0
              }]
            }""",
        )

        inputs = ExternalInputFiles(
            initial_cif=inputs0.initial_cif,
            initial_data=inputs0.initial_data,
            raspa_simulation_initial="simulation_initial.json",
            raspa_force_field=inputs0.raspa_force_field,
            raspa_molecule_files=inputs0.raspa_molecule_files,
            lammps_input=inputs0.lammps_input,
            lammps_force_field_files=inputs0.lammps_force_field_files,
            lammps_auxiliary_files=inputs0.lammps_auxiliary_files,
        )
        wp = WorkflowParams(
            base_dir=joinpath(root, "run"),
            temperature=423.0,
            pressure=1.5e5,
            raspa_n_init=20,
            raspa_n_equil=30,
            raspa_n_prod=40,
            raspa_print_every=5,
        )

        prepared = prepare_external_workflow_cycle(
            wp,
            inputs,
            1;
            input_base_dir=root,
        )
        document = JSON.parsefile(prepared.raspa_runtime_input)
        system = only(document["Systems"])

        @test system["Name"] == "framework"
        @test system["ExternalTemperature"] == 423.0
        @test system["ExternalPressure"] == 1.5e5
        @test document["NumberOfInitializationCycles"] == 20
        @test document["NumberOfEquilibrationCycles"] == 30
        @test document["NumberOfCycles"] == 40
        @test document["PrintEvery"] == 5

        # The staged source remains available for provenance.
        @test JSON.parsefile(prepared.staged.raspa_simulation)["NumberOfCycles"] == 10
    end
end

@testset "External workflow selects iterative RASPA input" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write(
            joinpath(root, "simulation.json.template_next"),
            """{
              \"marker\": \"iterative\",
              \"Systems\": [{
                \"Name\": \"__FRAMEWORK__\",
                \"ExternalTemperature\": __TEMPERATURE__,
                \"ExternalPressure\": __PRESSURE__
              }]
            }""",
        )
        write_text_file(joinpath(root, "current", "distorted.cif"),
                        "data_current\n")
        write_text_file(joinpath(root, "current", "npt_final.data"),
                        "current data\n")

        wp = WorkflowParams(base_dir=joinpath(root, "work"))
        prepared = prepare_external_workflow_cycle(
            wp,
            inputs,
            2;
            input_base_dir=root,
            current_cif=joinpath("current", "distorted.cif"),
            current_data=joinpath("current", "npt_final.data"),
        )

        document = JSON.parsefile(prepared.raspa_runtime_input)
        @test document["marker"] == "iterative"
        @test only(document["Systems"])["Name"] == "distorted"
        @test basename(prepared.staged.framework_cif) == "distorted.cif"
        @test basename(prepared.staged.framework_data) == "npt_final.data"
    end
end

@testset "External workflow runtime collision preflight" begin
    mktempdir() do root
        inputs0 = make_valid_external_inputs(root)
        write_text_file(joinpath(root, "simulation.json"),
                        "{\"unrelated\": true}\n")
        inputs = ExternalInputFiles(
            initial_cif=inputs0.initial_cif,
            initial_data=inputs0.initial_data,
            raspa_simulation_initial=inputs0.raspa_simulation_initial,
            raspa_simulation_iterative=inputs0.raspa_simulation_iterative,
            raspa_force_field=inputs0.raspa_force_field,
            raspa_molecule_files=inputs0.raspa_molecule_files,
            raspa_auxiliary_files=["simulation.json"],
            lammps_input=inputs0.lammps_input,
            lammps_force_field_files=inputs0.lammps_force_field_files,
            lammps_auxiliary_files=inputs0.lammps_auxiliary_files,
        )
        wp = WorkflowParams(base_dir=joinpath(root, "work"))

        err = try
            prepare_external_workflow_cycle(
                wp,
                inputs,
                1;
                input_base_dir=root,
            )
            nothing
        catch caught
            caught
        end
        @test err isa ExternalWorkflowError
        @test occursin("same filename", sprint(showerror, err))
        @test !isdir(joinpath(root, "work", "cycle_01"))
    end
end

@testset "External workflow API and argument checks" begin
    @test hasmethod(run_gcmc_md_workflow,
                    Tuple{WorkflowParams,ExternalInputFiles})

    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        wp = WorkflowParams(base_dir=joinpath(root, "work"))

        @test_throws ArgumentError prepare_external_workflow_cycle(
            wp, inputs, 0; input_base_dir=root,
        )
        @test_throws ArgumentError prepare_external_workflow_cycle(
            wp, inputs, 1;
            input_base_dir=root,
            raspa_system_index=0,
        )
        @test_throws ExternalWorkflowError prepare_external_workflow_cycle(
            wp, inputs, 1;
            input_base_dir=root,
            runtime_raspa_filename="inputs/simulation.json",
        )
    end
end

@testset "External command construction" begin
    mktempdir() do root
        raspa = JZeoMCMD._external_command(
            "raspa3",
            ["simulation.json"];
            dir=root,
        )
        @test raspa.exec == ["raspa3", "simulation.json"]
        @test raspa.dir == root
        @test eltype(raspa.exec) === String

        lammps = JZeoMCMD._external_command(
            "mpirun -np 4 lmp_mpi",
            ["-in", "run_npt.in"];
            dir=root,
        )
        @test lammps.exec == [
            "mpirun", "-np", "4", "lmp_mpi", "-in", "run_npt.in",
        ]
        @test lammps.dir == root
        @test eltype(lammps.exec) === String

        @test_throws ExternalWorkflowError JZeoMCMD._external_command(
            "   ",
            String[];
            dir=root,
        )
    end
end
