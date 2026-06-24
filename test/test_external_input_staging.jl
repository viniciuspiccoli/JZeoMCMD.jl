@testset "External input staging" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write(joinpath(root, "simulation.json.template"),
              "{\"Systems\": [], \"marker\": \"initial\"}")
        write(joinpath(root, "simulation.json.template_next"),
              "{\"Systems\": [], \"marker\": \"iterative\"}")

        cycle1 = joinpath(root, "work", "cycle_01")
        staged1 = stage_external_inputs(
            inputs,
            cycle1;
            cycle=1,
            base_dir=root,
        )

        @test staged1.cycle == 1
        @test staged1.cycle_dir == cycle1
        @test staged1.raspa_dir == joinpath(cycle1, "raspa")
        @test staged1.lammps_dir == joinpath(cycle1, "lammps")
        @test isdir(staged1.raspa_dir)
        @test isdir(staged1.lammps_dir)

        @test basename(staged1.framework_cif) == "framework.cif"
        @test basename(staged1.framework_data) == "framework.data"
        @test basename(staged1.raspa_simulation) ==
              "simulation.json.template"
        @test occursin("initial", read(staged1.raspa_simulation, String))
        @test basename(staged1.raspa_force_field) == "force_field.json"
        @test basename(only(staged1.raspa_molecule_files)) == "ethanol.json"
        @test basename(staged1.lammps_input) == "run_npt.in"
        @test basename(only(staged1.lammps_force_field_files)) ==
              "framework.ff"
        @test basename(only(staged1.lammps_auxiliary_files)) ==
              "framework.table"
        @test all(entry -> isfile(entry.destination), staged1.files)

        cycle2 = joinpath(root, "work", "cycle_02")
        staged2 = stage_external_inputs(
            inputs,
            cycle2;
            cycle=2,
            base_dir=root,
        )
        @test basename(staged2.raspa_simulation) ==
              "simulation.json.template_next"
        @test occursin("iterative", read(staged2.raspa_simulation, String))
    end
end

@testset "Current framework overrides" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write_text_file(joinpath(root, "current", "distorted.cif"),
                        "data_distorted\n")
        write_text_file(joinpath(root, "current", "npt_final.data"),
                        "relaxed framework\n")

        staged = stage_external_inputs(
            inputs,
            joinpath(root, "cycle_02");
            cycle=2,
            base_dir=root,
            framework_cif=joinpath("current", "distorted.cif"),
            framework_data=joinpath("current", "npt_final.data"),
        )

        @test basename(staged.framework_cif) == "distorted.cif"
        @test basename(staged.framework_data) == "npt_final.data"
        @test read(staged.framework_cif, String) == "data_distorted\n"
        @test read(staged.framework_data, String) == "relaxed framework\n"
    end
end

@testset "Staging overwrite policy" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        target = joinpath(root, "cycle_01")
        staged = stage_external_inputs(inputs, target; base_dir=root)

        err = try
            stage_external_inputs(inputs, target; base_dir=root)
            nothing
        catch caught
            caught
        end
        @test err isa ExternalInputStagingError
        @test occursin("already exists", sprint(showerror, err))

        write(joinpath(root, "run_npt.in"), "units metal\n")
        restaged = stage_external_inputs(
            inputs,
            target;
            base_dir=root,
            overwrite=true,
        )
        @test read(restaged.lammps_input, String) == "units metal\n"
        @test staged.lammps_input == restaged.lammps_input
    end
end

@testset "Staging preflight catches override collisions" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write_text_file(joinpath(root, "alternate", "force_field.json"),
                        "data_conflicting_name\n")
        target = joinpath(root, "cycle_01")

        err = try
            stage_external_inputs(
                inputs,
                target;
                base_dir=root,
                framework_cif=joinpath("alternate", "force_field.json"),
            )
            nothing
        catch caught
            caught
        end

        @test err isa ExternalInputStagingError
        @test occursin("already assigned", sprint(showerror, err))
        @test !isdir(joinpath(target, "raspa"))
        @test !isdir(joinpath(target, "lammps"))
    end
end

@testset "Staging argument and source checks" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)

        @test_throws ArgumentError stage_external_inputs(
            inputs,
            joinpath(root, "cycle_00");
            cycle=0,
            base_dir=root,
        )

        @test_throws ExternalInputStagingError stage_external_inputs(
            inputs,
            "";
            base_dir=root,
        )

        @test_throws ExternalInputStagingError stage_external_inputs(
            inputs,
            joinpath(root, "cycle_01");
            base_dir=root,
            framework_data="missing.data",
        )

        blocking_file = joinpath(root, "not_a_directory")
        write(blocking_file, "block")
        @test_throws ExternalInputStagingError stage_external_inputs(
            inputs,
            blocking_file;
            base_dir=root,
        )
    end
end
