function write_text_file(path::AbstractString, text::AbstractString="")
    mkpath(dirname(path))
    write(path, text)
    return path
end

function make_valid_external_inputs(root::AbstractString)
    write_text_file(joinpath(root, "framework.cif"), "data_framework\n")
    write_text_file(joinpath(root, "framework.data"), "LAMMPS data\n")

    # Bare numeric placeholders intentionally reproduce the current RASPA
    # template convention and are accepted by template-aware validation.
    write_text_file(
        joinpath(root, "simulation.json.template"),
        """{
          \"Systems\": [{
            \"Name\": \"__FRAMEWORK__\",
            \"ExternalTemperature\": __TEMPERATURE__,
            \"ExternalPressure\": __PRESSURE__
          }]
        }""",
    )
    write_text_file(
        joinpath(root, "simulation.json.template_next"),
        """{
          \"Systems\": [{
            \"Name\": \"__FRAMEWORK__\",
            \"ExternalTemperature\": __TEMPERATURE__,
            \"ExternalPressure\": __PRESSURE__
          }]
        }""",
    )
    write_text_file(joinpath(root, "force_field.json"), "{\"PseudoAtoms\": []}")
    write_text_file(joinpath(root, "ethanol.json"), "{\"Type\": \"flexible\"}")
    write_text_file(joinpath(root, "run_npt.in"), "units real\n")
    write_text_file(joinpath(root, "framework.ff"), "# force field\n")
    write_text_file(joinpath(root, "framework.table"), "# table\n")

    return ExternalInputFiles(
        initial_cif="framework.cif",
        initial_data="framework.data",
        raspa_simulation_initial="simulation.json.template",
        raspa_simulation_iterative="simulation.json.template_next",
        raspa_force_field="force_field.json",
        raspa_molecule_files=["ethanol.json"],
        lammps_input="run_npt.in",
        lammps_force_field_files=["framework.ff"],
        lammps_auxiliary_files=["framework.table"],
    )
end

@testset "External input validation" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)

        resolved = resolve_external_inputs(inputs; base_dir=root)
        @test isabspath(resolved.initial_cif)
        @test resolved.initial_cif == joinpath(root, "framework.cif")
        @test resolved.raspa_molecule_files == [joinpath(root, "ethanol.json")]
        @test resolved.raspa_simulation_iterative ==
              joinpath(root, "simulation.json.template_next")

        report = validate_external_inputs(inputs; base_dir=root)
        @test external_inputs_valid(report)
        @test isempty(validation_errors(report))
        @test isempty(validation_warnings(report))
        @test report.resolved_inputs.initial_cif == resolved.initial_cif
        @test report.resolved_inputs.lammps_input == resolved.lammps_input

        checked = assert_valid_external_inputs(inputs; base_dir=root)
        @test checked.initial_cif == resolved.initial_cif
        @test checked.raspa_force_field == resolved.raspa_force_field
    end
end

@testset "External input validation reports all errors" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        rm(joinpath(root, "framework.data"))
        write_text_file(joinpath(root, "force_field.json"), "{invalid json")
        mkpath(joinpath(root, "not_a_file.in"))

        broken = ExternalInputFiles(
            initial_cif="",
            initial_data="framework.data",
            raspa_simulation_initial="simulation.json.template",
            raspa_force_field="force_field.json",
            lammps_input="not_a_file.in",
        )

        report = validate_external_inputs(broken; base_dir=root)
        errors = validation_errors(report)
        codes = Set(issue.code for issue in errors)

        @test !external_inputs_valid(report)
        @test :empty_path in codes
        @test :missing_file in codes
        @test :not_a_file in codes
        @test :invalid_json in codes
        @test length(errors) == 4

        err = try
            assert_valid_external_inputs(broken; base_dir=root)
            nothing
        catch caught
            caught
        end
        @test err isa ExternalInputValidationError
        message = sprint(showerror, err)
        @test occursin("4 errors", message)
        @test occursin("missing_file", message)
        @test occursin("invalid_json", message)
    end
end

@testset "Extension policy" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        mv(joinpath(root, "framework.cif"), joinpath(root, "framework.structure"))

        unconventional = ExternalInputFiles(
            initial_cif="framework.structure",
            initial_data=inputs.initial_data,
            raspa_simulation_initial=inputs.raspa_simulation_initial,
            raspa_simulation_iterative=inputs.raspa_simulation_iterative,
            raspa_force_field=inputs.raspa_force_field,
            raspa_molecule_files=inputs.raspa_molecule_files,
            lammps_input=inputs.lammps_input,
            lammps_force_field_files=inputs.lammps_force_field_files,
            lammps_auxiliary_files=inputs.lammps_auxiliary_files,
        )

        report = validate_external_inputs(unconventional; base_dir=root)
        @test external_inputs_valid(report)
        @test length(validation_warnings(report)) == 1
        @test only(validation_warnings(report)).code == :unexpected_extension

        strict_report = validate_external_inputs(
            unconventional;
            base_dir=root,
            strict_extensions=true,
        )
        @test !external_inputs_valid(strict_report)
        @test only(validation_errors(strict_report)).code == :unexpected_extension

        unchecked = validate_external_inputs(
            unconventional;
            base_dir=root,
            check_extensions=false,
        )
        @test isempty(unchecked.issues)
    end
end

@testset "JSON checks can be disabled" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)
        write_text_file(joinpath(root, "ethanol.json"), "not JSON")

        checked = validate_external_inputs(inputs; base_dir=root)
        @test !external_inputs_valid(checked)
        @test any(issue -> issue.code == :invalid_json,
                  validation_errors(checked))

        unchecked = validate_external_inputs(
            inputs;
            base_dir=root,
            check_json=false,
        )
        @test external_inputs_valid(unchecked)
    end
end

@testset "Staging-name collisions" begin
    mktempdir() do root
        inputs = make_valid_external_inputs(root)

        dir_a = joinpath(root, "a")
        dir_b = joinpath(root, "b")
        write_text_file(joinpath(dir_a, "parameters.ff"), "A\n")
        write_text_file(joinpath(dir_b, "parameters.ff"), "B\n")

        collision = ExternalInputFiles(
            initial_cif=inputs.initial_cif,
            initial_data=inputs.initial_data,
            raspa_simulation_initial=inputs.raspa_simulation_initial,
            raspa_simulation_iterative=inputs.raspa_simulation_iterative,
            raspa_force_field=inputs.raspa_force_field,
            raspa_molecule_files=inputs.raspa_molecule_files,
            lammps_input=inputs.lammps_input,
            lammps_force_field_files=[
                joinpath("a", "parameters.ff"),
                joinpath("b", "parameters.ff"),
            ],
            lammps_auxiliary_files=inputs.lammps_auxiliary_files,
        )

        report = validate_external_inputs(collision; base_dir=root)
        @test !external_inputs_valid(report)
        collisions = filter(issue -> issue.code == :duplicate_destination_name,
                            validation_errors(report))
        @test length(collisions) == 1
        @test occursin("LAMMPS", uppercase(collisions[1].message))

        no_collision_check = validate_external_inputs(
            collision;
            base_dir=root,
            check_collisions=false,
        )
        @test external_inputs_valid(no_collision_check)

        # Referencing the exact same file twice is not considered an overwrite.
        reused = ExternalInputFiles(
            initial_cif=inputs.initial_cif,
            initial_data=inputs.initial_data,
            raspa_simulation_initial=inputs.raspa_simulation_initial,
            raspa_force_field=inputs.raspa_force_field,
            lammps_input=inputs.lammps_input,
            lammps_force_field_files=["framework.ff", "framework.ff"],
        )
        @test external_inputs_valid(
            validate_external_inputs(reused; base_dir=root),
        )
    end
end
