# External-input bridge for the existing GCMC/NPT-MD workflow.
#
# This file deliberately reuses the current chemistry-specific construction,
# CIF-writing, analysis, and convergence functions from `workflow.jl`. Its
# purpose is to remove package-internal scientific input files from execution
# without changing the established silicalite-1 / H-ZSM-5 behavior yet.

const DEFAULT_EXTERNAL_LAMMPS_OUTPUT_FILES = (
    "npt_final.data",
    "loaded_npt_final.data",
    "npt_final.lmp",
    "loaded_npt_final.lmp",
)

"""
    PreparedExternalWorkflowCycle

Files prepared for one cycle of the external-input workflow.

`staged` contains the exact source-to-destination staging record.
`raspa_runtime_input` is the rendered RASPA input that will be read by RASPA3
(default filename: `simulation.json`).
"""
struct PreparedExternalWorkflowCycle
    staged::StagedExternalInputs
    raspa_runtime_input::String
end

"""
    ExternalWorkflowRunResult

Result returned by
`run_gcmc_md_workflow(wp::WorkflowParams, inputs::ExternalInputFiles)`.

`stop_reason` is one of `:converged`, `:max_iterations`, `:gcmc_failed`, or
`:npt_failed`.
"""
struct ExternalWorkflowRunResult
    history::Vector{Dict{String,Any}}
    converged::Bool
    stop_reason::Symbol
    final_cif::String
    final_data::String
end

"""
    ExternalWorkflowError

Exception raised while preparing or executing the external-input workflow.
"""
struct ExternalWorkflowError <: Exception
    message::String
end

Base.showerror(io::IO, err::ExternalWorkflowError) = print(io, err.message)

function _workflow_runtime_filename(name::AbstractString,
                                    role::AbstractString)::String
    value = strip(String(name))
    isempty(value) &&
        throw(ExternalWorkflowError("$(role) filename cannot be empty"))
    basename(value) == value ||
        throw(ExternalWorkflowError(
            "$(role) must be a filename without directory components: $(value)",
        ))
    return value
end

function _external_cycle_directory(wp::WorkflowParams, cycle::Integer)::String
    cycle >= 1 || throw(ArgumentError("cycle must be greater than or equal to 1"))
    root = normpath(abspath(expanduser(wp.base_dir)))
    return joinpath(root, @sprintf("cycle_%02d", cycle))
end

function _current_external_path(path::Union{Nothing,AbstractString},
                                fallback::AbstractString,
                                base_dir::AbstractString)::String
    isnothing(path) && return String(fallback)
    return _resolve_stage_path(path, base_dir)
end

function _raspa_runtime_sources(inputs::ExternalInputFiles,
                                cycle::Integer,
                                framework_cif::AbstractString)
    selected_simulation = raspa_simulation_for_cycle(inputs, cycle)
    sources = Tuple{Symbol,String}[
        (:framework_cif, String(framework_cif)),
        (:raspa_force_field, inputs.raspa_force_field),
    ]
    append!(sources,
            ((:raspa_molecule_file, path)
             for path in inputs.raspa_molecule_files))
    append!(sources,
            ((:raspa_auxiliary_file, path)
             for path in inputs.raspa_auxiliary_files))
    return selected_simulation, sources
end

function _check_raspa_runtime_collision(inputs::ExternalInputFiles,
                                        cycle::Integer,
                                        framework_cif::AbstractString,
                                        runtime_filename::AbstractString)
    selected, sources = _raspa_runtime_sources(inputs, cycle, framework_cif)
    runtime_key = lowercase(String(runtime_filename))

    # The selected simulation may itself already be named `simulation.json`;
    # in that case it is intentionally rendered in place after staging.
    for (role, source) in sources
        if lowercase(basename(source)) == runtime_key &&
           normpath(source) != normpath(selected)
            throw(ExternalWorkflowError(
                "cannot create RASPA runtime input $(runtime_filename): " *
                "the staged $(role) file $(source) has the same filename",
            ))
        end
    end
    return nothing
end

function _replace_external_workflow_placeholders(text::AbstractString,
                                                 wp::WorkflowParams,
                                                 framework_name::AbstractString)
    rendered = String(text)
    substitutions = (
        "__FRAMEWORK__"   => String(framework_name),
        "__TEMPERATURE__" => string(wp.temperature),
        "__PRESSURE__"    => string(wp.pressure),
        "__N_INIT__"      => string(wp.raspa_n_init),
        "__N_EQUIL__"     => string(wp.raspa_n_equil),
        "__N_PROD__"      => string(wp.raspa_n_prod),
        "__PRINT_EVERY__" => string(wp.raspa_print_every),
    )
    for (placeholder, value) in substitutions
        rendered = replace(rendered, placeholder => value)
    end

    unresolved = unique([m.match for m in eachmatch(r"__[A-Za-z0-9_]+__", rendered)])
    isempty(unresolved) ||
        throw(ExternalWorkflowError(
            "unresolved placeholders in RASPA input: $(join(unresolved, ", "))",
        ))
    return rendered
end

function _patch_external_raspa_dictionary!(document::AbstractDict,
                                           wp::WorkflowParams,
                                           framework_name::AbstractString,
                                           system_index::Integer)
    haskey(document, "Systems") ||
        throw(ExternalWorkflowError(
            "RASPA simulation input does not contain a Systems array",
        ))
    systems = document["Systems"]
    systems isa AbstractVector ||
        throw(ExternalWorkflowError("RASPA Systems entry is not an array"))
    1 <= system_index <= length(systems) ||
        throw(ExternalWorkflowError(
            "RASPA system_index=$(system_index) is outside the available " *
            "range 1:$(length(systems))",
        ))
    system = systems[system_index]
    system isa AbstractDict ||
        throw(ExternalWorkflowError(
            "RASPA Systems[$(system_index)] is not a JSON object",
        ))

    system["Name"] = String(framework_name)
    system["ExternalTemperature"] = wp.temperature
    system["ExternalPressure"] = wp.pressure

    # These fields are made authoritative in the external-input entry point.
    # This corrects the legacy situation where fixed numeric values in a JSON
    # template silently overrode WorkflowParams.
    document["NumberOfInitializationCycles"] = wp.raspa_n_init
    document["NumberOfEquilibrationCycles"] = wp.raspa_n_equil
    document["NumberOfCycles"] = wp.raspa_n_prod
    document["PrintEvery"] = wp.raspa_print_every
    return document
end

function _render_external_raspa_input!(source::AbstractString,
                                       destination::AbstractString,
                                       wp::WorkflowParams,
                                       framework_cif::AbstractString;
                                       system_index::Integer=1)
    isfile(source) ||
        throw(ExternalWorkflowError(
            "staged RASPA simulation input does not exist: $(source)",
        ))

    framework_name = splitext(basename(framework_cif))[1]
    text = read(source, String)
    text = _replace_external_workflow_placeholders(text, wp, framework_name)

    document = try
        JSON.parse(text)
    catch err
        throw(ExternalWorkflowError(
            "failed to parse rendered RASPA input $(source): " *
            sprint(showerror, err),
        ))
    end
    document isa AbstractDict ||
        throw(ExternalWorkflowError(
            "RASPA simulation input must contain a top-level JSON object",
        ))

    _patch_external_raspa_dictionary!(document, wp, framework_name,
                                      system_index)
    mkpath(dirname(destination))
    open(destination, "w") do io
        JSON.print(io, document, 2)
        println(io)
    end
    return String(destination)
end

"""
    prepare_external_workflow_cycle(wp, inputs, cycle; kwargs...)
        -> PreparedExternalWorkflowCycle

Validate and stage externally supplied scientific inputs for one workflow
cycle, then render the selected RASPA simulation file as `simulation.json`.
No external program is executed.

This function is useful for inspecting the exact files that will be supplied to
RASPA3 and LAMMPS before launching a long calculation.

Keyword arguments:

- `input_base_dir=pwd()`: base directory for relative paths in `inputs` and
  relative current-framework overrides. This is independent of `wp.base_dir`,
  which remains the output directory.
- `current_cif=nothing`: current framework CIF; defaults to
  `inputs.initial_cif`.
- `current_data=nothing`: current LAMMPS framework data; defaults to
  `inputs.initial_data`.
- `runtime_raspa_filename="simulation.json"`: filename read by RASPA3.
- `raspa_system_index=1`: one-based entry in the RASPA `Systems` array to
  update.
- `validate=true`: validate the complete external input specification.
- `overwrite=false`: allow replacement of an already staged cycle.
- `validation_kwargs...`: forwarded to `assert_valid_external_inputs`.
"""
function prepare_external_workflow_cycle(
    wp::WorkflowParams,
    inputs::ExternalInputFiles,
    cycle::Integer;
    input_base_dir::AbstractString=pwd(),
    current_cif::Union{Nothing,AbstractString}=nothing,
    current_data::Union{Nothing,AbstractString}=nothing,
    runtime_raspa_filename::AbstractString="simulation.json",
    raspa_system_index::Integer=1,
    validate::Bool=true,
    overwrite::Bool=false,
    validation_kwargs...,
)::PreparedExternalWorkflowCycle
    cycle >= 1 || throw(ArgumentError("cycle must be greater than or equal to 1"))
    raspa_system_index >= 1 ||
        throw(ArgumentError("raspa_system_index must be at least one"))

    input_root = normpath(abspath(expanduser(String(input_base_dir))))
    resolved = if validate
        assert_valid_external_inputs(inputs;
                                     base_dir=input_root,
                                     validation_kwargs...)
    else
        resolve_external_inputs(inputs; base_dir=input_root)
    end

    current_cif_path = _current_external_path(current_cif,
                                              resolved.initial_cif,
                                              input_root)
    current_data_path = _current_external_path(current_data,
                                               resolved.initial_data,
                                               input_root)
    runtime_name = _workflow_runtime_filename(runtime_raspa_filename,
                                              "RASPA runtime input")
    _check_raspa_runtime_collision(resolved, cycle, current_cif_path,
                                   runtime_name)

    cycle_root = _external_cycle_directory(wp, cycle)
    runtime_path = joinpath(cycle_root, "raspa", runtime_name)
    if ispath(runtime_path) && !overwrite
        selected_source = raspa_simulation_for_cycle(resolved, cycle)
        selected_destination = joinpath(cycle_root, "raspa",
                                        basename(selected_source))
        if normpath(runtime_path) != normpath(selected_destination)
            throw(ExternalWorkflowError(
                "RASPA runtime input already exists: $(runtime_path). " *
                "Pass overwrite=true to replace it.",
            ))
        end
    end

    staged = stage_external_inputs(
        resolved,
        cycle_root;
        cycle=cycle,
        base_dir=input_root,
        framework_cif=current_cif_path,
        framework_data=current_data_path,
        validate=false,
        overwrite=overwrite,
    )

    runtime_path = joinpath(staged.raspa_dir, runtime_name)
    _render_external_raspa_input!(staged.raspa_simulation,
                                  runtime_path,
                                  wp,
                                  staged.framework_cif;
                                  system_index=raspa_system_index)

    # An explicit overwrite request means that this cycle is being rebuilt.
    # Remove stale RASPA output so a previous restart cannot be mistaken for
    # the result of the new execution.
    raspa_output = joinpath(staged.raspa_dir, "output")
    overwrite && isdir(raspa_output) && rm(raspa_output; recursive=true)

    return PreparedExternalWorkflowCycle(staged, runtime_path)
end

function _external_command(
    command::AbstractString,
    extra_arguments::AbstractVector{<:AbstractString}=String[];
    dir::AbstractString,
)
    # `split` returns `Vector{SubString{String}}`. Concatenating that vector
    # with `Vector{String}` promotes the result to `Vector{AbstractString}`,
    # but Julia's `Cmd` constructor requires a concrete `Vector{String}`.
    words = String.(split(strip(String(command))))
    isempty(words) &&
        throw(ExternalWorkflowError("external executable command is empty"))

    arguments = vcat(words, String.(extra_arguments))

    # The `dir` keyword is supported by the Cmd-to-Cmd constructor, not by
    # the Vector-to-Cmd constructor. Construct the command first, then attach
    # its working directory.
    return Cmd(Cmd(arguments); dir=String(dir))
end

function _external_restart_loading(restart_path::AbstractString,
                                   component_name::AbstractString,
                                   atoms_per_molecule::Integer)::Float64
    atoms_per_molecule > 0 ||
        throw(ExternalWorkflowError("atoms_per_mol must be greater than zero"))

    restart = JSON.parsefile(restart_path)
    sites = nothing

    if haskey(restart, component_name)
        sites = restart[component_name]
    else
        candidates = [(String(key), value) for (key, value) in restart
                      if value isa AbstractVector]
        if length(candidates) == 1
            fallback_name, sites = only(candidates)
            @warn "RASPA restart component not found; using the only vector-valued entry" requested=component_name fallback=fallback_name
        else
            @warn "RASPA restart component not found or ambiguous" requested=component_name vector_entries=first.(candidates)
            return NaN
        end
    end

    sites isa AbstractVector || begin
        @warn "RASPA restart component is not an array" component=component_name
        return NaN
    end
    isempty(sites) && return 0.0

    nsites = length(sites)
    if nsites % atoms_per_molecule != 0
        @warn "RASPA restart site count is not divisible by atoms_per_mol" component=component_name nsites=nsites atoms_per_molecule=atoms_per_molecule
        return NaN
    end
    return Float64(nsites ÷ atoms_per_molecule)
end

function _run_external_gcmc!(wp::WorkflowParams,
                             prepared::PreparedExternalWorkflowCycle)
    rdir = prepared.staged.raspa_dir
    println("  [GCMC] Running RASPA3 (P=$(wp.pressure) Pa)...")
    t0 = time()
    try
        # RASPA3 must receive the rendered runtime input explicitly.
        # `raspa_exe` should contain only the executable/launcher command,
        # for example `raspa3` or `mpirun -np 4 raspa3`.
        input_name = basename(prepared.raspa_runtime_input)
        command = _external_command(
            wp.raspa_exe,
            [input_name];
            dir=rdir,
        )
        run(command)
    catch err
        @warn "RASPA3 failed" exception=(err, catch_backtrace())
        return nothing
    end
    @printf("  [GCMC] Done (%.1f s)\n", time() - t0)

    output_dir = joinpath(rdir, "output")
    if !isdir(output_dir)
        @warn "No output/ directory in staged RASPA directory" directory=rdir
        return nothing
    end

    restarts = sort(filter(name -> startswith(name, "restart_") &&
                                  endswith(name, ".json"),
                           readdir(output_dir)))
    if isempty(restarts)
        @warn "No restart JSON found" directory=output_dir
        return nothing
    end
    restart_path = joinpath(output_dir, last(restarts))

    n_ads = try
        _external_restart_loading(restart_path,
                                  wp.adsorbate_component,
                                  wp.atoms_per_mol)
    catch err
        @warn "Failed to parse RASPA restart loading" exception=(err, catch_backtrace())
        NaN
    end
    @printf("  [GCMC] N_ads = %.0f molecules\n", n_ads)

    return (
        restart=restart_path,
        output=Dict{String,Any}("n_ads" => n_ads),
        dir=rdir,
    )
end

function _copy_external_runtime_data(source::AbstractString,
                                     destination::AbstractString)
    isfile(source) ||
        throw(ExternalWorkflowError(
            "LAMMPS runtime data source does not exist: $(source)",
        ))
    if normpath(abspath(source)) != normpath(abspath(destination))
        cp(source, destination; force=true)
    end
    return String(destination)
end

function _external_output_candidates(candidates)
    candidates isa AbstractString && return (String(candidates),)
    return candidates
end

function _clear_external_lammps_outputs!(ldir::AbstractString,
                                         candidates,
                                         protected_paths::Vector{String})
    protected = Set(normpath(abspath(path)) for path in protected_paths)
    for candidate in _external_output_candidates(candidates)
        filename = _workflow_runtime_filename(String(candidate),
                                              "LAMMPS output candidate")
        path = joinpath(ldir, filename)
        normalized = normpath(abspath(path))
        normalized in protected &&
            throw(ExternalWorkflowError(
                "LAMMPS output candidate $(filename) conflicts with a runtime input",
            ))
        isfile(path) && rm(path; force=true)
    end
    return nothing
end

function _external_lammps_output(ldir::AbstractString,
                                 candidates)
    for candidate in _external_output_candidates(candidates)
        filename = _workflow_runtime_filename(String(candidate),
                                              "LAMMPS output candidate")
        path = joinpath(ldir, filename)
        isfile(path) && return path
    end
    return nothing
end

function _run_external_npt!(
    wp::WorkflowParams,
    prepared::PreparedExternalWorkflowCycle,
    data_file::AbstractString;
    lammps_data_filename::AbstractString="loaded.lmp",
    lammps_output_files=DEFAULT_EXTERNAL_LAMMPS_OUTPUT_FILES,
)
    ldir = prepared.staged.lammps_dir
    data_name = _workflow_runtime_filename(lammps_data_filename,
                                           "LAMMPS runtime data")
    runtime_data = joinpath(ldir, data_name)
    _copy_external_runtime_data(data_file, runtime_data)

    input_name = basename(prepared.staged.lammps_input)
    input_path = joinpath(ldir, input_name)
    _clear_external_lammps_outputs!(
        ldir,
        lammps_output_files,
        [runtime_data, input_path],
    )

    println("  [NPT] Running LAMMPS...")
    t0 = time()
    try
        command = _external_command(
            wp.lammps_exe,
            ["-in", input_name];
            dir=ldir,
        )
        run(command)
    catch err
        @warn "LAMMPS failed" exception=(err, catch_backtrace())
        return nothing
    end
    @printf("  [NPT] Done (%.1f s)\n", time() - t0)

    output_data = _external_lammps_output(ldir, lammps_output_files)
    if isnothing(output_data)
        @warn "No configured LAMMPS output file found" directory=ldir candidates=collect(_external_output_candidates(lammps_output_files))
        return nothing
    end
    println("  [NPT] Output: $(basename(output_data))")

    log_file = joinpath(ldir, "log.lammps")
    log_data = isfile(log_file) ?
        parse_lammps_log(log_file) : Dict{String,Vector{Float64}}()
    cell = extract_cell_params(log_data)
    @printf("  [NPT] a=%.4f b=%.4f c=%.4f V=%.1f\n",
            cell.a, cell.b, cell.c, cell.volume)

    return (data=output_data, cell=cell, log=log_data, dir=ldir)
end

"""
    run_gcmc_md_workflow(wp::WorkflowParams, inputs::ExternalInputFiles; kwargs...)
        -> ExternalWorkflowRunResult

Run the existing coupled GCMC/NPT-MD algorithm using only user-supplied
scientific input files. The legacy one-argument method remains unchanged.

This is a migration bridge: file acquisition and simulator execution are
externalized, while framework topology construction, adsorbate rebuilding,
CIF atom-type mapping, MFI channel analysis, and the current convergence test
still use the existing silicalite-1 / H-ZSM-5 logic. Those scientific mappings
will be externalized in subsequent steps.

Important keywords:

- `input_base_dir=pwd()`: resolves relative paths in `inputs` independently of
  the output directory `wp.base_dir`.
- `overwrite_staged=false`: refuses to overwrite existing cycle inputs.
- `runtime_raspa_filename="simulation.json"`: runtime RASPA input filename.
- `raspa_system_index=1`: RASPA system entry updated with T, P, and framework.
- `lammps_data_filename="loaded.lmp"`: runtime data filename expected by the
  current LAMMPS scripts.
- `lammps_output_files`: ordered filenames accepted as the final NPT data.
- `validation_kwargs...`: forwarded to initial external-input validation.
"""
function run_gcmc_md_workflow(
    wp::WorkflowParams,
    inputs::ExternalInputFiles;
    input_base_dir::AbstractString=pwd(),
    overwrite_staged::Bool=false,
    runtime_raspa_filename::AbstractString="simulation.json",
    raspa_system_index::Integer=1,
    lammps_data_filename::AbstractString="loaded.lmp",
    lammps_output_files=DEFAULT_EXTERNAL_LAMMPS_OUTPUT_FILES,
    validation_kwargs...,
)::ExternalWorkflowRunResult
    wp.max_iterations >= 1 ||
        throw(ArgumentError("max_iterations must be greater than or equal to 1"))
    wp.atoms_per_mol >= 1 ||
        throw(ArgumentError("atoms_per_mol must be greater than or equal to 1"))

    input_root = normpath(abspath(expanduser(String(input_base_dir))))
    resolved = assert_valid_external_inputs(inputs;
                                            base_dir=input_root,
                                            validation_kwargs...)

    # Work on a copy so the caller's legacy WorkflowParams object is not
    # mutated by path normalization or by the external-input bridge.
    run_params = deepcopy(wp)
    run_params.base_dir = normpath(abspath(expanduser(wp.base_dir)))
    run_params.initial_cif = resolved.initial_cif
    run_params.initial_data = resolved.initial_data
    mkpath(run_params.base_dir)

    println("╔═══════════════════════════════════════════════════╗")
    println("║  JZeoMCMD — External-input GCMC/NPT-MD Workflow  ║")
    println("╚═══════════════════════════════════════════════════╝")
    @printf("  T = %.1f K   P = %.1e Pa\n",
            run_params.temperature, run_params.pressure)
    println("  Workdir:      $(run_params.base_dir)")
    println("  Initial CIF:  $(resolved.initial_cif)")
    println("  Initial data: $(resolved.initial_data)")
    println("  RASPA input:  $(resolved.raspa_simulation_initial)")
    println("  LAMMPS input: $(resolved.lammps_input)")
    println("  Internal scientific input folders: not used")
    println()

    current_cif = resolved.initial_cif
    current_data = resolved.initial_data
    history = Dict{String,Any}[]
    converged = false
    stop_reason = :max_iterations

    # Preserve the current reference-cell calculation in this migration step.
    # Triclinic reference-cell handling is scheduled as a separate tested
    # change so that file externalization does not alter numerical behavior.
    reference = read_lammps_data(current_data; verbose=false)
    bounds = reference.box_dimensions
    ref_cell = (
        a=bounds[1,2] - bounds[1,1],
        b=bounds[2,2] - bounds[2,1],
        c=bounds[3,2] - bounds[3,1],
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        volume=prod(bounds[d,2] - bounds[d,1] for d in 1:3),
    )
    @printf("  Reference cell: a=%.3f b=%.3f c=%.3f V=%.1f\n\n",
            ref_cell.a, ref_cell.b, ref_cell.c, ref_cell.volume)

    for cycle in 1:run_params.max_iterations
        println("═" ^ 60)
        @printf("  CYCLE %d / %d\n", cycle, run_params.max_iterations)
        println("═" ^ 60)

        prepared = prepare_external_workflow_cycle(
            run_params,
            resolved,
            cycle;
            input_base_dir=input_root,
            current_cif=current_cif,
            current_data=current_data,
            runtime_raspa_filename=runtime_raspa_filename,
            raspa_system_index=raspa_system_index,
            validate=false,
            overwrite=overwrite_staged,
        )

        println("\n  ── Step 1: RASPA3 GCMC ──")
        println("    CIF: $(basename(prepared.staged.framework_cif))")
        gcmc = _run_external_gcmc!(run_params, prepared)
        if isnothing(gcmc)
            println("  ✗ GCMC failed — stopping.")
            stop_reason = :gcmc_failed
            break
        end

        println("\n  ── Step 2: Build LAMMPS data file ──")
        build_params = deepcopy(run_params)
        build_params.initial_data = prepared.staged.framework_data
        data_file = step_build_data!(
            build_params,
            cycle,
            gcmc.restart,
            prepared.staged.framework_data,
        )

        println("\n  ── Step 3: LAMMPS NPT-MD ──")
        npt = _run_external_npt!(
            run_params,
            prepared,
            data_file;
            lammps_data_filename=lammps_data_filename,
            lammps_output_files=lammps_output_files,
        )
        if isnothing(npt)
            println("  ✗ NPT failed — stopping.")
            stop_reason = :npt_failed
            break
        end

        println("\n  ── Step 4: Extract distorted CIF ──")
        new_cif = step_write_cif!(run_params, cycle, npt.data)

        println("\n  ── Step 5: Analysis ──")
        summary = step_analyze!(run_params, cycle, npt.cell, gcmc.output,
                                ref_cell, npt.data)
        push!(history, summary)

        current_cif = new_cif
        current_data = npt.data

        if check_convergence(history, run_params)
            println("\n★ Converged after $(cycle) iterations ★")
            converged = true
            stop_reason = :converged
            break
        end
        println()
    end

    println("\n" * "═" ^ 60)
    println("  WORKFLOW COMPLETE: $(length(history)) completed cycles")
    println("  Stop reason: $(stop_reason)")
    println("  Results: $(joinpath(run_params.base_dir, "convergence.csv"))")
    println("═" ^ 60)

    return ExternalWorkflowRunResult(
        history,
        converged,
        stop_reason,
        current_cif,
        current_data,
    )
end
