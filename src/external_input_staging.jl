# Generic staging of externally supplied scientific input files.
#
# This layer copies validated inputs into cycle-specific RASPA and LAMMPS
# directories. It does not patch templates and does not execute either
# simulator.

"""
    StagedFile

One source-to-destination copy performed by [`stage_external_inputs`](@ref).
`role` identifies the purpose of the file in the coupled workflow.
"""
struct StagedFile
    role::Symbol
    source::String
    destination::String
end

"""
    StagedExternalInputs

Resolved paths to all files staged for one GCMC/NPT-MD cycle.

The selected RASPA simulation file is available as `raspa_simulation`. Cycle 1
uses `raspa_simulation_initial`; later cycles use
`raspa_simulation_iterative` when one was supplied. Filenames are preserved so
that include statements and component references remain valid.
"""
struct StagedExternalInputs
    cycle::Int
    cycle_dir::String
    raspa_dir::String
    lammps_dir::String

    framework_cif::String
    framework_data::String

    raspa_simulation::String
    raspa_force_field::String
    raspa_molecule_files::Vector{String}
    raspa_auxiliary_files::Vector{String}

    lammps_input::String
    lammps_force_field_files::Vector{String}
    lammps_auxiliary_files::Vector{String}

    files::Vector{StagedFile}
end

"""
    ExternalInputStagingError

Exception raised when external inputs cannot be staged safely. Typical causes
are a destination filename collision, an existing destination file when
`overwrite=false`, or a missing framework override.
"""
struct ExternalInputStagingError <: Exception
    message::String
end

Base.showerror(io::IO, err::ExternalInputStagingError) = print(io, err.message)

function _resolve_stage_path(path::AbstractString,
                             base_dir::AbstractString)::String
    stripped = strip(String(path))
    isempty(stripped) &&
        throw(ExternalInputStagingError("a staging source path cannot be empty"))
    expanded = expanduser(stripped)
    return normpath(abspath(isabspath(expanded) ? expanded :
                           joinpath(base_dir, expanded)))
end

function _require_staging_source(path::AbstractString, role::Symbol)::String
    source = String(path)
    ispath(source) ||
        throw(ExternalInputStagingError(
            "staging source for $(role) does not exist: $(source)",
        ))
    isfile(source) ||
        throw(ExternalInputStagingError(
            "staging source for $(role) is not a regular file: $(source)",
        ))
    return source
end

function _stage_entry(role::Symbol,
                      source::AbstractString,
                      destination_dir::AbstractString)::StagedFile
    src = _require_staging_source(source, role)
    return StagedFile(role, src, joinpath(destination_dir, basename(src)))
end

function _deduplicate_and_check_plan(entries::Vector{StagedFile})
    planned = StagedFile[]
    by_destination = Dict{String,StagedFile}()

    for entry in entries
        key = lowercase(normpath(entry.destination))
        if haskey(by_destination, key)
            previous = by_destination[key]
            if normpath(previous.source) == normpath(entry.source)
                # Preserve the additional role in the returned staging record.
                # The copy loop will copy this source/destination pair once.
                push!(planned, entry)
                continue
            end
            throw(ExternalInputStagingError(
                "cannot stage $(entry.role) from $(entry.source): destination " *
                "$(entry.destination) is already assigned to $(previous.role) " *
                "from $(previous.source)",
            ))
        end
        by_destination[key] = entry
        push!(planned, entry)
    end
    return planned
end

function _preflight_destinations(entries::Vector{StagedFile};
                                 overwrite::Bool)
    for entry in entries
        destination = entry.destination
        if ispath(destination)
            if normpath(abspath(destination)) == normpath(abspath(entry.source))
                continue
            elseif isdir(destination)
                throw(ExternalInputStagingError(
                    "staging destination is a directory: $(destination)",
                ))
            elseif !overwrite
                throw(ExternalInputStagingError(
                    "staging destination already exists: $(destination). " *
                    "Pass overwrite=true to replace it.",
                ))
            end
        end
    end
    return nothing
end

function _copy_stage_plan!(entries::Vector{StagedFile}; overwrite::Bool)
    copied_destinations = Set{String}()
    for entry in entries
        key = lowercase(normpath(entry.destination))
        key in copied_destinations && continue
        push!(copied_destinations, key)

        if normpath(abspath(entry.source)) ==
           normpath(abspath(entry.destination))
            continue
        end
        mkpath(dirname(entry.destination))
        cp(entry.source, entry.destination; force=overwrite)
    end
    return nothing
end

_paths_for_role(entries::Vector{StagedFile}, role::Symbol) =
    [entry.destination for entry in entries if entry.role === role]

function _path_for_role(entries::Vector{StagedFile}, role::Symbol)::String
    paths = _paths_for_role(entries, role)
    length(paths) == 1 ||
        throw(ExternalInputStagingError(
            "expected exactly one staged file for $(role), found $(length(paths))",
        ))
    return only(paths)
end

"""
    stage_external_inputs(inputs, cycle_directory; kwargs...)
        -> StagedExternalInputs

Copy one cycle's externally supplied scientific inputs into:

```text
cycle_directory/
├── raspa/
└── lammps/
```

The function performs no template substitution and runs no external program.
It only stages files and returns their exact destination paths.

Keyword arguments:

- `cycle=1`: workflow cycle number. Must be at least one.
- `base_dir=pwd()`: base directory used to resolve relative input paths,
  relative framework overrides, and a relative `cycle_directory`.
- `framework_cif=nothing`: optional current-cycle CIF. Defaults to
  `inputs.initial_cif`.
- `framework_data=nothing`: optional current-cycle LAMMPS data file. Defaults
  to `inputs.initial_data`.
- `validate=true`: run [`assert_valid_external_inputs`](@ref) before staging.
- `overwrite=false`: refuse to replace any existing destination file.
- `validation_kwargs...`: forwarded to `assert_valid_external_inputs`, for
  example `check_json=false` or `strict_extensions=true`.

Filenames are preserved. Distinct source files that map to the same filename
inside one simulator directory are rejected before any file is copied.
"""
function stage_external_inputs(
    inputs::ExternalInputFiles,
    cycle_directory::AbstractString;
    cycle::Integer=1,
    base_dir::AbstractString=pwd(),
    framework_cif::Union{Nothing,AbstractString}=nothing,
    framework_data::Union{Nothing,AbstractString}=nothing,
    validate::Bool=true,
    overwrite::Bool=false,
    validation_kwargs...,
)::StagedExternalInputs
    cycle >= 1 || throw(ArgumentError("cycle must be greater than or equal to 1"))

    base = normpath(abspath(expanduser(String(base_dir))))
    resolved = if validate
        assert_valid_external_inputs(inputs;
                                     base_dir=base,
                                     validation_kwargs...)
    else
        resolve_external_inputs(inputs; base_dir=base)
    end

    cycle_text = strip(String(cycle_directory))
    isempty(cycle_text) &&
        throw(ExternalInputStagingError("cycle_directory cannot be empty"))
    cycle_root = normpath(abspath(expanduser(
        isabspath(cycle_text) ? cycle_text : joinpath(base, cycle_text),
    )))
    if ispath(cycle_root) && !isdir(cycle_root)
        throw(ExternalInputStagingError(
            "cycle_directory exists and is not a directory: $(cycle_root)",
        ))
    end

    current_cif = isnothing(framework_cif) ? resolved.initial_cif :
                  _resolve_stage_path(framework_cif, base)
    current_data = isnothing(framework_data) ? resolved.initial_data :
                   _resolve_stage_path(framework_data, base)
    _require_staging_source(current_cif, :framework_cif)
    _require_staging_source(current_data, :framework_data)

    raspa_dir = joinpath(cycle_root, "raspa")
    lammps_dir = joinpath(cycle_root, "lammps")

    raw_plan = StagedFile[
        _stage_entry(:framework_cif, current_cif, raspa_dir),
        _stage_entry(:raspa_simulation,
                     raspa_simulation_for_cycle(resolved, cycle),
                     raspa_dir),
        _stage_entry(:raspa_force_field, resolved.raspa_force_field, raspa_dir),
    ]
    append!(raw_plan,
            (_stage_entry(:raspa_molecule_file, path, raspa_dir)
             for path in resolved.raspa_molecule_files))
    append!(raw_plan,
            (_stage_entry(:raspa_auxiliary_file, path, raspa_dir)
             for path in resolved.raspa_auxiliary_files))

    push!(raw_plan, _stage_entry(:framework_data, current_data, lammps_dir))
    push!(raw_plan, _stage_entry(:lammps_input, resolved.lammps_input,
                                 lammps_dir))
    append!(raw_plan,
            (_stage_entry(:lammps_force_field_file, path, lammps_dir)
             for path in resolved.lammps_force_field_files))
    append!(raw_plan,
            (_stage_entry(:lammps_auxiliary_file, path, lammps_dir)
             for path in resolved.lammps_auxiliary_files))

    plan = _deduplicate_and_check_plan(raw_plan)
    _preflight_destinations(plan; overwrite=overwrite)

    # Directory creation and copying happen only after the complete plan has
    # passed collision and overwrite checks.
    mkpath(raspa_dir)
    mkpath(lammps_dir)
    _copy_stage_plan!(plan; overwrite=overwrite)

    return StagedExternalInputs(
        Int(cycle),
        cycle_root,
        raspa_dir,
        lammps_dir,
        _path_for_role(plan, :framework_cif),
        _path_for_role(plan, :framework_data),
        _path_for_role(plan, :raspa_simulation),
        _path_for_role(plan, :raspa_force_field),
        _paths_for_role(plan, :raspa_molecule_file),
        _paths_for_role(plan, :raspa_auxiliary_file),
        _path_for_role(plan, :lammps_input),
        _paths_for_role(plan, :lammps_force_field_file),
        _paths_for_role(plan, :lammps_auxiliary_file),
        plan,
    )
end
