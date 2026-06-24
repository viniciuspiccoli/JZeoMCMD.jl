# External scientific-input validation for `ExternalInputFiles`.
#
# This layer is intentionally separate from the input specification and from
# file staging. It performs no simulation and does not modify user files.


"""
    ExternalInputValidationIssue

One problem found while validating an [`ExternalInputFiles`](@ref)
specification.

`severity` is either `:error` or `:warning`. Errors prevent staging and
execution. Warnings identify suspicious but technically usable inputs, such as
an unconventional filename extension.
"""
struct ExternalInputValidationIssue
    severity::Symbol
    code::Symbol
    field::Symbol
    path::String
    message::String

    function ExternalInputValidationIssue(severity::Symbol,
                                          code::Symbol,
                                          field::Symbol,
                                          path::AbstractString,
                                          message::AbstractString)
        severity in (:error, :warning) ||
            throw(ArgumentError("severity must be :error or :warning"))
        return new(severity, code, field, String(path), String(message))
    end
end

"""
    ExternalInputValidationReport

Result returned by [`validate_external_inputs`](@ref).

`resolved_inputs` contains normalized absolute paths. `issues` contains all
validation errors and warnings detected in one pass, allowing the user to fix
multiple input problems at once.
"""
struct ExternalInputValidationReport
    resolved_inputs::ExternalInputFiles
    issues::Vector{ExternalInputValidationIssue}
end

"""
    ExternalInputValidationError

Exception thrown by [`assert_valid_external_inputs`](@ref) when one or more
validation errors are present.
"""
struct ExternalInputValidationError <: Exception
    report::ExternalInputValidationReport
end

"""
    validation_errors(report)

Return only the error-level issues in a validation report.
"""
validation_errors(report::ExternalInputValidationReport) =
    filter(issue -> issue.severity === :error, report.issues)

"""
    validation_warnings(report)

Return only the warning-level issues in a validation report.
"""
validation_warnings(report::ExternalInputValidationReport) =
    filter(issue -> issue.severity === :warning, report.issues)

"""
    external_inputs_valid(report) -> Bool

Return `true` when the report contains no error-level issues. Warnings do not
make an input specification invalid.
"""
external_inputs_valid(report::ExternalInputValidationReport) =
    isempty(validation_errors(report))

function Base.showerror(io::IO, err::ExternalInputValidationError)
    errors = validation_errors(err.report)
    print(io, "external input validation failed with ", length(errors),
          length(errors) == 1 ? " error" : " errors")
    for issue in errors
        print(io, "\n  - [", issue.code, "] ", issue.field)
        isempty(issue.path) || print(io, " (", issue.path, ")")
        print(io, ": ", issue.message)
    end
end

# Expand the common Unix `~` shorthand without depending on shell expansion.
function _expand_user_path(path::AbstractString)::String
    p = String(path)
    p == "~" && return homedir()
    p == "~/" && return homedir()
    startswith(p, "~/") && return joinpath(homedir(), p[3:end])
    return p
end

function _resolve_external_path(path::AbstractString,
                                base_dir::AbstractString)::String
    p = strip(String(path))
    isempty(p) && return p
    p = _expand_user_path(p)
    return normpath(isabspath(p) ? p : joinpath(base_dir, p))
end

"""
    resolve_external_inputs(inputs; base_dir=pwd()) -> ExternalInputFiles

Return a copy of `inputs` in which every nonempty path is normalized and made
absolute relative to `base_dir`. No filesystem or file-content checks are
performed.
"""
function resolve_external_inputs(inputs::ExternalInputFiles;
                                 base_dir::AbstractString=pwd())::ExternalInputFiles
    root = abspath(_expand_user_path(String(base_dir)))
    resolve_one(path) = _resolve_external_path(path, root)

    return ExternalInputFiles(
        initial_cif=resolve_one(inputs.initial_cif),
        initial_data=resolve_one(inputs.initial_data),
        raspa_simulation_initial=resolve_one(inputs.raspa_simulation_initial),
        raspa_simulation_iterative=isnothing(inputs.raspa_simulation_iterative) ?
            nothing : resolve_one(inputs.raspa_simulation_iterative),
        raspa_force_field=resolve_one(inputs.raspa_force_field),
        raspa_molecule_files=resolve_one.(inputs.raspa_molecule_files),
        raspa_auxiliary_files=resolve_one.(inputs.raspa_auxiliary_files),
        lammps_input=resolve_one(inputs.lammps_input),
        lammps_force_field_files=resolve_one.(inputs.lammps_force_field_files),
        lammps_auxiliary_files=resolve_one.(inputs.lammps_auxiliary_files),
    )
end

# Flatten the input structure while retaining the originating field name.
function _external_file_entries(inputs::ExternalInputFiles)
    entries = Tuple{Symbol,String}[
        (:initial_cif, inputs.initial_cif),
        (:initial_data, inputs.initial_data),
        (:raspa_simulation_initial, inputs.raspa_simulation_initial),
        (:raspa_force_field, inputs.raspa_force_field),
        (:lammps_input, inputs.lammps_input),
    ]

    if !isnothing(inputs.raspa_simulation_iterative)
        push!(entries,
              (:raspa_simulation_iterative,
               inputs.raspa_simulation_iterative))
    end

    append!(entries,
            ((:raspa_molecule_files, path)
             for path in inputs.raspa_molecule_files))
    append!(entries,
            ((:raspa_auxiliary_files, path)
             for path in inputs.raspa_auxiliary_files))
    append!(entries,
            ((:lammps_force_field_files, path)
             for path in inputs.lammps_force_field_files))
    append!(entries,
            ((:lammps_auxiliary_files, path)
             for path in inputs.lammps_auxiliary_files))

    return entries
end

function _extension_is_expected(field::Symbol, path::AbstractString)::Bool
    name = lowercase(basename(path))

    field === :initial_cif && return endswith(name, ".cif")
    field === :initial_data &&
        return any(endswith(name, suffix) for suffix in (".data", ".lmp", ".dat"))

    field in (:raspa_simulation_initial,
              :raspa_simulation_iterative,
              :raspa_force_field,
              :raspa_molecule_files) && return occursin(".json", name)

    field === :lammps_input &&
        return any(endswith(name, suffix) for suffix in (".in", ".lmp", ".lammps"))

    field === :lammps_force_field_files &&
        return any(endswith(name, suffix)
                   for suffix in (".ff", ".inc", ".in", ".lmp", ".lammps"))

    # Auxiliary files are intentionally unrestricted.
    return true
end

function _expected_extension_description(field::Symbol)::String
    field === :initial_cif && return ".cif"
    field === :initial_data && return ".data, .lmp, or .dat"
    field in (:raspa_simulation_initial,
              :raspa_simulation_iterative,
              :raspa_force_field,
              :raspa_molecule_files) && return "a filename containing .json"
    field === :lammps_input && return ".in, .lmp, or .lammps"
    field === :lammps_force_field_files &&
        return ".ff, .inc, .in, .lmp, or .lammps"
    return "an unrestricted extension"
end

# RASPA simulation templates may contain unquoted placeholders such as
# __TEMPERATURE__. Replacing placeholder tokens with zero validates the JSON
# structure without assigning physical values.
function _parse_json_or_template(path::AbstractString;
                                 allow_placeholders::Bool=false)
    text = read(path, String)
    if allow_placeholders
        text = replace(text, r"__[A-Za-z0-9_]+__" => "0")
    end
    return JSON.parse(text)
end

function _push_staging_collision!(issues::Vector{ExternalInputValidationIssue},
                                  destination::Symbol,
                                  field::Symbol,
                                  path::String,
                                  previous_field::Symbol,
                                  previous_path::String)
    push!(issues,
          ExternalInputValidationIssue(
              :error,
              :duplicate_destination_name,
              field,
              path,
              "$(basename(path)) would overwrite the file supplied by " *
              "$(previous_field) ($(previous_path)) in the $(destination) " *
              "staging directory",
          ))
end

function _check_staging_collisions!(
    issues::Vector{ExternalInputValidationIssue},
    destination::Symbol,
    entries::Vector{Tuple{Symbol,String}},
)
    seen = Dict{String,Tuple{Symbol,String}}()

    for (field, path) in entries
        isempty(path) && continue
        key = lowercase(basename(path))
        if haskey(seen, key)
            previous_field, previous_path = seen[key]
            # Reusing the exact same source file is intentional and safe.
            if normpath(path) != normpath(previous_path)
                _push_staging_collision!(issues, destination, field, path,
                                         previous_field, previous_path)
            end
        else
            seen[key] = (field, path)
        end
    end
    return issues
end

function _raspa_common_entries(inputs::ExternalInputFiles)
    entries = Tuple{Symbol,String}[
        (:initial_cif, inputs.initial_cif),
        (:raspa_force_field, inputs.raspa_force_field),
    ]
    append!(entries,
            ((:raspa_molecule_files, path)
             for path in inputs.raspa_molecule_files))
    append!(entries,
            ((:raspa_auxiliary_files, path)
             for path in inputs.raspa_auxiliary_files))
    return entries
end

function _check_entry_against_group!(
    issues::Vector{ExternalInputValidationIssue},
    destination::Symbol,
    entry::Tuple{Symbol,String},
    group::Vector{Tuple{Symbol,String}},
)
    field, path = entry
    isempty(path) && return issues
    key = lowercase(basename(path))

    for (previous_field, previous_path) in group
        isempty(previous_path) && continue
        if lowercase(basename(previous_path)) == key &&
           normpath(previous_path) != normpath(path)
            _push_staging_collision!(issues, destination, field, path,
                                     previous_field, previous_path)
            break
        end
    end
    return issues
end

function _check_all_staging_collisions!(
    issues::Vector{ExternalInputValidationIssue},
    inputs::ExternalInputFiles,
)
    raspa_common = _raspa_common_entries(inputs)
    _check_staging_collisions!(issues, :raspa, raspa_common)

    # Initial and iterative simulation templates are staged in different
    # cycles. Check each against the common RASPA files, but not against one
    # another.
    _check_entry_against_group!(
        issues,
        :raspa,
        (:raspa_simulation_initial, inputs.raspa_simulation_initial),
        raspa_common,
    )
    if !isnothing(inputs.raspa_simulation_iterative)
        _check_entry_against_group!(
            issues,
            :raspa,
            (:raspa_simulation_iterative,
             inputs.raspa_simulation_iterative),
            raspa_common,
        )
    end

    lammps_entries = Tuple{Symbol,String}[
        (:initial_data, inputs.initial_data),
        (:lammps_input, inputs.lammps_input),
    ]
    append!(lammps_entries,
            ((:lammps_force_field_files, path)
             for path in inputs.lammps_force_field_files))
    append!(lammps_entries,
            ((:lammps_auxiliary_files, path)
             for path in inputs.lammps_auxiliary_files))
    _check_staging_collisions!(issues, :lammps, lammps_entries)

    return issues
end

"""
    validate_external_inputs(inputs; kwargs...) -> ExternalInputValidationReport

Validate an external scientific-input specification and collect all detected
problems in one report.

Keyword arguments:

- `base_dir=pwd()`: base directory for resolving relative paths.
- `check_json=true`: parse RASPA JSON files. Known `__PLACEHOLDER__`
  tokens are accepted in simulation templates.
- `check_extensions=true`: check conventional filename extensions.
- `strict_extensions=false`: treat unconventional extensions as errors rather
  than warnings.
- `check_collisions=true`: detect distinct files that would have the same
  basename in a RASPA or LAMMPS staging directory.

The returned report always contains normalized absolute paths in
`resolved_inputs`. Missing files, directories supplied where files are
required, invalid JSON, empty paths, and staging collisions are errors.
Filename extensions are warnings by default because external simulation files
may legitimately use local naming conventions.
"""
function validate_external_inputs(
    inputs::ExternalInputFiles;
    base_dir::AbstractString=pwd(),
    check_json::Bool=true,
    check_extensions::Bool=true,
    strict_extensions::Bool=false,
    check_collisions::Bool=true,
)::ExternalInputValidationReport
    resolved = resolve_external_inputs(inputs; base_dir=base_dir)
    issues = ExternalInputValidationIssue[]

    for (field, path) in _external_file_entries(resolved)
        if isempty(strip(path))
            push!(issues,
                  ExternalInputValidationIssue(
                      :error,
                      :empty_path,
                      field,
                      path,
                      "an input path cannot be empty",
                  ))
            continue
        end

        if !ispath(path)
            push!(issues,
                  ExternalInputValidationIssue(
                      :error,
                      :missing_file,
                      field,
                      path,
                      "file does not exist",
                  ))
            continue
        elseif !isfile(path)
            push!(issues,
                  ExternalInputValidationIssue(
                      :error,
                      :not_a_file,
                      field,
                      path,
                      "expected a regular file",
                  ))
            continue
        end

        if check_extensions && !_extension_is_expected(field, path)
            severity = strict_extensions ? :error : :warning
            push!(issues,
                  ExternalInputValidationIssue(
                      severity,
                      :unexpected_extension,
                      field,
                      path,
                      "expected $(_expected_extension_description(field))",
                  ))
        end
    end

    if check_collisions
        _check_all_staging_collisions!(issues, resolved)
    end

    if check_json
        json_entries = Tuple{Symbol,String,Bool}[
            (:raspa_simulation_initial,
             resolved.raspa_simulation_initial,
             true),
            (:raspa_force_field, resolved.raspa_force_field, false),
        ]
        if !isnothing(resolved.raspa_simulation_iterative)
            push!(json_entries,
                  (:raspa_simulation_iterative,
                   resolved.raspa_simulation_iterative,
                   true))
        end
        append!(json_entries,
                ((:raspa_molecule_files, path, false)
                 for path in resolved.raspa_molecule_files))
        append!(json_entries,
                ((:raspa_auxiliary_files, path, false)
                 for path in resolved.raspa_auxiliary_files
                 if occursin(".json", lowercase(basename(path)))))

        for (field, path, allow_placeholders) in json_entries
            # Existence/type errors have already been reported above.
            isfile(path) || continue
            try
                _parse_json_or_template(path;
                                        allow_placeholders=allow_placeholders)
            catch err
                push!(issues,
                      ExternalInputValidationIssue(
                          :error,
                          :invalid_json,
                          field,
                          path,
                          sprint(showerror, err),
                      ))
            end
        end
    end

    return ExternalInputValidationReport(resolved, issues)
end

"""
    assert_valid_external_inputs(inputs; kwargs...) -> ExternalInputFiles

Validate `inputs`, throw [`ExternalInputValidationError`](@ref) when any error
is found, and otherwise return the resolved absolute-path specification.
Warnings do not prevent the function from returning.
"""
function assert_valid_external_inputs(inputs::ExternalInputFiles; kwargs...)
    report = validate_external_inputs(inputs; kwargs...)
    external_inputs_valid(report) || throw(ExternalInputValidationError(report))
    return report.resolved_inputs
end
