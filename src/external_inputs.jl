"""
    ExternalInputFiles

Paths to the scientific input files supplied by the user for a coupled
GCMC/NPT-MD calculation.

This type only stores the input-file specification. It deliberately does not
check whether the files exist and does not copy files into a workflow
directory. Validation and staging are separate operations so they can be
introduced and tested without changing the current simulation workflow.

Required files:

- `initial_cif`: initial framework CIF used by RASPA.
- `initial_data`: initial, topologized LAMMPS framework data file.
- `raspa_simulation_initial`: RASPA simulation JSON or template for cycle 1.
- `raspa_force_field`: RASPA force-field JSON.
- `lammps_input`: LAMMPS input script used for NPT-MD.

Optional collections allow each force field to use as many supporting files as
needed. `raspa_simulation_iterative` can point to a separate RASPA input for
cycles 2 and later. If it is `nothing`, the cycle-1 input is reused.

All paths may be absolute or relative. Their interpretation is deferred to the
validation/staging layer.
"""
Base.@kwdef struct ExternalInputFiles
    # Framework structures
    initial_cif::String
    initial_data::String

    # RASPA inputs
    raspa_simulation_initial::String
    raspa_force_field::String

    # LAMMPS inputs
    lammps_input::String

    # Optional RASPA inputs
    raspa_simulation_iterative::Union{Nothing,String} = nothing
    raspa_molecule_files::Vector{String} = String[]
    raspa_auxiliary_files::Vector{String} = String[]

    # Optional LAMMPS inputs
    lammps_force_field_files::Vector{String} = String[]
    lammps_auxiliary_files::Vector{String} = String[]
end

"""
    raspa_simulation_for_cycle(inputs, cycle) -> String

Return the RASPA simulation file assigned to `cycle`.

Cycle 1 always uses `raspa_simulation_initial`. Later cycles use
`raspa_simulation_iterative` when supplied; otherwise they reuse the initial
file.
"""
function raspa_simulation_for_cycle(inputs::ExternalInputFiles,
                                    cycle::Integer)::String
    cycle >= 1 || throw(ArgumentError("cycle must be greater than or equal to 1"))

    if cycle == 1 || isnothing(inputs.raspa_simulation_iterative)
        return inputs.raspa_simulation_initial
    end
    return inputs.raspa_simulation_iterative
end
