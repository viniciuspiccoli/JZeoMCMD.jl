# alumino calculation

using JZeoMCMD

inputs = ExternalInputFiles(
    initial_cif="mfi_sial_19.cif",
    initial_data="mfi_sial_19.data",

    raspa_simulation_initial="simulation.json.template",
    raspa_simulation_iterative="simulation.json.template_next",
    raspa_force_field="force_field.json",
    raspa_molecule_files=["ethanol.json"],

    lammps_input="run_npt.in",
    lammps_force_field_files=["hillsauer_alumino.ff"],
    lammps_auxiliary_files=["hillsauer_alumino.table"],
)

wp = WorkflowParams(
    base_dir="/home/viniciusp/.julia/dev/JZeoMCMD/run_simulation_alumino/results",
    temperature=373.0,
    pressure=1.0e5,
    raspa_n_init=10_000,
    raspa_n_equil=10_000,
    raspa_n_prod=50_000,
    raspa_print_every=1_000,
    max_iterations=15,
)



result = run_gcmc_md_workflow(
    wp,
    inputs;
    input_base_dir="/home/viniciusp/.julia/dev/JZeoMCMD/run_simulation_alumino",
)

println(result.stop_reason)
println(result.converged)
println(result.final_cif)
println(result.final_data)
