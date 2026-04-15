# JZeoMCMD

[![Build Status](https://github.com/viniciuspiccoli/JZeoMCMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/viniciuspiccoli/JZeoMCMD.jl/actions/workflows/CI.yml?query=branch%3Amain)


### check parameters and conversions - add original papers for the force fields
### check analysis and meaning of them
### check scripts that handle the transformations

### add different Si/Al structure with 2x2x2 supercell and their data file from ovito
### need to add aluminossilicate to reload_adsorbate.jl



┌─────────────────────────────────────────────────┐
│  Cycle 1 (initial)                              │
│                                                 │
│  test.data (Ovito)                              │
│       ↓  build_loaded_zeolite.jl                │
│  loaded_zeolite.lmp + run_loaded.in             │
│       ↓  LAMMPS NPT-MD                         │
│  loaded_npt_final.lmp                           │
│       ↓  write_cif.jl                           │
│  relaxed_MFI.cif                                │
└────────────────┬────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────┐
│  Cycle N (iterative)                            │
│                                                 │
│  relaxed_MFI.cif → RASPA3 GCMC → restart.json  │
│       ↓                              ↓          │
│       │    reload_adsorbate.jl ←─────┘          │
│       │    (reads loaded_npt_final.lmp          │
│       │     + new restart.json)                 │
│       ↓                                         │
│  cycleN_loaded.lmp                              │
│       ↓  LAMMPS NPT-MD (same run_loaded.in)    │
│  cycleN_npt_final.lmp                           │
│       ↓  write_cif.jl                           │
│  cycleN_relaxed.cif → RASPA3 GCMC → ...        │
└─────────────────────────────────────────────────┘




# what to do:

check: force_field.json
check charges and lennard jones (check why using lj parameters using Bai.). Why Si for silicate and Si for aluminossilicate are different?
check params.toml and params_loader.jl
adjust the write_lammps_data functions in: reload_adsorbate.jl, build_loaded_zeolite.jl, add_zeolite_topology.jl
check the paramters in generate_nb_tables.jl

add this: nohup ./your_script.sh > my_output.log 2>&1 &




base_dir/                         ← user provides these:
├── MFI_SI.data                   ← supercell data (any name)
├── MFI_SI.cif                    ← unit cell CIF (any name)
├── hillsauer_nb.table            ← NB table
├── run_npt.in                    ← LAMMPS input
└── cycle_01/...                  ← workflow creates these

the package provides (from raspa_inputs/):
force_field.json                  ← H-S charges + Bai LJ + TraPPE
ethanol.json                      ← TraPPE-UA molecule definition
simulation.json.template          ← cycle 1 (2×2×2 unit cells)
simulation.json.template_next     ← cycle 2+ (1×1×1 supercell)


wp = WorkflowParams(...)
wp.initial_cif  = "MFI_SI.cif"      # ← in base_dir
wp.initial_data = "MFI_SI.data"      # ← in base_dir
run_gcmc_md_workflow(wp)