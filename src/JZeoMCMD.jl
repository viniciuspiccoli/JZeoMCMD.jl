module JZeoMCMD
	# pack of functions to parse lammps data file
	using Printf, LinearAlgebra, Statistics, JSON
	import TOML
	include("read_lammps_data.jl")
	export LammpsData, read_lammps_data, write_lammps_data
	export RaspaForceField, read_raspa3_forcefield, map_raspa3_forcefield!
	export make_supercell, box_to_matrix
	export compute_rdf, compute_msd, read_lammps_dump, parse_raspa3_isotherm
       
	# function to add full bonded topology to a lammps data file from OVITO
	# IT ONLY WORKS FOR OVITO data files - THE CODE CANNOT CREATE A DATA FILE FROM THE CIF!!!!
	include("add_zeolite_topology.jl")
	export add_zeolite_topology! 

	# function to build lammps data file for a zeolite loaded with adsorbate
	include("build_loaded_zeolite.jl")
	export ZeoliteConfig
	export add_framework_topology!, read_and_remap_framework
	export read_raspa3_ethanol, merge_framework_ethanol!, write_complete_data
	export write_input_script

	# write a cif file using the lammps data file from the previous NPT run
	include("write_cif.jl")
	export write_cif, write_framework_cif

	# read params.toml and populate ZeoliteConfig — single source of truth for FF parameters
	include("params_loader.jl")
	export load_config, load_ff_params, write_complete_data_toml

	# auto-generate Hill-Sauer A/r^9 non-bonded tables from params.toml
	include("generate_nb_tables.jl")
	export generate_table, compute_A_cross

	# LAMMPS log parser, strain tracker, RASPA3 output parser, MFI channel occupancy
	include("analysis.jl")
	export parse_lammps_log, extract_cell_params, compute_strain
	export parse_raspa3_output, analyze_channel_occupancy, write_cycle_summary

	# reload adsorbate: strip old ethanol from NPT data, insert new from RASPA3 JSON
	# uses write_complete_data from build_loaded_zeolite.jl
	include("reload_adsorbate.jl")
	export reload_adsorbate

	# master GCMC/MD iterative workflow + pressure sweep setup
	# all simulation files go to a user-chosen directory, NOT inside this package
	include("workflow.jl")
	include("setup_pressure_sweep.jl")
	export WorkflowParams, load_workflow_params, run_gcmc_md_workflow
	export setup_pressure_sweep, ensure_tables!

	# ── Package resource paths ──────────────────────────────────
	# Use these to find default files shipped with the package:
	#   JZeoMCMD.resource_dir()       → .../JZeoMCMD/src/
	#   JZeoMCMD.params_path()        → .../JZeoMCMD/src/ff/params.toml
	#   JZeoMCMD.raspa_inputs_path()  → .../JZeoMCMD/src/raspa_inputs/
	#   JZeoMCMD.structures_path()    → .../JZeoMCMD/src/structures_data_files/
	"""Directory containing package source files."""
	resource_dir() = @__DIR__
	"""Path to the default params.toml shipped with the package."""
	params_path() = joinpath(@__DIR__, "ff", "params.toml")
	"""Path to the default RASPA3 input files."""
	raspa_inputs_path() = joinpath(@__DIR__, "raspa_inputs")
	"""Path to the default structure files (CIF, Ovito data)."""
	structures_path() = joinpath(@__DIR__, "structures_data_files")
	export resource_dir, params_path, raspa_inputs_path, structures_path
end






