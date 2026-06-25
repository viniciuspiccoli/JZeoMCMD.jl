## script to running initial tests of the package

using JZeoMCMD
using Test
using LinearAlgebra
using JSON

include("helpers.jl")
include("test_lammps_io.jl")
include("test_analysis.jl")
include("test_statistics_types.jl")
include("test_cell_parameter_series.jl")
include("test_observable_summaries.jl")
include("test_cycle_statistics_reporting.jl")
include("test_cif.jl")
include("test_external_inputs.jl")

include("test_external_input_validation.jl")
include("test_external_input_staging.jl")

include("test_external_workflow.jl")
