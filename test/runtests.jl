## script to running initial tests of the package

using JZeoMCMD
using Test
using LinearAlgebra

include("helpers.jl")
include("test_lammps_io.jl")
include("test_analysis.jl")
include("test_cif.jl")
include("test_external_inputs.jl")
