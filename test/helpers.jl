"""
Create a small, fully populated `LammpsData` object for unit tests.
This script does not depend on RASPA, LAMMPS, OVITO, VMD, etc.
"""
function synthetic_lammps_data(; triclinic::Bool=false)
    data = LammpsData()

    data.masses = Dict(
        1 => 28.0855,
        2 => 15.9994,
        3 => 15.0350,
    )

    data.coords = [
        0.50  0.75  1.00
        2.00  2.50  3.00
        4.00  4.50  5.00
        5.00  4.50  5.00
    ]
    data.molecule_labels = [0, 0, 1, 1]
    data.atom_charges = [0.5236, -0.2618, 0.0, 0.0]
    data.atom_labels = [1, 2, 3, 3]
    data.atom_ids = [1, 2, 3, 4]

    data.box_dimensions = [
        0.0  10.0
        0.0  11.0
        0.0  12.0
    ]
    data.tilt_factors = triclinic ? [1.25, -0.50, 0.75] : [0.0, 0.0, 0.0]
    data.image_flags = [
         0  0  0
         0  0  0
         1  0 -1
         1  0 -1
    ]

    data.bonds = reshape([3, 4], 1, 2)
    data.bond_labels = [1]
    data.nbond_types = 1

    data.angles = zeros(Int, 0, 3)
    data.angle_labels = Int[]
    data.nangle_types = 0

    data.dihedrals = zeros(Int, 0, 4)
    data.dihedral_labels = Int[]
    data.ndihedral_types = 0

    data.impropers = zeros(Int, 0, 4)
    data.improper_labels = Int[]
    data.nimproper_types = 0

    data.atom_style = :full
    return data
end
