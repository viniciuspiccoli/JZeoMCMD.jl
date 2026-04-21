using Printf

include("read_lammps_data.jl")
include("add_zeolite_topology.jl")
include("refine_topology_types.jl")

# 1) Read OVITO-exported data
data = read_lammps_data("corrected.data"; verbose=false)

# Save original labels and masses
orig_labels = copy(data.atom_labels)
orig_masses = deepcopy(data.masses)

# corrected.data type map:
# 1 = Al     4 = Oas    7 = Si_a
# 2 = Hb     5 = Ob     8 = Si_b
# 3 = O(Oss) 6 = Si

# 2) Collapse labels for distance-based bond detection
for i in eachindex(data.atom_labels)
    t = data.atom_labels[i]
    if t in (6, 7, 8);     data.atom_labels[i] = 6;  end   # all Si → 6
    if t in (3, 4, 5);     data.atom_labels[i] = 3;  end   # all O  → 3
end

# 3) Build coarse topology
add_zeolite_topology!(data;
    si_type=6, o_type=3, al_type=1, h_type=2,
    si_o_cutoff=1.85, al_o_cutoff=2.00, o_h_cutoff=1.05,
    verbose=true)

# 4) Restore original 8-type labels
data.atom_labels = orig_labels
data.masses = orig_masses

# 5) Refine topology types to match .ff numbering
refine_topology_types!(data; verbose=true)

# 6) Zero image flags (avoid LAMMPS "inconsistent image flags" warning)
data.image_flags = zeros(Int, size(data.coords, 1), 3)

# 7) Write output
write_lammps_data("corrected_with_topology.data", data;
    comment="H-ZSM-5 with refined Hill-Sauer topology (refine_topology_types v2)")

println("\n✓ corrected_with_topology.data written")
println("  Use with:  include hillsauer_alumsil_empty.ff")
println("  Bond types 1-6, Angle types 1-10, Dihedral types 1-10, Improper type 1")
