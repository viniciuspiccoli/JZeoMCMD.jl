# ════════════════════════════════════════════════════════════════
# alumino_support.jl — Aluminosilicate H-ZSM-5 support (8-type)
#
# Atom type convention (from corrected.data / OVITO):
#   1=Al  2=Hb  3=O(Oss)  4=Oas  5=Ob  6=Si  7=Si_a  8=Si_b
#   --- adsorbate (when loaded) ---
#   9=CH3_eth  10=CH2_eth  11=O_eth  12=H_eth
#
# Provides:
#   - AluminoConfig()       → ZeoliteConfig with correct type map
#   - build_alumino_data!() → full topology pipeline for corrected.data
#   - alumino_workflow_params() → WorkflowParams for aluminosilicate
# ════════════════════════════════════════════════════════════════

# ── Constants ──
const ALUMSIL_FW_TYPES = [1, 2, 3, 4, 5, 6, 7, 8]
const ALUMSIL_ADS_TYPES = [9, 10, 11, 12]
const ALUMSIL_NFW_ATOMS_2x2x2 = 2344   # H-ZSM-5 2×2×2 supercell
const ALUMSIL_TYPE_ELEMENTS = Dict(
    1=>"Al", 2=>"Hb", 3=>"O", 4=>"Oas", 5=>"Ob",
    6=>"Si", 7=>"Si_a", 8=>"Si_b",
    9=>"C", 10=>"C", 11=>"O", 12=>"H")


const ALUMSIL_CHARGES = Dict(
    1 =>  0.5366,   # Al
    2 =>  0.0839,   # Hb
    3 => -0.2618,   # Oss
    4 => -0.2959,   # Oas
    5 => -0.2515,   # Ob
    6 =>  0.5236,   # Si
    7 =>  0.5192,   # Si_a
    8 =>  0.5319)   # Si_b

const ALUMSIL_MASSES = Dict(
    1 => 26.981538, 2 => 1.00794,  3 => 15.9994,
    4 => 15.9994,   5 => 15.9994,  6 => 28.0855,
    7 => 28.0855,   8 => 28.0855)

# Ethanol (TraPPE-UA)
const ETH_CHARGES = Dict(
    "CH3" => 0.000, "CH2" => 0.265,
    "O_eth" => -0.700, "H_eth" => 0.435)
const ETH_MASSES = Dict(
    "CH3" => 15.035, "CH2" => 14.027,
    "O_eth" => 15.999, "H_eth" => 1.008)
const ETH_TYPES = Dict(
    "CH3" => 9, "CH2" => 10, "O_eth" => 11, "H_eth" => 12)

# ════════════════════════════════════════════════════════════════
# Configuration constructors
# ════════════════════════════════════════════════════════════════

"""
    AluminoConfig(; kwargs...) -> ZeoliteConfig

Create a ZeoliteConfig for aluminosilicate H-ZSM-5 with the 8-type
atom convention from corrected.data.

Type map:
  1=Al 2=Hb 3=Oss 4=Oas 5=Ob 6=Si 7=Si_a 8=Si_b
  9=CH3 10=CH2 11=O_eth 12=H_eth
"""
function AluminoConfig(; kwargs...)
    ZeoliteConfig(;
        si_type     = 6,
        o_type      = 3,
        al_type     = 1,
        h_acid_type = 2,

        # No remapping needed — corrected.data is already in final convention
        ovito_type_remap = Dict{Int,Int}(),

        eth_atoms_per_mol = 4,
        eth_atom_names = ["CH3", "CH2", "O_eth", "H_eth"],
        eth_types   = ETH_TYPES,
        eth_charges = ETH_CHARGES,
        eth_masses  = ETH_MASSES,

        # Ethanol bonded topology (continues after 6 fw bond types, 10 angle, 10 dihedral)
        eth_bond_defs     = [(7,1,2), (8,2,3), (9,3,4)],
        eth_angle_defs    = [(11,1,2,3), (12,2,3,4)],
        eth_dihedral_defs = [(11,1,2,3,4)],

        kwargs...
    )
end

"""
    alumino_workflow_params(; base_dir=pwd(), kwargs...) -> WorkflowParams

Create WorkflowParams preconfigured for aluminosilicate H-ZSM-5.
"""
function alumino_workflow_params(; base_dir::String=pwd(), kwargs...)
    WorkflowParams(;
        base_dir        = base_dir,
        nfw_atoms       = ALUMSIL_NFW_ATOMS_2x2x2,
        atoms_per_mol   = 4,
        table_file      = "hillsauer_alumsil_nb.table",
        is_alumino      = true,
        ff_include      = "hillsauer_alumsil_loaded.ff",
        kwargs...
    )
end

# ════════════════════════════════════════════════════════════════
# Topology pipeline for aluminosilicate
# ════════════════════════════════════════════════════════════════

"""
    build_alumino_topology!(data; verbose=true)

Full topology pipeline for aluminosilicate data from OVITO:
  1. Collapse types (Si*→6, O*→3) for distance-based detection
  2. Build coarse topology via add_zeolite_topology!
  3. Restore original 8-type labels
  4. Refine bond/angle/dihedral types via refine_topology_types!
  5. Zero image flags

Input `data` must have the 8-type atom labels from corrected.data.
After this call, the data file is ready for use with hillsauer_alumsil.ff.
"""
function build_alumino_topology!(data; verbose::Bool=true)
    # Save originals
    orig_labels = copy(data.atom_labels)
    orig_masses = deepcopy(data.masses)

    # Step 1: Collapse for topology detection
    for i in eachindex(data.atom_labels)
        t = data.atom_labels[i]
        if t in (6, 7, 8);  data.atom_labels[i] = 6;  end  # all Si
        if t in (3, 4, 5);  data.atom_labels[i] = 3;  end  # all O
    end

    # Step 2: Build coarse topology
    add_zeolite_topology!(data;
        si_type=6, o_type=3, al_type=1, h_type=2,
        si_o_cutoff=1.85, al_o_cutoff=2.00, o_h_cutoff=1.05,
        verbose=verbose)

    # Step 3: Restore labels
    data.atom_labels = orig_labels
    data.masses = orig_masses

    # Step 4: Refine types to match .ff numbering
    refine_topology_types!(data; verbose=verbose)

    # Step 5: Zero image flags
    data.image_flags = zeros(Int, size(data.coords, 1), 3)

    return data
end

"""
    assign_alumino_charges!(data)

Assign Hill-Sauer bond increment charges to framework atoms.
Only touches types 1-8; leaves adsorbate charges untouched.
"""
function assign_alumino_charges!(data)
    for j in eachindex(data.atom_labels)
        t = data.atom_labels[j]
        if haskey(ALUMSIL_CHARGES, t)
            data.atom_charges[j] = ALUMSIL_CHARGES[t]
        end
    end
    total = sum(data.atom_charges[j] for j in eachindex(data.atom_labels)
                if data.atom_labels[j] in ALUMSIL_FW_TYPES)
    @printf("  Framework charge: %.6f e\n", total)
    return data
end

# ════════════════════════════════════════════════════════════════
# Adsorbate reload helper
# ════════════════════════════════════════════════════════════════

"""
    AluminoReloadConfig(; kwargs...) -> ReloadConfig

Configuration for reloading adsorbate into an aluminosilicate framework.
Framework types = 1-8, adsorbate types = 9-12.
"""
function AluminoReloadConfig(; kwargs...)
    ReloadConfig(;
        adsorbate_types = [9, 10, 11, 12],
        framework_types = [1, 2, 3, 4, 5, 6, 7, 8],
        si_type = 6,
        o_type  = 3,
        eth_atoms_per_mol = 4,
        eth_atom_names = ["CH3", "CH2", "O_eth", "H_eth"],
        eth_types   = Dict("CH3"=>9, "CH2"=>10, "O_eth"=>11, "H_eth"=>12),
        eth_charges = Dict("CH3"=>0.0, "CH2"=>0.265, "O_eth"=>-0.700, "H_eth"=>0.435),
        eth_masses  = Dict("CH3"=>15.035, "CH2"=>14.027, "O_eth"=>15.999, "H_eth"=>1.008),
        eth_bond_defs     = [(7,1,2), (8,2,3), (9,3,4)],
        eth_angle_defs    = [(11,1,2,3), (12,2,3,4)],
        eth_dihedral_defs = [(11,1,2,3,4)],
        kwargs...
    )
end

# ════════════════════════════════════════════════════════════════
# CIF writer helper
# ════════════════════════════════════════════════════════════════

"""
    write_alumino_cif(filename, data; comment="")

Write CIF for aluminosilicate, mapping all 8 atom types correctly.
Framework only (strips adsorbate atoms with types > 8).
"""
function write_alumino_cif(filename::String, data; comment::String="")
    write_cif(filename, data;
        framework_types = ALUMSIL_FW_TYPES,
        type_elements   = ALUMSIL_TYPE_ELEMENTS,
        comment         = comment)
end
