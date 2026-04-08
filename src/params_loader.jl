# ════════════════════════════════════════════════════════════════
# params_loader.jl — Read params.toml and populate ZeoliteConfig
#
# This replaces all hardcoded FF values. Every script reads from
# one params.toml so changing a parameter propagates everywhere.
#
# Usage:
#   include("params_loader.jl")
#   cfg = load_config("ff/params.toml")
#   # or:
#   cfg = load_config("ff/params.toml"; ovito_data="test.data",
#                     raspa_restart="restart.json")
# ════════════════════════════════════════════════════════════════

#import TOML

"""
    load_config(toml_path; overrides...) -> ZeoliteConfig

Read `params.toml` and return a fully populated `ZeoliteConfig`.
Any keyword argument overrides the TOML value.

Example:
    cfg = load_config("ff/params.toml";
              ovito_data = "test.data",
              raspa_restart = "restart.json",
              output_data = "loaded.lmp")
"""
function load_config(toml_path::String; kwargs...)
    !isfile(toml_path) && error("params.toml not found: $toml_path")
    p = TOML.parsefile(toml_path)

    fw  = p["framework"]
    eth = p["ethanol"]
    ps  = p["pair_style"]

    # ── Build ethanol dicts from TOML arrays ──
    names   = eth["atom_names"]
    types_v = eth["types"]
    q_v     = eth["charges"]
    m_v     = eth["masses"]

    eth_types   = Dict(names[i] => types_v[i] for i in 1:length(names))
    eth_charges = Dict(names[i] => q_v[i]     for i in 1:length(names))
    eth_masses  = Dict(names[i] => m_v[i]     for i in 1:length(names))

    # ── Framework charges ──
    fw_charges = fw["charges"]

    # ── Determine table file ──
    al_type = get(fw, "Al_type", 0)
    table_file = al_type > 0 ? ps["table_file_alumsil"] : ps["table_file_silica"]

    # ── Build config ──
    cfg = ZeoliteConfig(
        # Files (overridable)
        ovito_data    = get(kwargs, :ovito_data, "test.data"),
        raspa_restart = get(kwargs, :raspa_restart, ""),
        table_file    = get(kwargs, :table_file, table_file),
        output_data   = get(kwargs, :output_data, "loaded_zeolite.lmp"),
        output_input  = get(kwargs, :output_input, "run_loaded.in"),
        fw_check_xyz  = get(kwargs, :fw_check_xyz, "framework_check.xyz"),
        eth_check_xyz = get(kwargs, :eth_check_xyz, "ethanol_check.xyz"),

        # Ovito remap
        ovito_type_remap = Dict(1 => 2, 2 => 1),  # O→2, Si→1

        # Framework types
        si_type     = fw["Si_type"],
        o_type      = fw["O_type"],
        al_type     = al_type,
        h_acid_type = get(fw, "H_acid_type", 0),

        # Topology cutoffs
        si_o_cutoff = get(kwargs, :si_o_cutoff, 1.85),
        al_o_cutoff = get(kwargs, :al_o_cutoff, 2.00),
        o_h_cutoff  = get(kwargs, :o_h_cutoff, 1.05),

        # Ethanol
        eth_atoms_per_mol = eth["atoms_per_mol"],
        eth_atom_names    = names,
        eth_types         = eth_types,
        eth_charges       = eth_charges,
        eth_masses        = eth_masses,

        # Pair
        pair_cutoff = ps["cutoff"],
        coul_cutoff = ps["coul_cutoff"],

        box_tolerance = get(kwargs, :box_tolerance, 0.5),
    )

    return cfg
end

"""
    load_ff_params(toml_path) -> Dict

Load the raw TOML as a nested Dict for use by write_complete_data
and other functions that need specific bonded coefficients.
"""
function load_ff_params(toml_path::String)
    return TOML.parsefile(toml_path)
end

# ════════════════════════════════════════════════════════════════
# write_complete_data_from_toml — TOML-driven data file writer
# ════════════════════════════════════════════════════════════════

"""
    write_complete_data(fname, data, cfg; params_toml="ff/params.toml")

Write a LAMMPS data file with all class2 coefficients.
Coefficients are read from params.toml instead of being hardcoded.
Falls back to hardcoded values if TOML is not found.
"""
function write_complete_data_toml(fname::String, d, cfg::ZeoliteConfig;
                                   params_toml::String = "ff/params.toml")
    # Try to load TOML; fall back to original hardcoded writer if not found
    if !isfile(params_toml)
        @warn "params.toml not found at $params_toml — using hardcoded coefficients"
        write_complete_data(fname, d, cfg)
        return
    end

    p = TOML.parsefile(params_toml)
    fw_bond = p["framework"]["bonded"]
    eth_bond = p["ethanol"]["bonds"]
    eth_angle = p["ethanol"]["angles"]
    eth_torsion = p["ethanol"]["torsion"]

    natoms = size(d.coords, 1)
    nbt = d.nbond_types; nat = d.nangle_types; ndt = d.ndihedral_types
    natypes = maximum(keys(d.masses))
    is_tri = any(d.tilt_factors .!= 0.0)
    has_eth = natypes > 2

    # Convert ethanol bond K from RASPA convention (K/Å²) to class2 (kcal/mol/Å²)
    kB = 0.0019872041
    K_eth_bond = first(values(eth_bond))[2] * kB  # all same K in TraPPE

    # Convert ethanol angle K from K/rad² to kcal/mol/rad²
    K_eth_ang_1 = eth_angle["CH3_CH2_O"][2] * kB
    K_eth_ang_2 = eth_angle["CH2_O_H"][2] * kB

    # Convert torsion from K to kcal/mol
    torsion_K = eth_torsion["CH3_CH2_O_H"]
    torsion_kcal = [c * kB for c in torsion_K]

    # Si-O bond params
    sio = fw_bond["Si_O"]
    # Angle params
    siosi = fw_bond["Si_O_Si"]
    osio  = fw_bond["O_Si_O"]
    # Dihedral
    osiosi = fw_bond["O_Si_O_Si"]

    open(fname, "w") do io
        println(io, "LAMMPS data — loaded zeolite (params.toml-driven)\n")
        println(io, "$(natoms) atoms")
        println(io, "$(size(d.bonds,1)) bonds")
        println(io, "$(size(d.angles,1)) angles")
        println(io, "$(size(d.dihedrals,1)) dihedrals")
        println(io, "0 impropers\n")
        println(io, "$natypes atom types")
        println(io, "$nbt bond types")
        println(io, "$nat angle types")
        println(io, "$ndt dihedral types")
        println(io, "0 improper types\n")
        @printf(io, "%.10f %.10f xlo xhi\n", d.box_dimensions[1,1], d.box_dimensions[1,2])
        @printf(io, "%.10f %.10f ylo yhi\n", d.box_dimensions[2,1], d.box_dimensions[2,2])
        @printf(io, "%.10f %.10f zlo zhi\n", d.box_dimensions[3,1], d.box_dimensions[3,2])
        is_tri && @printf(io, "%.10f %.10f %.10f xy xz yz\n",
                          d.tilt_factors[1], d.tilt_factors[2], d.tilt_factors[3])

        mn = Dict(1=>"Si",2=>"O",3=>"CH3_eth",4=>"CH2_eth",5=>"O_eth",6=>"H_eth")
        println(io, "\nMasses\n")
        for t in sort(collect(keys(d.masses)))
            @printf(io, "  %d  %.6f  # %s\n", t, d.masses[t], get(mn,t,""))
        end

        # ── Bond Coeffs (from TOML) ──
        println(io, "\nBond Coeffs # class2\n")
        @printf(io, "  1  %.4f  %.4f  %.4f  %.4f  # Si-O\n", sio...)
        if has_eth
            r1 = eth_bond["CH3_CH2"][1]
            r2 = eth_bond["CH2_O_eth"][1]
            r3 = eth_bond["O_eth_H"][1]
            @printf(io, "  2  %.3f  %.2f  0.0  0.0  # CH3-CH2\n", r1, K_eth_bond)
            @printf(io, "  3  %.3f  %.2f  0.0  0.0  # CH2-O_eth\n", r2, K_eth_bond)
            @printf(io, "  4  %.3f  %.2f  0.0  0.0  # O_eth-H_eth\n", r3, K_eth_bond)
        end

        # ── Angle Coeffs ──
        println(io, "\nAngle Coeffs # class2\n")
        @printf(io, "  1  %.1f  %.4f  %.4f  %.4f  # Si-O-Si\n", siosi...)
        @printf(io, "  2  %.1f  %.4f  %.4f  %.4f  # O-Si-O\n", osio...)
        if has_eth
            @printf(io, "  3  %.2f  %.2f  0.0  0.0  # CH3-CH2-O_eth\n",
                    eth_angle["CH3_CH2_O"][1], K_eth_ang_1)
            @printf(io, "  4  %.2f  %.2f  0.0  0.0  # CH2-O_eth-H_eth\n",
                    eth_angle["CH2_O_H"][1], K_eth_ang_2)
        end

        # ── BondBond Coeffs ──
        println(io, "\nBondBond Coeffs\n")
        println(io, "  1  151.8742  $(sio[1])  $(sio[1])")
        println(io, "  2  0.0  $(sio[1])  $(sio[1])")
        if has_eth
            r1 = eth_bond["CH3_CH2"][1]; r2 = eth_bond["CH2_O_eth"][1]; r3 = eth_bond["O_eth_H"][1]
            println(io, "  3  0.0  $r1  $r2")
            println(io, "  4  0.0  $r2  $r3")
        end

        # ── BondAngle Coeffs ──
        println(io, "\nBondAngle Coeffs\n")
        println(io, "  1  9.2390  9.2390  $(sio[1])  $(sio[1])")
        println(io, "  2  78.1239  78.1239  $(sio[1])  $(sio[1])")
        if has_eth
            r1 = eth_bond["CH3_CH2"][1]; r2 = eth_bond["CH2_O_eth"][1]; r3 = eth_bond["O_eth_H"][1]
            println(io, "  3  0.0  0.0  $r1  $r2")
            println(io, "  4  0.0  0.0  $r2  $r3")
        end

        # ── Dihedral Coeffs ──
        println(io, "\nDihedral Coeffs # class2\n")
        @printf(io, "  1  %.4f  %.1f  %.4f  %.1f  %.4f  %.1f  # O-Si-O-Si\n", osiosi...)
        if has_eth
            # Convert TRAPPE [c0,c1,c2,c3] to class2 [V1,0,V2,0,V3,0]
            # c1→V1, c2→V2, c3→V3 (all ×kB, phase=0)
            @printf(io, "  2  %.4f  0.0  %.4f  0.0  %.4f  0.0  # CH3-CH2-O-H\n",
                    torsion_kcal[2], torsion_kcal[3], torsion_kcal[4])
        end

        # ── MiddleBondTorsion ──
        println(io, "\nMiddleBondTorsion Coeffs\n")
        println(io, "  1  0.0  0.0  0.0  $(sio[1])")
        has_eth && println(io, "  2  0.0  0.0  0.0  $(eth_bond["CH2_O_eth"][1])")

        # ── EndBondTorsion ──
        println(io, "\nEndBondTorsion Coeffs\n")
        println(io, "  1  0.0  0.0  0.0  0.0  0.0  0.0  $(sio[1])  $(sio[1])")
        if has_eth
            r1 = eth_bond["CH3_CH2"][1]; r3 = eth_bond["O_eth_H"][1]
            println(io, "  2  0.0  0.0  0.0  0.0  0.0  0.0  $r1  $r3")
        end

        # ── AngleTorsion ──
        println(io, "\nAngleTorsion Coeffs\n")
        @printf(io, "  1  0.0  0.0  0.0  0.0  0.0  0.0  %.1f  %.1f\n", osio[1], siosi[1])
        if has_eth
            @printf(io, "  2  0.0  0.0  0.0  0.0  0.0  0.0  %.2f  %.2f\n",
                    eth_angle["CH3_CH2_O"][1], eth_angle["CH2_O_H"][1])
        end

        # ── AngleAngleTorsion ──
        println(io, "\nAngleAngleTorsion Coeffs\n")
        @printf(io, "  1  -4.5150  %.1f  %.1f\n", osio[1], siosi[1])
        if has_eth
            @printf(io, "  2  0.0  %.2f  %.2f\n",
                    eth_angle["CH3_CH2_O"][1], eth_angle["CH2_O_H"][1])
        end

        # ── BondBond13 ──
        println(io, "\nBondBond13 Coeffs\n")
        println(io, "  1  0.0  $(sio[1])  $(sio[1])")
        if has_eth
            r1 = eth_bond["CH3_CH2"][1]; r3 = eth_bond["O_eth_H"][1]
            println(io, "  2  0.0  $r1  $r3")
        end

        # ── Atoms ──
        println(io, "\nAtoms # full\n")
        for j in 1:natoms
            @printf(io, "  %d %d %d %.6f %.10f %.10f %.10f %d %d %d\n",
                    d.atom_ids[j], d.molecule_labels[j], d.atom_labels[j],
                    d.atom_charges[j], d.coords[j,1], d.coords[j,2], d.coords[j,3],
                    d.image_flags[j,1], d.image_flags[j,2], d.image_flags[j,3])
        end

        # ── Topology ──
        for (secname, labels, atoms) in [
            ("Bonds", d.bond_labels, d.bonds),
            ("Angles", d.angle_labels, d.angles),
            ("Dihedrals", d.dihedral_labels, d.dihedrals)]
            size(atoms,1) == 0 && continue
            println(io, "\n$secname\n")
            for k in 1:size(atoms,1)
                print(io, "  $k $(labels[k])")
                for c in 1:size(atoms,2); print(io, " $(atoms[k,c])"); end
                println(io)
            end
        end
    end
    println("  Wrote $fname (params from TOML)")
end
