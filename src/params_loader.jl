# ════════════════════════════════════════════════════════════════
# params_loader.jl — Read params.toml → ZeoliteConfig
# Updated to match the corrected params.toml (v3) structure
# ════════════════════════════════════════════════════════════════

import TOML

"""
    load_config(toml_path; is_alumino=false, kwargs...) -> ZeoliteConfig

Read params.toml and return a fully populated ZeoliteConfig.
Set `is_alumino=true` for aluminosilicate type numbering.
"""
function load_config(toml_path::String; is_alumino::Bool=false, kwargs...)
    !isfile(toml_path) && error("params.toml not found: $toml_path")
    p = TOML.parsefile(toml_path)

    eth = p["ethanol"]
    ps  = p["pair_style"]
    kB  = get(get(p, "units", Dict()), "kB_kcal", 0.0019872041)

    names = eth["atom_names"]

    if is_alumino
        al = p["alumino"]
        types_v = eth["types_alumino"]
        cfg = ZeoliteConfig(;
            si_type      = al["Si_type"],
            o_type       = al["Oss_type"],
            al_type      = al["Al_type"],
            h_acid_type  = al["Hb_type"],
            ovito_type_remap = Dict(1 => 2, 2 => 1),
            eth_atoms_per_mol = eth["atoms_per_mol"],
            eth_atom_names    = names,
            eth_types   = Dict(names[i] => types_v[i] for i in 1:length(names)),
            eth_charges = Dict(names[i] => eth["charges"][i] for i in 1:length(names)),
            eth_masses  = Dict(names[i] => eth["masses"][i] for i in 1:length(names)),
            pair_cutoff = ps["cutoff"],
            coul_cutoff = ps["coul_cutoff"],
            kwargs...
        )
    else
        si = p["silica"]
        types_v = eth["types_silica"]
        cfg = ZeoliteConfig(;
            si_type      = si["Si_type"],
            o_type       = si["O_type"],
            al_type      = 0,
            h_acid_type  = 0,
            ovito_type_remap = Dict(1 => 2, 2 => 1),
            eth_atoms_per_mol = eth["atoms_per_mol"],
            eth_atom_names    = names,
            eth_types   = Dict(names[i] => types_v[i] for i in 1:length(names)),
            eth_charges = Dict(names[i] => eth["charges"][i] for i in 1:length(names)),
            eth_masses  = Dict(names[i] => eth["masses"][i] for i in 1:length(names)),
            pair_cutoff = ps["cutoff"],
            coul_cutoff = ps["coul_cutoff"],
            kwargs...
        )
    end

    return cfg
end

"""
    load_ff_params(toml_path) -> Dict

Load the raw TOML as a nested Dict.
"""
function load_ff_params(toml_path::String)
    return TOML.parsefile(toml_path)
end
