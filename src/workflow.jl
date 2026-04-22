# ════════════════════════════════════════════════════════════════
# workflow.jl — Master GCMC/MD loop (included by JZeoMCMD.jl)
#
# Directory structure (matches run_test.jl):
#   base_dir/
#   ├── hillsauer_nb.table
#   ├── run_npt.in
#   ├── MFI_SI.cif
#   ├── convergence.csv
#   ├── cycle_01/
#   │   ├── raspa/
#   │   │   ├── simulation.json
#   │   │   └── output/           ← RASPA3 creates this
#   │   │       └── restart_*.json
#   │   ├── lammps/
#   │   │   ├── loaded.lmp
#   │   │   ├── run_npt.in
#   │   │   └── npt_final.data
#   │   ├── loaded.lmp
#   │   └── distorted.cif
#   └── cycle_02/
#       └── ...
# ════════════════════════════════════════════════════════════════

import TOML, Dates

"""
    WorkflowParams

Configuration for the GCMC/MD iterative workflow.
All paths are absolute or relative to the current working directory.
Package defaults come from JZeoMCMD resource paths.
"""
Base.@kwdef mutable struct WorkflowParams
    # ── Simulation directory (user chooses this) ──
    base_dir::String = pwd()

    # ── Initial files (in user's project) ──
    initial_cif::String  = ""
    initial_data::String = ""

    # ── Package resources ──
    params_toml::String     = params_path()
    raspa_input_dir::String = raspa_inputs_path()
    lammps_input::String    = ""    # path to run_npt.in in base_dir

    # ── Executables ──
    raspa_exe::String  = "raspa3"
    lammps_exe::String = "lmp_mpi"

    # ── Physics ──
    temperature::Float64 = 300.0
    pressure::Float64    = 1e5

    # ── RASPA3 GCMC settings ──
    raspa_n_init::Int    = 10000    # NumberOfInitializationCycles
    raspa_n_equil::Int   = 10000   # NumberOfEquilibrationCycles
    raspa_n_prod::Int    = 50000   # NumberOfCycles (production)
    raspa_print_every::Int = 1000

    # ── LAMMPS NPT settings ──
    lammps_minimize_steps::Int  = 5000    # minimize iterations
    lammps_nvt_ramp_steps::Int  = 240000  # total NVT ramp (40k+80k+120k)
    lammps_nvt_hold_steps::Int  = 200000  # NVT hold at T
    lammps_npt_steps::Int       = 500000  # NPT production
    lammps_timestep::Float64    = 0.25    # fs

    # ── Iteration control ──
    max_iterations::Int     = 10
    convergence_tol::Float64 = 0.005
    convergence_window::Int  = 5

    # ── Derived ──
    table_file::String = "hillsauer_nb.table"
    adsorbate_component::String = "ethanol"
    atoms_per_mol::Int = 4         # ethanol = CH3 + CH2 + O + H
    n_unit_cells::Int = 8          # 2×2×2
    nfw_atoms::Int = 2304          # MFI 2×2×2

    # ── Aluminosilicate support ──
    is_alumino::Bool = false       # true → use 8-type convention + .ff includes
    ff_include::String = ""        # e.g. "hillsauer_alumsil_loaded.ff" (copied to lammps dir)

end

"""
    load_workflow_params(toml_path; base_dir=pwd(), kwargs...)

Load workflow settings from params.toml.
"""
function load_workflow_params(toml_path::String; base_dir::String=pwd(), kwargs...)
    p = TOML.parsefile(toml_path)
    wf = get(p, "workflow", Dict())

    wp = WorkflowParams(
        base_dir            = base_dir,
        params_toml         = toml_path,
        raspa_input_dir     = get(wf, "raspa_input_dir", raspa_inputs_path()),
        lammps_input        = get(wf, "lammps_input", ""),
        raspa_exe           = get(wf, "raspa_exe", "raspa3"),
        lammps_exe          = get(wf, "lammps_exe", "lmp_mpi"),
        table_file          = get(wf, "table_file", "hillsauer_nb.table"),
        adsorbate_component = get(wf, "adsorbate_component", "ethanol"),
    )

    for (k, v) in kwargs
        hasproperty(wp, k) && setproperty!(wp, k, v)
    end
    return wp
end

# ════════════════════════════════════════════════════════════════
# Auto-generate table in base_dir if missing
# ════════════════════════════════════════════════════════════════

function ensure_tables!(wp::WorkflowParams)
    full = joinpath(wp.base_dir, wp.table_file)
    isfile(full) && return full

    # Try copying from package
    pkg_table = joinpath(resource_dir(), "ff", "hillsauer_silica.table")
    if isfile(pkg_table)
        cp(pkg_table, full; force=true)
        println("  Copied table from package: $full")
        return full
    end

    # Auto-generate from params.toml
    p = TOML.parsefile(wp.params_toml)
    nb = p["framework"]["nonbonded"]
    ps = p["pair_style"]
    A_Si = nb["A_Si"]; A_O = nb["A_O"]
    A_Al = get(nb, "A_Al", A_Si)
    rmin = ps["table_rmin"]; rmax = ps["table_rmax"]; N = ps["table_npoints"]
    al_type = get(p["framework"], "Al_type", 0)

    pairs = [("Si_Si",A_Si,A_Si), ("Si_O",A_Si,A_O), ("O_O",A_O,A_O)]
    al_type > 0 && append!(pairs, [("Si_Al",A_Si,A_Al), ("O_Al",A_O,A_Al), ("Al_Al",A_Al,A_Al)])

    open(full, "w") do io
        println(io, "# Hill-Sauer A/r^9 table (auto-generated)\n")
        for (label, Ai, Aj) in pairs
            A = sqrt(Ai*Aj)
            println(io, label); println(io, "N $N R $rmin $rmax\n")
            dr = (rmax-rmin)/(N-1)
            for k in 1:N
                r = rmin+(k-1)*dr
                @printf(io, "%d  %.6f  %.6f  %.6f\n", k, r, A/r^9, 9A/r^10)
            end
            println(io)
        end
    end
    println("  Auto-generated table: $full")
    return full
end

# ════════════════════════════════════════════════════════════════
# Folder management
# ════════════════════════════════════════════════════════════════

function cycle_dir(wp::WorkflowParams, n::Int)
    d = joinpath(wp.base_dir, @sprintf("cycle_%02d", n))
    mkpath(d)
    return d
end

# ════════════════════════════════════════════════════════════════
# Step 1: RASPA3 GCMC
#   - Copy raspa_input_dir → cycle_NN/raspa/
#   - Cycle 1: uses simulation.json.template (unit cell CIF)
#   - Cycle 2+: uses simulation.json.template_next (supercell CIF, 1×1×1)
#   - Patches __FRAMEWORK__, __TEMPERATURE__, __PRESSURE__
#   - RASPA3 output goes to cycle_NN/raspa/output/
# ════════════════════════════════════════════════════════════════

function step_gcmc!(wp::WorkflowParams, cycle::Int, cif_path::String)
    rdir = joinpath(cycle_dir(wp, cycle), "raspa")
    mkpath(rdir)

    # Copy all RASPA input files from package
    pkg_raspa = wp.raspa_input_dir
    if !isdir(pkg_raspa)
        @warn "RASPA input dir not found: $pkg_raspa"
        return nothing
    end
    for f in readdir(pkg_raspa; join=true)
        isfile(f) && cp(f, joinpath(rdir, basename(f)); force=true)
    end

    # Copy CIF
    cif_name = basename(cif_path)
    cp(cif_path, joinpath(rdir, cif_name); force=true)
    fw_name = replace(cif_name, ".cif" => "")

    # Select template based on cycle
    sim_file = joinpath(rdir, "simulation.json")
    if cycle == 1
        sim_template = joinpath(rdir, "simulation.json.template")
    else
        sim_template = joinpath(rdir, "simulation.json.template_next")
    end

    if isfile(sim_template)
        txt = read(sim_template, String)
        txt = replace(txt,
            "__FRAMEWORK__"   => fw_name,
            "__TEMPERATURE__" => string(wp.temperature),
            "__PRESSURE__"    => string(wp.pressure),
            "__N_INIT__"      => string(wp.raspa_n_init),
            "__N_EQUIL__"     => string(wp.raspa_n_equil),
            "__N_PROD__"      => string(wp.raspa_n_prod),
            "__PRINT_EVERY__" => string(wp.raspa_print_every),
        )
        write(sim_file, txt)
        # Clean up unused templates
        for t in ["simulation.json.template", "simulation.json.template_next"]
            tf = joinpath(rdir, t)
            isfile(tf) && rm(tf)
        end
    elseif isfile(sim_file)
        # Patch existing JSON directly
        content = JSON.parsefile(sim_file)
        if haskey(content, "Systems") && !isempty(content["Systems"])
            content["Systems"][1]["ExternalPressure"] = wp.pressure
            content["Systems"][1]["ExternalTemperature"] = wp.temperature
            content["Systems"][1]["Name"] = fw_name
        end
        haskey(content, "NumberOfCycles") &&
            (content["NumberOfCycles"] = wp.raspa_n_prod)
        haskey(content, "NumberOfInitializationCycles") &&
            (content["NumberOfInitializationCycles"] = wp.raspa_n_init)
        haskey(content, "NumberOfEquilibrationCycles") &&
            (content["NumberOfEquilibrationCycles"] = wp.raspa_n_equil)
        haskey(content, "PrintEvery") &&
            (content["PrintEvery"] = wp.raspa_print_every)
        open(sim_file, "w") do io
            JSON.print(io, content, 2)
        end
    else
        error("No simulation.json or template found in $rdir")
    end

    # Run RASPA3
    println("  [GCMC] Running RASPA3 (P=$(wp.pressure) Pa)...")
    t0 = time()
    try
        run(Cmd(Cmd(split(wp.raspa_exe)); dir=rdir))
    catch e
        @warn "RASPA3 failed" exception=e
        return nothing
    end
    @printf("  [GCMC] Done (%.1f s)\n", time()-t0)

    # Find restart JSON in output/ subdirectory (RASPA3 convention)
    output_dir = joinpath(rdir, "output")
    if !isdir(output_dir)
        @warn "No output/ directory in $rdir"
        return nothing
    end

    restarts = sort(filter(f -> startswith(f, "restart_") && endswith(f, ".json"),
                           readdir(output_dir)))
    if isempty(restarts)
        @warn "No restart JSON found in $output_dir"
        return nothing
    end
    restart_path = joinpath(output_dir, restarts[end])

    # Parse N_ads from restart JSON (most reliable)
    raspa_out = Dict{String,Any}("n_ads" => NaN)
    try
        restart_data = JSON.parsefile(restart_path)
        for (k, v) in restart_data
            #if v isa AbstractVector && !isempty(v) && v[1] isa AbstractVector
            #    raspa_out["n_ads"] = Float64(length(v) ÷ wp.atoms_per_mol)
            #    break
            #end

            if v isa AbstractVector
                raspa_out["n_ads"] = isempty(v) ? 0.0 : Float64(length(v) ÷ wp.atoms_per_mol)
                break
            end

        end
    catch e
        @warn "Failed to parse restart JSON" exception=e
    end
    @printf("  [GCMC] N_ads = %.0f molecules\n", raspa_out["n_ads"])

    return (restart=restart_path, output=raspa_out, dir=rdir)
end

# ════════════════════════════════════════════════════════════════
# Step 2: Build / reload LAMMPS data file
# ════════════════════════════════════════════════════════════════

function step_build_data!(wp::WorkflowParams, cycle::Int,
                           restart_json::String,
                           prev_data::String)
    cdir = cycle_dir(wp, cycle)
    data_file = joinpath(cdir, "loaded.lmp")


    if wp.is_alumino
        # ── Aluminosilicate path ──
        if cycle == 1
            println("  [BUILD] Cycle 1: Alumino framework + RASPA3 ethanol...")
            data = read_lammps_data(wp.initial_data; verbose=false)
            build_alumino_topology!(data; verbose=true)

            # Read ethanol
            acfg = AluminoConfig(raspa_restart=restart_json)
            mols = read_raspa3_ethanol(acfg, data.box_dimensions)
            merge_framework_ethanol!(data, mols, acfg)
            write_lammps_data(data_file, data;
                comment="H-ZSM-5 + ethanol (build_alumino_topology)")
        else
            println("  [BUILD] Cycle $cycle: reload adsorbate (alumino)...")
            rcfg = ReloadConfig(
                adsorbate_types = [9,10,11,12],
                framework_types = [1,2,3,4,5,6,7,8],
                si_type = 6, o_type = 3,
                eth_types   = Dict("CH3"=>9,"CH2"=>10,"O_eth"=>11,"H_eth"=>12),
                eth_charges = Dict("CH3"=>0.0,"CH2"=>0.265,"O_eth"=>-0.700,"H_eth"=>0.435),
                eth_masses  = Dict("CH3"=>15.035,"CH2"=>14.027,"O_eth"=>15.999,"H_eth"=>1.008),
                eth_bond_defs     = [(7,1,2),(8,2,3),(9,3,4)],
                eth_angle_defs    = [(11,1,2,3),(12,2,3,4)],
                eth_dihedral_defs = [(11,1,2,3,4)],
            )
            reload_adsorbate(prev_data, restart_json, data_file; cfg=rcfg)
        end
    else
        if cycle == 1
            # First cycle: build from Ovito data + RASPA3 ethanol
            println("  [BUILD] Cycle 1: Ovito framework + RASPA3 ethanol...")
            # src = structures_path()
            cfg = ZeoliteConfig(
                ovito_data    = wp.initial_data,    #joinpath(src, "MFI_SI.data"),
                raspa_restart = restart_json,
                output_data   = data_file,
            )
            fw = read_and_remap_framework(cfg)
            add_framework_topology!(fw, cfg)
            mols = read_raspa3_ethanol(cfg, fw.box_dimensions)
            merge_framework_ethanol!(fw, mols, cfg)
            #write_complete_data(data_file, fw, cfg)
            write_lammps_data(data_file, fw; comment="H-ZSM-5 + ethanol (alumino cycle 1)")
        else
            # Subsequent cycles: reload adsorbate into previous NPT output
            println("  [BUILD] Cycle $cycle: reload adsorbate into NPT framework...")
            reload_adsorbate(prev_data, restart_json, data_file)
        end
    end

    # Verify
    check = read_lammps_data(data_file; verbose=false)
    ntot = size(check.coords, 1)
    neth = ntot - wp.nfw_atoms
    @printf("  [BUILD] Atoms: %d (%d fw + %d ethanol = %d molecules)\n",
            ntot, wp.nfw_atoms, neth, neth ÷ wp.atoms_per_mol)

    return data_file
end

# ════════════════════════════════════════════════════════════════
# Step 3: LAMMPS NPT-MD
#   - Copies run_npt.in + data + table to cycle_NN/lammps/
#   - Looks for npt_final.data in output
# ════════════════════════════════════════════════════════════════

function step_npt!(wp::WorkflowParams, cycle::Int, data_file::String)
    ldir = joinpath(cycle_dir(wp, cycle), "lammps")
    mkpath(ldir)

    # Copy data file as loaded.lmp
    data_name = "loaded.lmp"
    cp(data_file, joinpath(ldir, data_name); force=true)

    # Copy table
    table_src = joinpath(wp.base_dir, wp.table_file)
    cp(table_src, joinpath(ldir, wp.table_file); force=true)

    # Copy .ff include file (aluminosilicate)
  #  if wp.is_alumino && !isempty(wp.ff_include)
  #      ff_src = joinpath(wp.base_dir, wp.ff_include)
  #      isfile(ff_src) && cp(ff_src, joinpath(ldir, basename(wp.ff_include)); force=true)
  #  end

    # Copy .ff include file (if set)
    if !isempty(wp.ff_include)
	ff_src = joinpath(wp.base_dir, wp.ff_include)
        isfile(ff_src) && cp(ff_src, joinpath(ldir, basename(wp.ff_include)); force=true)
    end


    # Copy LAMMPS input script
    input_src = wp.lammps_input
    if !isfile(input_src)
        input_src = joinpath(wp.base_dir, "run_npt.in")
    end
    if !isfile(input_src)
        @warn "LAMMPS input not found: $input_src"
        return nothing
    end
    cp(input_src, joinpath(ldir, "run_npt.in"); force=true)

    # Run LAMMPS
    println("  [NPT] Running LAMMPS...")
    t0 = time()
    try
        #run(Cmd(Cmd(split(wp.lammps_exe)), `-in`, `run_npt.in`; dir=ldir))
        #lammps_cmd = Cmd(`$(split(wp.lammps)) -in run_npt.in`; dit=1dir)
        lammps_cmd = Cmd(`$(split(wp.lammps_exe)) -in run_npt.in`; dir=ldir)
        run(lammps_cmd)
    catch e
        @warn "LAMMPS failed" exception=e
        return nothing
    end
    @printf("  [NPT] Done (%.1f s)\n", time()-t0)

    # Find output data file
    output_data = ""
    for name in ["npt_final.data", "loaded_npt_final.data",
                  "npt_final.lmp", "loaded_npt_final.lmp"]
        f = joinpath(ldir, name)
        isfile(f) && (output_data = f; break)
    end
    if isempty(output_data)
        for f in readdir(ldir)
            fp = joinpath(ldir, f)
            if isfile(fp) && f != data_name && f != "run_npt.in" &&
               (endswith(f, ".data") || endswith(f, ".lmp"))
                output_data = fp
                break
            end
        end
    end
    if isempty(output_data)
        @warn "No LAMMPS output found. Files: $(readdir(ldir))"
        return nothing
    end
    println("  [NPT] Output: $(basename(output_data))")

    # Parse log
    log_file = joinpath(ldir, "log.lammps")
    log_data = isfile(log_file) ? parse_lammps_log(log_file) : Dict{String,Vector{Float64}}()
    cell = extract_cell_params(log_data)
    @printf("  [NPT] a=%.4f b=%.4f c=%.4f V=%.1f\n",
            cell.a, cell.b, cell.c, cell.volume)

    return (data=output_data, cell=cell, log=log_data, dir=ldir)
end

# ════════════════════════════════════════════════════════════════
# Step 4: Extract distorted CIF
# ════════════════════════════════════════════════════════════════

function step_write_cif!(wp::WorkflowParams, cycle::Int, data_file::String)
    cdir = cycle_dir(wp, cycle)
    cif_out = joinpath(cdir, "distorted.cif")

    data = read_lammps_data(data_file; verbose=false)

    if wp.is_alumino 
        fw_types = ALUMSIL_FW_TYPES
        type_elem = ALUMSIL_TYPE_ELEMENTS
    else
        #fw_types = [1, 2]
        #p = TOML.parsefile(wp.params_toml)
        #fw = get(p, "silica", get(p, "framework", Dict()))
        #al = get(fw, "Al_type", 0)
        #al > 0 && push!(fw_types, al)
        fw_types = [1, 2]
        type_elem = Dict(1=>"Si", 2=>"O")


    end
    #p = TOML.parsefile(wp.params_toml)
    #al = get(p["framework"], "Al_type", 0)

    #write_cif(cif_out, data;
    #          framework_types=fw_types,
    #          type_elements=Dict(1=>"Si", 2=>"O", 3=>"Al"),
    #          comment=@sprintf("cycle %d P=%.1e Pa", cycle, wp.pressure))

     write_cif(cif_out, data;
              framework_types=fw_types,
              type_elements=type_elem,
              comment=@sprintf("cycle %d P=%.1e Pa", cycle, wp.pressure))

    return cif_out
end

# ════════════════════════════════════════════════════════════════
# Step 5: Analysis — strain, loading, channel occupancy
# ════════════════════════════════════════════════════════════════

function step_analyze!(wp::WorkflowParams, cycle::Int,
                        cell::NamedTuple, raspa_out::Dict,
                        ref_cell::NamedTuple, data_file::String)
    strain = compute_strain(cell, ref_cell)

    #data = read_lammps_data(data_file; verbose=false)
    #occ = analyze_channel_occupancy(
    #    data.coords, data.atom_labels,
    #    Set([3, 4, 5, 6]),   # adsorbate types (silica)
    #    data.box_dimensions)

    data = read_lammps_data(data_file; verbose=false)
    ads_types = wp.is_alumino ? Set([9, 10, 11, 12]) : Set([3, 4, 5, 6])
    occ = analyze_channel_occupancy(
        data.coords, data.atom_labels,
        ads_types,
        data.box_dimensions)

    n_ads = get(raspa_out, "n_ads", NaN)

    @printf("  [ANALYSIS] N_ads=%.0f  ε_V=%.4e\n", n_ads, strain.εV)
    @printf("  [ANALYSIS] straight=%d  sinusoidal=%d  intersection=%d\n",
            occ["straight"], occ["sinusoidal"], occ["intersection"])

    summary = Dict{String,Any}(
        "iteration" => cycle, "pressure" => wp.pressure,
        "n_ads" => n_ads,
        "n_straight" => occ["straight"],
        "n_sinusoidal" => occ["sinusoidal"],
        "n_intersection" => occ["intersection"],
        "a" => cell.a, "b" => cell.b, "c" => cell.c,
        "alpha" => cell.alpha, "beta" => cell.beta,
        "gamma" => cell.gamma, "volume" => cell.volume,
        "strain_a" => strain.εa, "strain_b" => strain.εb,
        "strain_c" => strain.εc, "strain_V" => strain.εV,
    )

    write_cycle_summary(joinpath(wp.base_dir, "convergence.csv"), summary)
    return summary
end

# ════════════════════════════════════════════════════════════════
# Convergence check
# ════════════════════════════════════════════════════════════════

function check_convergence(history::Vector{Dict{String,Any}}, wp::WorkflowParams)
    n = wp.convergence_window
    length(history) < n && return false
    recent = history[end-n+1:end]
    for key in ["volume", "n_ads"]
        vals = Float64[r[key] for r in recent if !isnan(get(r, key, NaN))]
        length(vals) < n && continue
        mean(vals) > 0 || continue
        std(vals) / mean(vals) > wp.convergence_tol && return false
    end
    return true
end

# ════════════════════════════════════════════════════════════════
# Main workflow loop
# ════════════════════════════════════════════════════════════════

"""
    run_gcmc_md_workflow(wp::WorkflowParams)

Run the iterative GCMC/MD workflow. All files are created under
`wp.base_dir`.

Each cycle:
  1. RASPA3 GCMC → ethanol positions (JSON restart)
  2. Build/reload LAMMPS data file
  3. LAMMPS NPT-MD → relaxed structure
  4. Extract distorted CIF → feeds next GCMC
  5. Analysis → convergence.csv
"""
function run_gcmc_md_workflow(wp::WorkflowParams)
    mkpath(wp.base_dir)

    println("╔═══════════════════════════════════════════════════╗")
    println("║  JZeoMCMD — GCMC/NPT-MD Workflow                 ║")
    println("╚═══════════════════════════════════════════════════╝")
    @printf("  T = %.1f K   P = %.1e Pa\n", wp.temperature, wp.pressure)
    println("  Workdir:      $(wp.base_dir)")
    println("  RASPA inputs: $(wp.raspa_input_dir)")
    println("  LAMMPS input: $(wp.lammps_input)")
    println("  RASPA steps:  init=$(wp.raspa_n_init) equil=$(wp.raspa_n_equil) prod=$(wp.raspa_n_prod)")
    println("  LAMMPS steps: NVT=$(wp.lammps_nvt_ramp_steps)+$(wp.lammps_nvt_hold_steps) NPT=$(wp.lammps_npt_steps)")
    println("  Convergence:  CV < $(wp.convergence_tol*100)% over $(wp.convergence_window) cycles")
    println()

    # Ensure table exists
    table_path = ensure_tables!(wp)

    # Initial state
    current_cif  = wp.initial_cif
    current_data = wp.initial_data
    history = Dict{String,Any}[]

    # Reference cell for strain
    ref = read_lammps_data(current_data; verbose=false)
    bd = ref.box_dimensions
    ref_cell = (a=bd[1,2]-bd[1,1], b=bd[2,2]-bd[2,1], c=bd[3,2]-bd[3,1],
                alpha=90.0, beta=90.0, gamma=90.0,
                volume=prod([bd[d,2]-bd[d,1] for d in 1:3]))
    @printf("  Reference cell: a=%.3f b=%.3f c=%.3f V=%.1f\n\n",
            ref_cell.a, ref_cell.b, ref_cell.c, ref_cell.volume)

    for cycle in 1:wp.max_iterations
        println("═" ^ 60)
        @printf("  CYCLE %d / %d   \n", cycle, wp.max_iterations)
        println("═" ^ 60)

        # 1. GCMC
        println("\n  ── Step 1: RASPA3 GCMC ──")
        println("    CIF: $(basename(current_cif))")
        gcmc = step_gcmc!(wp, cycle, current_cif)
        if gcmc === nothing
            println("  ✗ GCMC failed — stopping.")
            break
        end

        # 2. Build/reload data file
        println("\n  ── Step 2: Build LAMMPS data file ──")
        data_file = step_build_data!(wp, cycle, gcmc.restart, current_data)

        # 3. NPT-MD
        println("\n  ── Step 3: LAMMPS NPT-MD ──")
        npt = step_npt!(wp, cycle, data_file)
        if npt === nothing
            println("  ✗ NPT failed — stopping.")
            break
        end

        # 4. Distorted CIF
        println("\n  ── Step 4: Extract distorted CIF ──")
        new_cif = step_write_cif!(wp, cycle, npt.data)

        # 5. Analysis
        println("\n  ── Step 5: Analysis ──")
        summary = step_analyze!(wp, cycle, npt.cell, gcmc.output,
                                 ref_cell, npt.data)
        push!(history, summary)

        # Update for next cycle
        current_cif  = new_cif      # next GCMC uses distorted CIF
        current_data = npt.data     # next reload uses NPT output

        # Convergence
        if check_convergence(history, wp)
            println("\n★ Converged after $cycle iterations ★")
            break
        end
        println()
    end

    println("\n" * "═" ^ 60)
    println("  WORKFLOW COMPLETE: $(length(history)) cycles")
    println("  Results: $(joinpath(wp.base_dir, "convergence.csv"))")
    println("═" ^ 60)
end

# ════════════════════════════════════════════════════════════════
# Pressure sweep setup
# ════════════════════════════════════════════════════════════════

"""
    setup_pressure_sweep(pressures, initial_data, initial_cif;
                          base_dir=".", temperature=300.0)

Create folder structure for multiple pressure points.
"""
function setup_pressure_sweep(pressures::Vector{Float64},
                               initial_data::String, initial_cif::String;
                               base_dir::String = ".",
                               temperature::Float64 = 300.0,
                               max_iterations::Int = 50)
    mkpath(base_dir)
    for p in pressures
        pname = @sprintf("P_%.1e", p)
        pdir = joinpath(base_dir, pname)
        mkpath(pdir)
        cp(initial_data, joinpath(pdir, basename(initial_data)); force=true)
        cp(initial_cif, joinpath(pdir, basename(initial_cif)); force=true)
        open(joinpath(pdir, "config.toml"), "w") do io
            println(io, "pressure = $p")
            println(io, "temperature = $temperature")
            println(io, "max_iterations = $max_iterations")
            println(io, "initial_data = \"$(basename(initial_data))\"")
            println(io, "initial_cif = \"$(basename(initial_cif))\"")
        end
        println("  Created $pname/")
    end
end
