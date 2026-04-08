# ════════════════════════════════════════════════════════════════
# workflow.jl — Master GCMC/MD loop (included by JZeoMCMD.jl)
#
# All simulation files go to wp.base_dir (user's project folder).
# Package resources (params.toml, raspa_inputs/) are found via
# JZeoMCMD.resource_dir().
# ════════════════════════════════════════════════════════════════


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

    # ── Package resources (defaults to package paths) ──
    params_toml::String     = params_path()         # from JZeoMCMD
    raspa_input_dir::String = raspa_inputs_path()    # from JZeoMCMD
    lammps_input::String    = ""                     # user must provide

    # ── Executables ──
    raspa_exe::String  = "raspa3"
    lammps_exe::String = "lmp"

    # ── Physics ──
    temperature::Float64 = 300.0
    pressure::Float64    = 1e5

    # ── Iteration control ──
    max_iterations::Int     = 50
    convergence_tol::Float64 = 0.005
    convergence_window::Int  = 5

    # ── Derived (set automatically) ──
    table_file::String = ""   # auto-generated if empty
    adsorbate_component::String = "ethanol"
end

"""
    load_workflow_params(toml_path; base_dir=pwd(), kwargs...)

Load workflow settings from params.toml. The `base_dir` is where
simulations will run (NOT the package directory).
"""
function load_workflow_params(toml_path::String; base_dir::String=pwd(), kwargs...)
    p = TOML.parsefile(toml_path)
    wf = get(p, "workflow", Dict())

    wp = WorkflowParams(
        base_dir        = base_dir,
        params_toml     = toml_path,
        raspa_input_dir = get(wf, "raspa_input_dir", raspa_inputs_path()),
        lammps_input    = get(wf, "lammps_input", ""),
        raspa_exe       = get(wf, "raspa_exe", "raspa3"),
        lammps_exe      = get(wf, "lammps_exe", "lmp"),
        table_file      = get(wf, "table_file", ""),
        adsorbate_component = get(wf, "adsorbate_component", "ethanol"),
    )

    for (k, v) in kwargs
        hasproperty(wp, k) && setproperty!(wp, k, v)
    end
    return wp
end

# ════════════════════════════════════════════════════════════════
# Auto-generate table in the user's base_dir
# ════════════════════════════════════════════════════════════════

function ensure_tables!(wp::WorkflowParams)
    # Check if table already exists in base_dir
    if !isempty(wp.table_file)
        full = joinpath(wp.base_dir, wp.table_file)
        isfile(full) && return full
    end

    # Auto-generate
    p = TOML.parsefile(wp.params_toml)
    nb = p["framework"]["nonbonded"]
    ps = p["pair_style"]

    A_Si = nb["A_Si"]; A_O = nb["A_O"]
    A_Al = get(nb, "A_Al", A_Si)
    rmin = ps["table_rmin"]; rmax = ps["table_rmax"]; N = ps["table_npoints"]
    al_type = get(p["framework"], "Al_type", 0)

    pairs = [("Si_Si",A_Si,A_Si), ("Si_O",A_Si,A_O), ("O_O",A_O,A_O)]
    if al_type > 0
        append!(pairs, [("Si_Al",A_Si,A_Al), ("O_Al",A_O,A_Al), ("Al_Al",A_Al,A_Al)])
    end

    fname = isempty(wp.table_file) ? "hillsauer_nb.table" : wp.table_file
    out = joinpath(wp.base_dir, fname)

    open(out, "w") do io
        println(io, "# Auto-generated Hill-Sauer A/r^9 table\n# E=A/r^9  F=9A/r^10\n# real units\n")
        for (label, Ai, Aj) in pairs
            A = sqrt(Ai * Aj)
            println(io, label); println(io, "N $N R $rmin $rmax\n")
            dr = (rmax - rmin) / (N - 1)
            for k in 1:N
                r = rmin + (k-1)*dr
                @printf(io, "%d  %.6f  %.6f  %.6f\n", k, r, A/r^9, 9A/r^10)
            end
            println(io)
        end
    end
    wp.table_file = fname
    println("  Auto-generated table: $out")
    return out
end

# ════════════════════════════════════════════════════════════════
# Folder management — all under base_dir
# ════════════════════════════════════════════════════════════════

function iter_dir(wp::WorkflowParams, n::Int)
    d = joinpath(wp.base_dir, @sprintf("iter_%03d", n))
    mkpath(d)
    return d
end

# ════════════════════════════════════════════════════════════════
# Step 1: GCMC — copy raspa_input_dir, patch pressure + CIF
# ════════════════════════════════════════════════════════════════

function step_gcmc!(wp::WorkflowParams, iteration::Int, cif_path::String)
    rdir = joinpath(iter_dir(wp, iteration), "raspa")
    mkpath(rdir)

    # Copy user's (or package's) RASPA input files
    src_dir = wp.raspa_input_dir
    if !isdir(src_dir)
        @warn "RASPA input dir not found: $src_dir"
        return nothing
    end
    for f in readdir(src_dir; join=true)
        isfile(f) && cp(f, joinpath(rdir, basename(f)); force=true)
    end

    # Copy CIF
    cp(cif_path, joinpath(rdir, basename(cif_path)); force=true)

    # Patch simulation.json: pressure + framework name
    local _JSON = Base.require(Base.PkgId(
        Base.UUID("682c06a0-de6a-54ab-a142-c8b1cf79cde6"), "JSON"))

    sim_json = joinpath(rdir, "simulation.json")
    if isfile(sim_json)
        content = _JSON.parsefile(sim_json)
        if haskey(content, "Systems") && !isempty(content["Systems"])
            content["Systems"][1]["ExternalPressure"] = wp.pressure
            content["Systems"][1]["Name"] = replace(basename(cif_path), ".cif"=>"")
        end
        open(sim_json, "w") do io; _JSON.print(io, content, 2); end
    else
        # Maybe it's a template — rename
        tmpl = joinpath(rdir, "simulation.json.template")
        if isfile(tmpl)
            txt = read(tmpl, String)
            txt = replace(txt,
                "__PRESSURE__"  => string(wp.pressure),
                "__FRAMEWORK__" => replace(basename(cif_path), ".cif"=>""),
                "__TEMPERATURE__" => string(wp.temperature))
            write(sim_json, txt)
        end
    end

    println("  [GCMC] Running RASPA3 (P=$(wp.pressure) Pa)...")
    try
        run(Cmd(`$(wp.raspa_exe)`, dir=rdir); wait=true)
    catch e
        @warn "RASPA3 failed" exception=e
        return nothing
    end

    # Find restart JSON
    restarts = filter(f -> startswith(f,"restart_") && endswith(f,".json"), readdir(rdir))
    isempty(restarts) && (@warn "No restart JSON"; return nothing)
    restart_path = joinpath(rdir, sort(restarts)[end])

    raspa_out = parse_raspa3_output(rdir)
    @printf("  [GCMC] N_ads = %.1f\n", raspa_out["n_ads"])

    return (restart=restart_path, output=raspa_out, dir=rdir)
end

# ════════════════════════════════════════════════════════════════
# Step 2: Reload adsorbate
# ════════════════════════════════════════════════════════════════

function step_reload!(wp::WorkflowParams, iteration::Int,
                       prev_data::String, restart_json::String)
    combined = joinpath(iter_dir(wp, iteration), "combined.lmp")
    println("  [RELOAD] Swapping adsorbate...")
    reload_adsorbate(prev_data, restart_json, combined)
    return combined
end

# ════════════════════════════════════════════════════════════════
# Step 3: LAMMPS NPT — copy user's input + data + table
# ════════════════════════════════════════════════════════════════

function step_npt!(wp::WorkflowParams, iteration::Int,
                    data_file::String, table_path::String)
    ldir = joinpath(iter_dir(wp, iteration), "lammps")
    mkpath(ldir)

    cp(data_file, joinpath(ldir, basename(data_file)); force=true)
    cp(table_path, joinpath(ldir, basename(table_path)); force=true)

    if isfile(wp.lammps_input)
        cp(wp.lammps_input, joinpath(ldir, basename(wp.lammps_input)); force=true)
    else
        @warn "LAMMPS input not found: $(wp.lammps_input)"
        return nothing
    end

    println("  [NPT] Running LAMMPS...")
    try
        run(Cmd(`$(wp.lammps_exe) -in $(basename(wp.lammps_input))`, dir=ldir); wait=true)
    catch e
        @warn "LAMMPS failed" exception=e
        return nothing
    end

    # Find output data
    outputs = filter(f -> (endswith(f,".data")||endswith(f,".lmp")) &&
                          f != basename(data_file), readdir(ldir))
    isempty(outputs) && (@warn "No LAMMPS output"; return nothing)
    output_data = joinpath(ldir, outputs[end])

    log_file = joinpath(ldir, "log.lammps")
    log_data = isfile(log_file) ? parse_lammps_log(log_file) : Dict{String,Vector{Float64}}()
    cell = extract_cell_params(log_data)
    @printf("  [NPT] a=%.3f b=%.3f c=%.3f V=%.1f\n", cell.a, cell.b, cell.c, cell.volume)

    return (data=output_data, cell=cell, log=log_data, dir=ldir)
end

# ════════════════════════════════════════════════════════════════
# Step 4: CIF
# ════════════════════════════════════════════════════════════════

function step_write_cif!(wp::WorkflowParams, iteration::Int, data_file::String)
    cif_out = joinpath(iter_dir(wp, iteration), "framework_relaxed.cif")
    data = LammpsDataReader.read_lammps_data(data_file; verbose=false)

    fw_types = [1, 2]
    p = TOML.parsefile(wp.params_toml)
    al = get(p["framework"], "Al_type", 0)
    al > 0 && push!(fw_types, al)

    write_cif(cif_out, data;
              framework_types=fw_types,
              type_elements=Dict(1=>"Si",2=>"O",3=>"Al"),
              comment=@sprintf("cycle %d P=%.1e Pa", iteration, wp.pressure))
    return cif_out
end

# ════════════════════════════════════════════════════════════════
# Step 5: Analysis
# ════════════════════════════════════════════════════════════════

function step_analyze!(wp::WorkflowParams, iteration::Int,
                        npt_result, gcmc_out::Dict,
                        ref_cell::NamedTuple, data_file::String)
    strain = compute_strain(npt_result.cell, ref_cell)

    data = LammpsDataReader.read_lammps_data(data_file; verbose=false)
    occ = analyze_channel_occupancy(data.coords, data.atom_labels,
                                     Set([3,4,5,6]), data.box_dimensions)
    n_ads = get(gcmc_out, "n_ads", NaN)

    @printf("  [ANALYSIS] N=%d ε_V=%.2e str=%d sin=%d int=%d\n",
            round(Int, n_ads), strain.εV,
            occ["straight"], occ["sinusoidal"], occ["intersection"])

    summary = Dict{String,Any}(
        "iteration"=>iteration, "pressure"=>wp.pressure, "n_ads"=>n_ads,
        "n_straight"=>occ["straight"], "n_sinusoidal"=>occ["sinusoidal"],
        "n_intersection"=>occ["intersection"],
        "a"=>npt_result.cell.a, "b"=>npt_result.cell.b, "c"=>npt_result.cell.c,
        "alpha"=>npt_result.cell.alpha, "beta"=>npt_result.cell.beta,
        "gamma"=>npt_result.cell.gamma, "volume"=>npt_result.cell.volume,
        "strain_a"=>strain.εa, "strain_b"=>strain.εb,
        "strain_c"=>strain.εc, "strain_V"=>strain.εV)

    write_cycle_summary(joinpath(wp.base_dir, "convergence.csv"), summary)
    return summary
end

# ════════════════════════════════════════════════════════════════
# Convergence
# ════════════════════════════════════════════════════════════════

function check_convergence(history::Vector{Dict{String,Any}}, wp::WorkflowParams)
    n = wp.convergence_window
    length(history) < n && return false
    recent = history[end-n+1:end]
    for key in ["volume", "n_ads"]
        vals = Float64[r[key] for r in recent if !isnan(get(r,key,NaN))]
        length(vals) < n && continue
        mean(vals) > 0 || continue
        std(vals)/mean(vals) > wp.convergence_tol && return false
    end
    return true
end

# ════════════════════════════════════════════════════════════════
# Main loop
# ════════════════════════════════════════════════════════════════

"""
    run_gcmc_md_workflow(wp::WorkflowParams)

Run the iterative GCMC/MD workflow. All files are created under
`wp.base_dir`, which should be the user's project directory.
"""
function run_gcmc_md_workflow(wp::WorkflowParams)
    mkpath(wp.base_dir)

    println("╔═══════════════════════════════════════════════════╗")
    println("║  JZeoMCMD — GCMC/NPT-MD Workflow                 ║")
    println("╚═══════════════════════════════════════════════════╝")
    @printf("  T = %.1f K   P = %.1e Pa\n", wp.temperature, wp.pressure)
    println("  Working dir:  $(wp.base_dir)")
    println("  RASPA inputs: $(wp.raspa_input_dir)")
    println("  LAMMPS input: $(wp.lammps_input)")
    println("  Params:       $(wp.params_toml)\n")

    table_path = ensure_tables!(wp)

    current_cif  = wp.initial_cif
    current_data = wp.initial_data
    history = Dict{String,Any}[]

    # Reference cell
    ref = LammpsDataReader.read_lammps_data(current_data; verbose=false)
    bd = ref.box_dimensions
    ref_cell = (a=bd[1,2]-bd[1,1], b=bd[2,2]-bd[2,1], c=bd[3,2]-bd[3,1],
                alpha=90.0, beta=90.0, gamma=90.0,
                volume=prod([bd[d,2]-bd[d,1] for d in 1:3]))

    for iter in 1:wp.max_iterations
        println("\n═══ Iteration $iter / $(wp.max_iterations) ═══")

        gcmc = step_gcmc!(wp, iter, current_cif)
        gcmc_out = gcmc !== nothing ? gcmc.output : Dict{String,Any}("n_ads"=>NaN)

        if gcmc !== nothing
            combined = step_reload!(wp, iter, current_data, gcmc.restart)
        else
            combined = joinpath(iter_dir(wp, iter), "combined.lmp")
            cp(current_data, combined; force=true)
        end

        npt = step_npt!(wp, iter, combined, table_path)
        npt === nothing && (println("  ✗ NPT failed"); break)

        new_cif = step_write_cif!(wp, iter, npt.data)
        summary = step_analyze!(wp, iter, npt, gcmc_out, ref_cell, npt.data)
        push!(history, summary)

        current_cif  = new_cif
        current_data = npt.data

        if check_convergence(history, wp)
            println("\n★ Converged after $iter iterations ★")
            break
        end
    end

    println("\n═══ Done: $(length(history)) iterations ═══")
    println("  $(joinpath(wp.base_dir, "convergence.csv"))")
end

# ════════════════════════════════════════════════════════════════
# Pressure sweep setup
# ════════════════════════════════════════════════════════════════

"""
    setup_pressure_sweep(pressures, initial_data, initial_cif;
                          base_dir=".", temperature=300.0)

Create folder structure for multiple pressure points.
Each gets its own config.toml that can be passed to the workflow.
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
