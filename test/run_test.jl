#!/usr/bin/env julia
"""
run_test.jl — Full GCMC/MD test for MFI + ethanol

Run from any directory:
    julia ~/.julia/dev/JZeoMCMD/test/run_test.jl \\
        --workdir ~/simulations/MFI_test \\
        --pressure 1e6 \\
        --ncycles 3

Requires: LAMMPS (lmp) and RASPA3 (raspa3) in PATH.
"""

using Pkg
Pkg.activate(joinpath(homedir(), ".julia", "dev", "JZeoMCMD"))
using JZeoMCMD
using Printf
using Dates

import TOML
import JSON

# ════════════════════════════════════════════════════════════════
# Configuration
# ════════════════════════════════════════════════════════════════

Base.@kwdef mutable struct TestConfig
    workdir::String      = joinpath(homedir(), "simulations", "MFI_ethanol_test")
    pressure::Float64    = 1e6          # Pa
    temperature::Float64 = 300.0        # K
    ncycles::Int         = 3            # GCMC/MD iterations
    lammps_exe::String   = "lmp_mpi"
    raspa_exe::String    = "raspa3"
    # LAMMPS NPT settings
    npt_steps::Int       = 500000       # production steps
    timestep::Float64    = 0.25         # fs
end

function parse_test_args(args)
    tc = TestConfig()
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--workdir"     && i<length(args); tc.workdir = args[i+1]; i+=2
        elseif a == "--pressure"    && i<length(args); tc.pressure = parse(Float64,args[i+1]); i+=2
        elseif a == "--temperature" && i<length(args); tc.temperature = parse(Float64,args[i+1]); i+=2
        elseif a == "--ncycles"     && i<length(args); tc.ncycles = parse(Int,args[i+1]); i+=2
        elseif a == "--lammps"      && i<length(args); tc.lammps_exe = args[i+1]; i+=2
        elseif a == "--raspa"       && i<length(args); tc.raspa_exe = args[i+1]; i+=2
        elseif a == "--npt-steps"   && i<length(args); tc.npt_steps = parse(Int,args[i+1]); i+=2
        elseif a in ["--help", "-h"]
            println("""
            Usage: julia run_test.jl [options]

              --workdir DIR       Simulation directory (default: ~/simulations/MFI_ethanol_test)
              --pressure P        Pressure in Pa (default: 1e6)
              --temperature T     Temperature in K (default: 300)
              --ncycles N         Number of GCMC/MD cycles (default: 3)
              --npt-steps N       LAMMPS NPT production steps (default: 500000)
              --lammps EXE        LAMMPS executable (default: lmp)
              --raspa EXE         RASPA3 executable (default: raspa3)
            """)
            exit(0)
        else
            i += 1
        end
    end
    return tc
end

# ════════════════════════════════════════════════════════════════
# Helper: run RASPA3 GCMC
# ════════════════════════════════════════════════════════════════

function run_raspa_gcmc(tc::TestConfig, cycle::Int, cif_path::String)
    rdir = joinpath(tc.workdir, @sprintf("cycle_%02d", cycle), "raspa")
    mkpath(rdir)

    # Copy RASPA3 input files from package
    pkg_raspa = JZeoMCMD.raspa_inputs_path()
    for f in readdir(pkg_raspa; join=true)
        isfile(f) && cp(f, joinpath(rdir, basename(f)); force=true)
    end

    # Copy CIF into raspa directory
    cif_name = basename(cif_path)
    cp(cif_path, joinpath(rdir, cif_name); force=true)
    fw_name = replace(cif_name, ".cif" => "")

    # Handle simulation.json
    if cycle == 1 # the first gcmc will use the unicell
    	sim_file = joinpath(rdir, "simulation.json")
    	sim_template = joinpath(rdir, "simulation.json.template")
    else # when cycle is greater than 1 we will use the supercell extracted from lammps simulation, thus raspa will use 1x1x1 - because the structure alreasy has a 2x2x2 uc size!
	sim_file = joinpath(rdir, "simulation.json")
	sim_template = joinpath(rdir, "simulation.json.template_next")
    end

    if isfile(sim_template)
        # Fill template
        txt = read(sim_template, String)
        txt = replace(txt,
            "__FRAMEWORK__"   => fw_name,
            "__TEMPERATURE__" => string(tc.temperature),
            "__PRESSURE__"    => string(tc.pressure),
        )
        write(sim_file, txt)
        rm(sim_template)  # clean up template from run dir
    elseif isfile(sim_file)
        # Patch existing JSON
        content = JSON.parsefile(sim_file)
        if haskey(content, "Systems") && !isempty(content["Systems"])
            content["Systems"][1]["ExternalPressure"] = tc.pressure
            content["Systems"][1]["ExternalTemperature"] = tc.temperature
            content["Systems"][1]["Name"] = fw_name
        end
        open(sim_file, "w") do io
            JSON.print(io, content, 2)
        end
    else
        error("No simulation.json or simulation.json.template found")
    end

    # Run RASPA3
    println("    Running RASPA3 in $rdir ...")
    t0 = time()
    run(Cmd(`$(tc.raspa_exe)`, dir=rdir))
    @printf("    RASPA3 done (%.1f s)\n", time()-t0)

    # Find restart JSON
    #  the output and restart file will be inside the output folder created by raspa3
    output_dir = joinpath(rdir, "output")
    # ideally it will be only one, since each bach of cycles will be responsible to simulate one pressue and one temperature
    restarts = sort(filter(f -> startswith(f,"restart_") && endswith(f,".json"),
                           readdir(output_dir)))
    if isempty(restarts)
        error("No RASPA3 restart JSON found in $rdir")
    end
    restart_path = joinpath(output_dir, restarts[end])

    # Parse N_ads directly from the restart JSON file
    
    raspa_out = parse_raspa3_output(restart_path)
    @printf("    N_ads = %.1f molecules\n", raspa_out["n_ads"])
   # raspa_out = Dict{String,Any}()
   # restart_data = JSON.parsefile(restart_path)
   # for (k, v) in restart_data
   # 	if v isa AbstractVector && !isempty(v) && v[1] isa AbstractVector
   #     	raspa_out["n_ads"] = length(v) / tc.atoms_per_mol  # 4 for ethanol
   #     	break
   # 	end
    #end
    
   # @printf("    N_ads = %.0f molecules\n", raspa_out["n_ads"])

    return restart_path, raspa_out
end

# ════════════════════════════════════════════════════════════════
# Helper: run LAMMPS NPT-MD
# ════════════════════════════════════════════════════════════════

function run_lammps_npt(tc::TestConfig, cycle::Int, data_file::String)
    ldir = joinpath(tc.workdir, @sprintf("cycle_%02d", cycle), "lammps")
    mkpath(ldir)

    # Copy data file
    data_name = "loaded.lmp"
    cp(data_file, joinpath(ldir, data_name); force=true)

    # Copy table file
    table_src = joinpath(tc.workdir, "hillsauer_nb.table")
    table_name = basename(table_src)
    cp(table_src, joinpath(ldir, table_name); force=true)

    # Copy LAMMPS input script
    input_src = joinpath(tc.workdir, "run_npt.in")
    cp(input_src, joinpath(ldir, "run_npt.in"); force=true)

    # Run LAMMPS
    println("    Running LAMMPS in $ldir ...")
    t0 = time()
    run(Cmd(`$(tc.lammps_exe) -in run_npt.in`, dir=ldir))
    @printf("    LAMMPS done (%.1f s)\n", time()-t0)

    # Find output data
    output_data = joinpath(ldir,"npt_final.data")
    if !isfile(output_data)
        # Try other common names from write_data
        for candidate in readdir(ldir)
		#if endswith(candidate, (".data", ".lmp")) && candidate != data_name
	    if (endswith(candidate, ".data") || endswith(candidate, ".lmp")) && candidate != data_name
		output_data = joinpath(ldir, candidate)
                break
            end
        end
    end
    !isfile(output_data) && error("No LAMMPS output data found in $ldir")

    # Parse log
    log_file = joinpath(ldir, "log.lammps")
    log_data = isfile(log_file) ? parse_lammps_log(log_file) : Dict{String,Vector{Float64}}()
    cell = extract_cell_params(log_data)

    return output_data, cell
end

# ════════════════════════════════════════════════════════════════
# Main test
# ════════════════════════════════════════════════════════════════

function run_test(tc::TestConfig)
    println("╔════════════════════════════════════════════════════╗")
    println("║  JZeoMCMD — GCMC/NPT-MD Test                     ║")
    println("╚════════════════════════════════════════════════════╝")
    @printf("  Workdir:     %s\n", tc.workdir)
    @printf("  Pressure:    %.1e Pa\n", tc.pressure)
    @printf("  Temperature: %.1f K\n", tc.temperature)
    @printf("  Cycles:      %d\n", tc.ncycles)
    println()

    mkpath(tc.workdir)
    src = JZeoMCMD.structures_path()

    # ──────────────────────────────────────────────────────────
    # PREPARATION: Generate table file in workdir (if missing)
    # ──────────────────────────────────────────────────────────
    table_path = joinpath(tc.workdir, "hillsauer_nb.table")
    if !isfile(table_path)
        println("Generating non-bonded table...")
        # Copy from package ff/ or auto-generate
        pkg_table = joinpath(JZeoMCMD.resource_dir(), "ff", "hillsauer_silica.table")
        if isfile(pkg_table)
            cp(pkg_table, table_path; force=true)
            println("  Copied from package: $table_path")
        else
            # Auto-generate from params.toml
            p = TOML.parsefile(JZeoMCMD.params_path())
            nb = p["framework"]["nonbonded"]
            ps = p["pair_style"]
            A_Si = nb["A_Si"]; A_O = nb["A_O"]
            rmin = ps["table_rmin"]; rmax = ps["table_rmax"]; N = ps["table_npoints"]
            pairs = [("Si_Si",A_Si,A_Si), ("Si_O",A_Si,A_O), ("O_O",A_O,A_O)]
            open(table_path, "w") do io
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
            println("  Auto-generated: $table_path")
        end
    end

    # ──────────────────────────────────────────────────────────
    # PREPARATION: Copy initial CIF to workdir
    # ──────────────────────────────────────────────────────────
    initial_cif = joinpath(tc.workdir, "MFI_SI.cif")
    if !isfile(initial_cif)
        cp(joinpath(src, "MFI_SI.cif"), initial_cif; force=true)
        println("Copied initial CIF: $initial_cif")
    end

    # ──────────────────────────────────────────────────────────
    # PREPARATION: Build LAMMPS input script template in workdir
    # (only needs to be done once — user can edit afterwards)
    # ──────────────────────────────────────────────────────────
    lammps_input_path = joinpath(tc.workdir, "run_npt.in")
    if !isfile(lammps_input_path)
        println("Generating LAMMPS input script...")
        # Use the write_input_script from build_loaded_zeolite
        # but we write it to workdir, referencing generic filenames
        cfg_tmp = ZeoliteConfig(pair_cutoff=12.0, coul_cutoff=12.0,
                                 table_file="hillsauer_nb.table")
        write_input_script(lammps_input_path, "loaded.lmp", cfg_tmp)
        println("  Wrote: $lammps_input_path")
        println("  ⚠  Review and edit this file before production runs!")
    end

    # Reference cell (for strain calculation)
    ref_data = read_lammps_data(joinpath(src, "MFI_SI.data"); verbose=false)
    bd = ref_data.box_dimensions
    ref_cell = (a=bd[1,2]-bd[1,1], b=bd[2,2]-bd[2,1], c=bd[3,2]-bd[3,1],
                alpha=90.0, beta=90.0, gamma=90.0,
                volume=prod([bd[d,2]-bd[d,1] for d in 1:3]))
    @printf("  Reference cell: a=%.3f b=%.3f c=%.3f V=%.1f\n\n",
            ref_cell.a, ref_cell.b, ref_cell.c, ref_cell.volume)

    # ──────────────────────────────────────────────────────────
    # Tracking CSV
    # ──────────────────────────────────────────────────────────
    csv_path = joinpath(tc.workdir, "convergence.csv")

    # State variables for the loop
    current_cif  = initial_cif
    current_data = ""  # set after cycle 1

    for cycle in 1:tc.ncycles
        println("═" ^ 60)
        @printf("  CYCLE %d / %d    [%s]\n", cycle, tc.ncycles,
                Dates.format(now(), "HH:MM:SS"))
        println("═" ^ 60)

        cycle_dir = joinpath(tc.workdir, @sprintf("cycle_%02d", cycle))
        mkpath(cycle_dir)

        # ══════════════════════════════════════════════════════
        # STEP 1: RASPA3 GCMC with current CIF
        # ══════════════════════════════════════════════════════
        println("\n  ── Step 1: RASPA3 GCMC ──")
        println("    CIF: $(basename(current_cif))")

        restart_json, raspa_out = run_raspa_gcmc(tc, cycle, current_cif)

        # ══════════════════════════════════════════════════════
        # STEP 2: Build / reload LAMMPS data file
        # ══════════════════════════════════════════════════════
        println("\n  ── Step 2: Build LAMMPS data file ──")
        data_file = joinpath(cycle_dir, "loaded.lmp")

        if cycle == 1
            # First cycle: build from Ovito data + RASPA3 ethanol
            println("    Building from Ovito data + RASPA3 ethanol...")
            cfg = ZeoliteConfig(
                ovito_data    = joinpath(src, "MFI_SI.data"),
                raspa_restart = restart_json,
                output_data   = data_file,
            )
            fw = read_and_remap_framework(cfg)
            add_framework_topology!(fw, cfg)
            mols = read_raspa3_ethanol(cfg, fw.box_dimensions)
            merge_framework_ethanol!(fw, mols, cfg)
            write_complete_data(data_file, fw, cfg)
        else
            # Subsequent cycles: reload adsorbate into previous NPT output
            println("    Reloading adsorbate into previous NPT framework...")
            reload_adsorbate(current_data, restart_json, data_file)
        end

        nfw = 2304  # MFI 2×2×2
        data_check = read_lammps_data(data_file; verbose=false)
        neth = size(data_check.coords, 1) - nfw
        @printf("    Atoms: %d (%d fw + %d ethanol = %d molecules)\n",
                size(data_check.coords,1), nfw, neth, neth÷4)

        # ══════════════════════════════════════════════════════
        # STEP 3: LAMMPS NPT-MD
        # ══════════════════════════════════════════════════════
        println("\n  ── Step 3: LAMMPS NPT-MD ──")
        npt_output, cell = run_lammps_npt(tc, cycle, data_file)
        @printf("    Cell: a=%.4f b=%.4f c=%.4f\n", cell.a, cell.b, cell.c)
        @printf("    Angles: α=%.3f° β=%.3f° γ=%.3f°\n", cell.alpha, cell.beta, cell.gamma)
        @printf("    Volume: %.2f ų\n", cell.volume)

        # ══════════════════════════════════════════════════════
        # STEP 4: Extract distorted CIF
        # ══════════════════════════════════════════════════════
        println("\n  ── Step 4: Extract distorted CIF ──")
        distorted_cif = joinpath(cycle_dir, "distorted.cif")
        npt_data = read_lammps_data(npt_output; verbose=false)
        write_cif(distorted_cif, npt_data;
                  framework_types=[1,2],
                  type_elements=Dict(1=>"Si", 2=>"O"),
                  comment=@sprintf("MFI cycle %d P=%.1e Pa", cycle, tc.pressure))

        # ══════════════════════════════════════════════════════
        # STEP 5: Analysis
        # ══════════════════════════════════════════════════════
        println("\n  ── Step 5: Analysis ──")
        strain = compute_strain(cell, ref_cell)
        @printf("    Strain: ε_a=%.4e  ε_b=%.4e  ε_c=%.4e  ε_V=%.4e\n",
                strain.εa, strain.εb, strain.εc, strain.εV)

        occ = analyze_channel_occupancy(
            npt_data.coords, npt_data.atom_labels,
            Set([3,4,5,6]), npt_data.box_dimensions)
        @printf("    Channels: straight=%d  sinusoidal=%d  intersection=%d\n",
                occ["straight"], occ["sinusoidal"], occ["intersection"])

        # Write to convergence CSV
        summary = Dict{String,Any}(
            "iteration" => cycle, "pressure" => tc.pressure,
            "n_ads" => get(raspa_out, "n_ads", NaN),
            "n_straight" => occ["straight"],
            "n_sinusoidal" => occ["sinusoidal"],
            "n_intersection" => occ["intersection"],
            "a" => cell.a, "b" => cell.b, "c" => cell.c,
            "alpha" => cell.alpha, "beta" => cell.beta, "gamma" => cell.gamma,
            "volume" => cell.volume,
            "strain_a" => strain.εa, "strain_b" => strain.εb,
            "strain_c" => strain.εc, "strain_V" => strain.εV,
        )
        write_cycle_summary(csv_path, summary)

        # ══════════════════════════════════════════════════════
        # Update state for next cycle
        # ══════════════════════════════════════════════════════
        current_cif  = distorted_cif   # next GCMC uses the distorted structure
        current_data = npt_output       # next reload uses the NPT output

        println()
    end

    # ══════════════════════════════════════════════════════════
    # Final summary
    # ══════════════════════════════════════════════════════════
    println("═" ^ 60)
    println("  TEST COMPLETE: $(tc.ncycles) cycles")
    println("  Results: $csv_path")
    println()
    println("  Files per cycle:")
    for c in 1:tc.ncycles
        cdir = joinpath(tc.workdir, @sprintf("cycle_%02d", c))
        println("    cycle_$(lpad(c,2,'0'))/")
        println("      raspa/    → GCMC output + restart JSON")
        println("      lammps/   → NPT output + trajectory")
        println("      distorted.cif → framework for next GCMC")
        println("      loaded.lmp    → input data file")
    end
    println("═" ^ 60)
end

# ════════════════════════════════════════════════════════════════
# Entry point
# ════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    tc = parse_test_args(ARGS)
    run_test(tc)
end
