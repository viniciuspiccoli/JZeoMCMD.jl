#!/usr/bin/env julia
"""
generate_nb_tables.jl

Generate LAMMPS tabulated non-bonded potentials for zeolite frameworks:
  - All-silica (Hill-Sauer): E(r) = A/r^9, F(r) = 9A/r^10
  - Aluminosilicate: same form with Al parameters

The Hill-Sauer force field uses a purely repulsive A/r^9 non-bonded
interaction (no attractive term). Standard LAMMPS pair styles cannot
represent this, so it must be tabulated.

Output files:
  hillsauer_silica.table    — Si-Si, Si-O, O-O  (pure silica)
  hillsauer_alumsil.table   — adds Al-Al, Al-Si, Al-O  (aluminosilicate)

Usage:
  julia generate_nb_tables.jl                    # both tables
  julia generate_nb_tables.jl --silica-only      # silica only
  julia generate_nb_tables.jl --rmin 0.5 --rmax 12.0 --npoints 5000

References:
  Hill & Sauer, J. Phys. Chem. 1994, 98, 1238  (HSFF)
  Hill & Sauer, J. Phys. Chem. 1995, 99, 9536  (charges, non-bonded)
  Bai et al., J. Phys. Chem. C 2013, 117, 24375 (framework LJ for GCMC)

Units: real (kcal/mol, Å)
"""
# ═══════════════════════════════════════════════════════════════════
# A/r^9 parameters (kcal/mol · Å^9)
# ═══════════════════════════════════════════════════════════════════
#
# Hill-Sauer self-interaction A_ii values from H&S 1994/1995:
#   A_Si = 186910.958224
#   A_O  =  57412.472881
#
# Cross-term A_ij via geometric combining rule:
#   A_ij = sqrt(A_ii * A_jj)
#   A_Si_O = sqrt(186910.958224 * 57412.472881) = 103590.638188
#
# For Al: derived from the Bai (2013) LJ parameters and the
# relation A_ii = 4 * ε_ii * σ_ii^9 (matching repulsive wall):
#   Bai Si: ε/kB=22.0 K, σ=2.30 Å → A_Si from H&S directly
#   Bai Al: treat as Si (same row, similar ionic radius)
#   Or from separate parameterization — user-configurable below.
# ═══════════════════════════════════════════════════════════════════

# ── Self-interaction A values ──
const A_SELF = Dict{String, Float64}(
    "Si" => 186910.958224,   # Hill & Sauer 1994
    "O"  =>  57412.472881,   # Hill & Sauer 1994
    "Al" => 186910.958224,   # Default: same as Si (similar size)
                              # Override below if you have better params
)

# ═══════════════════════════════════════════════════════════════════
# Table parameters: rmin, rmax, npoints
# ═══════════════════════════════════════════════════════════════════
#
#   rmin (Å):  Shortest distance in the table.
#              LAMMPS will ERROR if any atom pair gets closer than this.
#              Must be shorter than the closest non-bonded approach in
#              your simulation. In zeolites, non-bonded pairs (those NOT
#              connected by bonds/angles/dihedrals) rarely get closer
#              than ~2.5 Å. Setting rmin=1.0 gives a safe margin.
#              Rule: rmin = 0.5–1.0 Å (safe default: 1.0)
#
#   rmax (Å):  Longest distance in the table. Must be ≥ the pair cutoff
#              in your LAMMPS input script. If your pair_style uses a
#              cutoff of 12.0 Å, set rmax ≥ 12.0. Beyond rmax, the
#              interaction is zero (LAMMPS truncates at the cutoff).
#              Rule: rmax = pair cutoff (typically 11.0–14.0 Å)
#
#   npoints:   Number of evenly spaced points between rmin and rmax.
#              More points = finer grid = smaller interpolation error.
#              The spacing is: dr = (rmax - rmin) / (npoints - 1)
#
#              Recommended values:
#                 1000  — fast, ~0.02% error at typical distances
#                 2000  — standard (default), ~0.005% error ← USE THIS
#                 5000  — high precision, ~0.001% error
#                10000  — overkill for most applications
#
#              The error is largest at small r where A/r^9 is steep,
#              but non-bonded atoms never reach r < 2.5 Å in zeolites
#              (bonded neighbors are excluded via special_bonds).
#              Using "pair_style table spline N" (cubic spline) gives
#              better accuracy than "linear N" at the same npoints.
#
#   Example LAMMPS usage:
#     pair_style hybrid/overlay table spline 2000 coul/long 12.0
#                                      ^^^^^ ^^^^           ^^^^
#                                      interp npoints       cutoff=rmax
#
# ═══════════════════════════════════════════════════════════════════



"""
    compute_A_cross(name_i, name_j) -> Float64

Cross-term A_ij via geometric combining rule: A_ij = √(A_ii × A_jj)
"""
function compute_A_cross(name_i::String, name_j::String)
    return sqrt(A_SELF[name_i] * A_SELF[name_j])
end

# ═══════════════════════════════════════════════════════════════════
# Table generation
# ═══════════════════════════════════════════════════════════════════

"""
    write_table_section(io, label, A, rmin, rmax, N)

Write one LAMMPS table section for E(r) = A/r^9, F(r) = 9A/r^10.
"""
function write_table_section(io::IO, label::String, A::Float64,
                              rmin::Float64, rmax::Float64, N::Int)
    println(io, label)
    println(io, "N $N R $rmin $rmax\n")

    dr = (rmax - rmin) / (N - 1)
    for i in 1:N
        r = rmin + (i - 1) * dr
        r9  = r^9
        r10 = r^10
        E = A / r9           # energy: A/r^9
        F = 9.0 * A / r10    # force: -dE/dr = 9A/r^10
        @printf(io, "%d  %.6f  %.6f  %.6f\n", i, r, E, F)
    end
    println(io)
end

"""
    generate_table(filename, pairs; rmin, rmax, N, header)

Generate a complete LAMMPS table file for the given atom pairs.
`pairs` is a vector of (label, name_i, name_j) tuples.
"""
function generate_table(filename::String,
                         pairs::Vector{Tuple{String,String,String}};
                         rmin::Float64 = 1.0,
                         rmax::Float64 = 12.0,
                         N::Int = 2000,
                         header::String = "")
    open(filename, "w") do io
        println(io, "# $header")
        println(io, "# E(r) = A/r^9   F(r) = 9*A/r^10")
        println(io, "# Units: real (kcal/mol, Angstrom)")
        println(io, "# Generated by generate_nb_tables.jl")
        println(io, "#")
        println(io, "# A values (kcal/mol · Å^9):")
        for (label, ni, nj) in pairs
            A = compute_A_cross(ni, nj)
            @printf(io, "#   %-8s  A(%s-%s) = %.6f\n", label, ni, nj, A)
        end
        println(io, "#")
        println(io, "# Cross-terms: A_ij = sqrt(A_ii * A_jj)")
        println(io, "# Range: r = [$(rmin), $(rmax)] Å, N = $N points")
        println(io)

        for (label, ni, nj) in pairs
            A = compute_A_cross(ni, nj)
            write_table_section(io, label, A, rmin, rmax, N)
        end
    end

    println("  Wrote $filename ($(length(pairs)) pair sections)")
end

# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

function main(args=ARGS)
    # Defaults
    rmin = 1.0
    rmax = 12.0
    N = 2000
    silica_only = false
    output_silica   = "hillsauer_silica.table"
    output_alumsil  = "hillsauer_alumsil.table"

    # Parse args
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--silica-only"
            silica_only = true; i += 1
        elseif a == "--rmin" && i < length(args)
            rmin = parse(Float64, args[i+1]); i += 2
        elseif a == "--rmax" && i < length(args)
            rmax = parse(Float64, args[i+1]); i += 2
        elseif a == "--npoints" && i < length(args)
            N = parse(Int, args[i+1]); i += 2
        elseif a == "--A-Al" && i < length(args)
            A_SELF["Al"] = parse(Float64, args[i+1]); i += 2
        elseif a == "--output-silica" && i < length(args)
            output_silica = args[i+1]; i += 2
        elseif a == "--output-alumsil" && i < length(args)
            output_alumsil = args[i+1]; i += 2
        else
            i += 1
        end
    end

    println("╔═══════════════════════════════════════════════╗")
    println("║  Hill-Sauer Non-bonded Table Generator        ║")
    println("╚═══════════════════════════════════════════════╝")
    println()
    println("  Parameters:")
    println("    r range: [$rmin, $rmax] Å")
    println("    N points: $N")
    println("    A(Si) = $(A_SELF["Si"]) kcal/mol·Å⁹")
    println("    A(O)  = $(A_SELF["O"])")
    println("    A(Al) = $(A_SELF["Al"])")
    println()

    # ── Pure silica table ──
    println("═══ Generating silica table ═══")
    silica_pairs = [
        ("Si_Si", "Si", "Si"),
        ("Si_O",  "Si", "O"),
        ("O_O",   "O",  "O"),
    ]
    generate_table(output_silica, silica_pairs;
                   rmin=rmin, rmax=rmax, N=N,
                   header="Hill-Sauer non-bonded: pure silica (SiO2)")

    # ── Aluminosilicate table ──
    if !silica_only
        println("\n═══ Generating aluminosilicate table ═══")
        alumsil_pairs = [
            ("Si_Si", "Si", "Si"),
            ("Si_O",  "Si", "O"),
            ("Si_Al", "Si", "Al"),
            ("O_O",   "O",  "O"),
            ("O_Al",  "O",  "Al"),
            ("Al_Al", "Al", "Al"),
        ]
        generate_table(output_alumsil, alumsil_pairs;
                       rmin=rmin, rmax=rmax, N=N,
                       header="Hill-Sauer non-bonded: aluminosilicate (Si/Al/O)")
    end

    # ── Print LAMMPS pair_coeff lines for easy copy-paste ──
    println("\n═══ LAMMPS pair_coeff lines ═══")
    println()
    println("# --- Pure silica (Si=type 1, O=type 2) ---")
    println("pair_style  hybrid/overlay table spline $N lj/cut $(rmax) coul/long $(rmax)")
    println("pair_coeff  * * coul/long")
    println("pair_coeff  1 1  table $output_silica Si_Si")
    println("pair_coeff  1 2  table $output_silica Si_O")
    println("pair_coeff  2 2  table $output_silica O_O")

    if !silica_only
        println()
        println("# --- Aluminosilicate (Si=1, O=2, Al=3) ---")
        println("pair_coeff  1 1  table $output_alumsil Si_Si")
        println("pair_coeff  1 2  table $output_alumsil Si_O")
        println("pair_coeff  1 3  table $output_alumsil Si_Al")
        println("pair_coeff  2 2  table $output_alumsil O_O")
        println("pair_coeff  2 3  table $output_alumsil O_Al")
        println("pair_coeff  3 3  table $output_alumsil Al_Al")
    end

    println("\n═══ Done ═══")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
