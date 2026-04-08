# ═══════════════════════════════════════════════════════════════════
# write_cif.jl — Convert LAMMPS data to CIF for RASPA3
#
# Add to your workflow with:
#   include("read_lammps_data.jl")
#   include("write_cif.jl")
#   data = read_lammps_data("loaded_npt_final.lmp")
#   write_cif("relaxed_framework.cif", data; framework_types=[1,2])
# ═══════════════════════════════════════════════════════════════════

"""
    write_cif(filename, data; framework_types, type_elements, space_group, comment)

Write a CIF file from a LammpsData object.

# How it works — step by step:

## 1. Cell matrix → cell parameters (a, b, c, α, β, γ)
   LAMMPS stores the box as:
     - Orthogonal: xlo, xhi, ylo, yhi, zlo, zhi
     - Triclinic: adds xy, xz, yz tilt factors

   The cell matrix (rows = lattice vectors) is:
     | lx   0   0  |       lx = xhi - xlo
     | xy  ly   0  |       ly = yhi - ylo
     | xz  yz  lz |       lz = zhi - zlo

   From this we compute:
     a = |row1| = lx
     b = |row2| = √(xy² + ly²)
     c = |row3| = √(xz² + yz² + lz²)
     cos(α) = (row2 · row3) / (b × c)
     cos(β) = (row1 · row3) / (a × c)
     cos(γ) = (row1 · row2) / (a × b)

## 2. Cartesian → fractional coordinates
   CIF uses fractional coordinates (0 to 1). To convert:
     [fx, fy, fz] = inverse(cell_matrix) × [x, y, z]

   Then wrap into [0, 1):
     fx = mod(fx, 1.0)

## 3. Filter framework atoms
   RASPA3 only needs the framework (Si, O, Al), not adsorbate
   molecules. The `framework_types` argument selects which LAMMPS
   atom types to include (default: [1, 2] = Si and O).

## 4. Write CIF format
   The output follows the standard CIF 1.0 format with:
     - _cell_length_a/b/c and _cell_angle_alpha/beta/gamma
     - _atom_site_label, _atom_site_fract_x/y/z, _atom_site_type_symbol
   This is the format RASPA3 reads for framework definitions.

# Arguments
- `filename`: output CIF path
- `data`: LammpsData object
- `framework_types`: atom types to include (default [1, 2])
- `type_elements`: map atom type → element symbol
                   (default: {1=>"Si", 2=>"O", 3=>"Al"})
- `space_group`: space group to write (default "P 1" — no symmetry,
                 all atoms listed explicitly)
- `comment`: name for the structure
"""
function write_cif(filename::String, data::LammpsData;
                   framework_types::Vector{Int} = [1, 2],
                   type_elements::Dict{Int,String} = Dict(
                       1 => "Si", 2 => "O", 3 => "Al", 4 => "H"),
                   space_group::String = "P 1",
                   comment::String = "Zeolite framework from LAMMPS NPT")

    # ──────────────────────────────────────────────────────────────
    # Step 1: Build cell matrix and compute cell parameters
    # ──────────────────────────────────────────────────────────────
    lx = data.box_dimensions[1,2] - data.box_dimensions[1,1]
    ly = data.box_dimensions[2,2] - data.box_dimensions[2,1]
    lz = data.box_dimensions[3,2] - data.box_dimensions[3,1]
    xy, xz, yz = data.tilt_factors

    # Cell matrix: each ROW is a lattice vector
    #   a_vec = [lx,  0,  0]
    #   b_vec = [xy, ly,  0]
    #   c_vec = [xz, yz, lz]
    cell = [lx  0.0  0.0;
            xy  ly   0.0;
            xz  yz   lz]

    a_vec = cell[1, :]
    b_vec = cell[2, :]
    c_vec = cell[3, :]

    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)

    # Angles from dot products
    cos_alpha = dot(b_vec, c_vec) / (b * c)
    cos_beta  = dot(a_vec, c_vec) / (a * c)
    cos_gamma = dot(a_vec, b_vec) / (a * b)

    # Clamp to [-1, 1] to avoid NaN from floating-point noise
    clamp_cos(x) = clamp(x, -1.0, 1.0)
    alpha = rad2deg(acos(clamp_cos(cos_alpha)))
    beta  = rad2deg(acos(clamp_cos(cos_beta)))
    gamma = rad2deg(acos(clamp_cos(cos_gamma)))

    volume = abs(dot(a_vec, cross(b_vec, c_vec)))

    println("  Cell parameters:")
    @printf("    a = %.6f Å\n", a)
    @printf("    b = %.6f Å\n", b)
    @printf("    c = %.6f Å\n", c)
    @printf("    α = %.4f°\n", alpha)
    @printf("    β = %.4f°\n", beta)
    @printf("    γ = %.4f°\n", gamma)
    @printf("    V = %.2f ų\n", volume)

    # ──────────────────────────────────────────────────────────────
    # Step 2: Compute inverse cell matrix for Cartesian → fractional
    # ──────────────────────────────────────────────────────────────
    #   frac = inv(cell) * cart
    # where cell has lattice vectors as ROWS, so:
    #   frac = cart * inv(cell')  ... or equivalently:
    #   frac = inv(cell') * cart  when cell' has vectors as COLUMNS
    cell_inv = inv(cell')  # cell' = transpose, columns = lattice vectors

    # ──────────────────────────────────────────────────────────────
    # Step 3: Filter framework atoms
    # ──────────────────────────────────────────────────────────────
    natoms = size(data.coords, 1)
    fw_set = Set(framework_types)

    fw_indices = [j for j in 1:natoms if data.atom_labels[j] in fw_set]
    nfw = length(fw_indices)

    # Count per element
    elem_counts = Dict{String,Int}()
    for j in fw_indices
        elem = get(type_elements, data.atom_labels[j], "X")
        elem_counts[elem] = get(elem_counts, elem, 0) + 1
    end
    println("  Framework atoms: $nfw")
    for (elem, count) in sort(collect(elem_counts))
        println("    $elem: $count")
    end

    excluded = natoms - nfw
    excluded > 0 && println("  Excluded: $excluded adsorbate atoms")

    # ──────────────────────────────────────────────────────────────
    # Step 4: Convert to fractional and wrap to [0, 1)
    # ──────────────────────────────────────────────────────────────
    # Shift coordinates so origin is at box lo corner
    origin = [data.box_dimensions[1,1],
              data.box_dimensions[2,1],
              data.box_dimensions[3,1]]

    frac_coords = zeros(nfw, 3)
    for (i, j) in enumerate(fw_indices)
        cart = data.coords[j, :] .- origin
        frac = cell_inv * cart
        # Wrap to [0, 1)
        frac_coords[i, :] = mod.(frac, 1.0)
    end

    # ──────────────────────────────────────────────────────────────
    # Step 5: Write CIF file
    # ──────────────────────────────────────────────────────────────
    open(filename, "w") do io
        println(io, "#" ^ 70)
        println(io, "# CIF generated by write_cif.jl")
        println(io, "# $comment")
        println(io, "#" ^ 70)
        println(io)
        println(io, "data_zeolite_from_lammps")
        println(io)
        println(io, "_chemical_name_common  '$comment'")
        println(io)

        # Cell parameters
        @printf(io, "_cell_length_a    %.6f\n", a)
        @printf(io, "_cell_length_b    %.6f\n", b)
        @printf(io, "_cell_length_c    %.6f\n", c)
        @printf(io, "_cell_angle_alpha %.4f\n", alpha)
        @printf(io, "_cell_angle_beta  %.4f\n", beta)
        @printf(io, "_cell_angle_gamma %.4f\n", gamma)
        @printf(io, "_cell_volume      %.2f\n", volume)
        println(io)

        # Space group (P1 — all atoms listed explicitly)
        println(io, "_space_group_name_H-M_alt  '$space_group'")
        println(io, "_space_group_IT_number     1")
        println(io)

        # Symmetry operations (just identity for P1)
        println(io, "loop_")
        println(io, "_space_group_symop_operation_xyz")
        println(io, "   'x, y, z'")
        println(io)

        # Atom sites
        println(io, "loop_")
        println(io, "   _atom_site_label")
        println(io, "   _atom_site_occupancy")
        println(io, "   _atom_site_fract_x")
        println(io, "   _atom_site_fract_y")
        println(io, "   _atom_site_fract_z")
        println(io, "   _atom_site_adp_type")
        println(io, "   _atom_site_U_iso_or_equiv")
        println(io, "   _atom_site_type_symbol")

        # Track element counts for unique labels
        label_count = Dict{String,Int}()
        for (i, j) in enumerate(fw_indices)
            elem = get(type_elements, data.atom_labels[j], "X")
            label_count[elem] = get(label_count, elem, 0) + 1
            label = "$(elem)$(label_count[elem])"

            @printf(io, "   %-8s  1.0  %12.8f  %12.8f  %12.8f  Uiso  0.01  %s\n",
                    label, frac_coords[i,1], frac_coords[i,2], frac_coords[i,3], elem)
        end
    end

    println("  Wrote $filename")
    return nothing
end


"""
    write_framework_cif(lammps_data_file, cif_file; kwargs...)

Convenience function: read a LAMMPS data file and write the framework as CIF.

# Example
    write_framework_cif("loaded_npt_final.lmp", "relaxed_MFI.cif";
                        framework_types = [1, 2],
                        type_elements = Dict(1=>"Si", 2=>"O"))
"""
function write_framework_cif(lammps_file::String, cif_file::String;
                              framework_types::Vector{Int} = [1, 2],
                              type_elements::Dict{Int,String} = Dict(
                                  1=>"Si", 2=>"O", 3=>"Al", 4=>"H"),
                              space_group::String = "P 1",
                              comment::String = "")
    println("Reading $lammps_file...")
    data = LammpsDataReader.read_lammps_data(lammps_file; verbose=false)

    if isempty(comment)
        comment = "Framework from $lammps_file"
    end

    println("Writing CIF...")
    write_cif(cif_file, data;
              framework_types = framework_types,
              type_elements = type_elements,
              space_group = space_group,
              comment = comment)
end


# ══════════════════════════════════════════════════════════════════
# Standalone usage
# ══════════════════════════════════════════════════════════════════

function main_cif(args=ARGS)
    if length(args) < 2
        println("""
        Usage: julia write_cif.jl input.lmp output.cif [options]

        Reads a LAMMPS data file and writes a CIF of the framework
        atoms only (stripping adsorbate molecules).

        Options:
          --fw-types 1,2      Framework atom types (comma-separated)
          --elements 1=Si,2=O Map atom types to element symbols
          --comment "text"    Structure name in CIF header

        Example (all-silica):
          julia write_cif.jl loaded_npt_final.lmp relaxed_MFI.cif

        Example (aluminosilicate):
          julia write_cif.jl final.lmp MFI_Al.cif \\
                --fw-types 1,2,3 --elements 1=Si,2=O,3=Al
        """)
        return
    end

    input_file = args[1]
    output_file = args[2]

    fw_types = [1, 2]
    type_elem = Dict(1=>"Si", 2=>"O", 3=>"Al", 4=>"H")
    comment = ""

    i = 3
    while i <= length(args)
        a = args[i]
        if a == "--fw-types" && i < length(args)
            fw_types = parse.(Int, split(args[i+1], ","))
            i += 2
        elseif a == "--elements" && i < length(args)
            type_elem = Dict{Int,String}()
            for pair in split(args[i+1], ",")
                k, v = split(pair, "=")
                type_elem[parse(Int, k)] = String(v)
            end
            i += 2
        elseif a == "--comment" && i < length(args)
            comment = args[i+1]; i += 2
        else
            i += 1
        end
    end

    write_framework_cif(input_file, output_file;
                         framework_types = fw_types,
                         type_elements = type_elem,
                         comment = comment)
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "read_lammps_data.jl"))
    main_cif()
end
