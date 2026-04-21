"""
refine_topology_types.jl  (v2 — matches .ff numbering)

After add_zeolite_topology! builds connectivity with collapsed types
(all Si→6, all O→3), this refines bond/angle/dihedral labels using
the restored 8-type atom labels. Numbering matches hillsauer_alumsil.ff.

Atom types (corrected.data / OVITO):
  1=Al  2=Hb  3=O(Oss)  4=Oas  5=Ob  6=Si  7=Si_a  8=Si_b

═══════════════════════════════════════════════════════════════
  TYPE MAP  (authoritative — matches .ff and LAMMPS input)
═══════════════════════════════════════════════════════════════

  Bonds (6):
    1 = Si*─Oss      any Si(6,7,8) + O(3)       b₀=1.6104
    2 = Si_b─Ob      (8)─(5)                    b₀=1.6581
    3 = Si_a─Oas     (7)─(4)                    b₀=1.6157
    4 = Al─Ob        (1)─(5)                    b₀=1.9698
    5 = Al─Oas       (1)─(4)                    b₀=1.7193
    6 = Ob─Hb        (5)─(2)                    b₀=0.9540

  Angles (10):
    1  = Oss─Si*─Oss     O─T─O                  θ₀=112.02°
    2  = Si*─Oss─Si*     T─O─T (mHSFF)          θ₀=150.00°
    3  = Si_b─Ob─Al                              θ₀=136.66°
    4  = Si_b─Ob─Hb                              θ₀=119.28°
    5  = Al─Ob─Hb                                θ₀=107.12°
    6  = Al─Oas─Si_a                             θ₀=162.40°
    7  = Oas─Al─Ob                               θ₀= 81.50°
    8  = Oas─Al─Oas                              θ₀=113.40°
    9  = Oss─Si_b─Ob                             θ₀=107.61°
    10 = Oas─Si_a─Oss                            θ₀=112.43°

  Dihedrals (10):
    1  = Oss─Si*─Oss─Si*   bulk
    2  = Oas─Si_a─Oss─Si*
    3  = Oss─Si_a─Oas─Al
    4  = Ob─Si_b─Oss─Si*
    5  = Oss─Si_b─Ob─Al
    6  = Oss─Si_b─Ob─Hb
    7  = Oas─Al─Oas─Si_a
    8  = Oas─Al─Ob─Si_b
    9  = Oas─Al─Ob─Hb
    10 = Ob─Al─Oas─Si_a

  Impropers (1):
    1  = Si_b─Ob─Hb─Al
═══════════════════════════════════════════════════════════════
"""

# ── Atom type constants ──────────────────────────────────────
const _AL  = 1;  const _HB  = 2
const _OSS = 3;  const _OAS = 4;  const _OB = 5
const _SI_TYPES = Set([6, 7, 8])       # Si, Si_a, Si_b
const _T_TYPES  = Set([1, 6, 7, 8])    # Al + all Si

# ════════════════════════════════════════════════════════════════
#  BONDS
# ════════════════════════════════════════════════════════════════
function _ref_bond(t1::Int, t2::Int)
    a, b = minmax(t1, t2)
    a == _OSS && b in _SI_TYPES && return 1   # Si*-Oss
    a == _OB  && b == 8         && return 2   # Si_b-Ob
    a == _OAS && b == 7         && return 3   # Si_a-Oas
    a == _AL  && b == _OB       && return 4   # Al-Ob
    a == _AL  && b == _OAS      && return 5   # Al-Oas
    a == _HB  && b == _OB       && return 6   # Ob-Hb
    @warn "Unknown bond pair: ($t1, $t2)"; return 0
end

const _BOND_NAMES = Dict(
    1=>"Si*-Oss", 2=>"Si_b-Ob", 3=>"Si_a-Oas",
    4=>"Al-Ob",   5=>"Al-Oas",  6=>"Ob-Hb")

# ════════════════════════════════════════════════════════════════
#  ANGLES
# ════════════════════════════════════════════════════════════════
function _ref_angle(te1::Int, tc::Int, te2::Int)
    # ── Center is an O-type: T─O─T or T─O─H ──
    if tc == _OSS                                        # center = Oss
        te1 in _SI_TYPES && te2 in _SI_TYPES && return 2 # Si*-Oss-Si*  θ₀=150°
    elseif tc == _OAS                                    # center = Oas
        s1, s2 = minmax(te1, te2)
        s1 == _AL && s2 in _SI_TYPES && return 6         # Al-Oas-Si_a  θ₀=162.4°
        # rare: Si-Oas-Si, treat as type 6
        s1 in _SI_TYPES && s2 in _SI_TYPES && return 6
    elseif tc == _OB                                     # center = Ob
        s1, s2 = minmax(te1, te2)
        s1 == _AL && s2 in _SI_TYPES && return 3          # Si_b-Ob-Al   θ₀=136.7°
        # Si_b-Ob-Hb: one end is Si, other is Hb
        (s1 == _HB && s2 in _SI_TYPES) && return 4        # Si_b-Ob-Hb   θ₀=119.3°
        (s1 in _SI_TYPES && s2 == _HB) && return 4
        # Al-Ob-Hb
        (s1 == _AL && s2 == _HB) && return 5              # Al-Ob-Hb     θ₀=107.1°
        (s1 == _HB && s2 == _AL) && return 5
    end

    # ── Center is a T-atom: O─T─O ──
    if tc in _SI_TYPES
        s1, s2 = minmax(te1, te2)
        s1 == _OSS && s2 == _OSS && return 1              # Oss-Si*-Oss  θ₀=112°
        (s1 == _OSS && s2 == _OB)  && return 9             # Oss-Si_b-Ob  θ₀=107.6°
        (s1 == _OB  && s2 == _OSS) && return 9
        (s1 == _OAS && s2 == _OSS) && return 10            # Oas-Si_a-Oss θ₀=112.4°
        (s1 == _OSS && s2 == _OAS) && return 10
        # rare fallbacks
        s1 == _OAS && s2 == _OAS && return 10
        (s1 == _OAS && s2 == _OB) && return 9
    elseif tc == _AL
        s1, s2 = minmax(te1, te2)
        s1 == _OAS && s2 == _OAS && return 8               # Oas-Al-Oas   θ₀=113.4°
        (s1 == _OAS && s2 == _OB) && return 7              # Oas-Al-Ob    θ₀=81.5°
        (s1 == _OB && s2 == _OAS) && return 7
        s1 == _OB && s2 == _OB && return 7                  # Ob-Al-Ob (rare)
    end

    @warn "Unknown angle: ($te1, $tc, $te2)"; return 0
end

const _ANGLE_NAMES = Dict(
    1=>"Oss-Si*-Oss",   2=>"Si*-Oss-Si*",  3=>"Si_b-Ob-Al",
    4=>"Si_b-Ob-Hb",    5=>"Al-Ob-Hb",     6=>"Al-Oas-Si_a",
    7=>"Oas-Al-Ob",     8=>"Oas-Al-Oas",   9=>"Oss-Si_b-Ob",
    10=>"Oas-Si_a-Oss")

# ════════════════════════════════════════════════════════════════
#  DIHEDRALS
# ════════════════════════════════════════════════════════════════
function _ref_dihedral(t1::Int, t2::Int, t3::Int, t4::Int)
    # Pattern: O_a(t1) ─ T_inner(t2) ─ O_b(t3) ─ X_outer(t4)

    if t2 in _SI_TYPES   # inner T = Si*
        t1 == _OSS && t3 == _OSS && t4 in _T_TYPES  && return 1   # Oss-Si*-Oss-Si*
        t1 == _OAS && t3 == _OSS && t4 in _T_TYPES  && return 2   # Oas-Si_a-Oss-Si*
        t1 == _OSS && t3 == _OAS && t4 == _AL        && return 3   # Oss-Si_a-Oas-Al
        t1 == _OB  && t3 == _OSS && t4 in _T_TYPES  && return 4   # Ob-Si_b-Oss-Si*
        t1 == _OSS && t3 == _OB  && t4 == _AL        && return 5   # Oss-Si_b-Ob-Al
        t1 == _OSS && t3 == _OB  && t4 == _HB        && return 6   # Oss-Si_b-Ob-Hb
        # catch-all for Si
        t3 == _OSS && t4 in _T_TYPES && return 1
        t3 == _OAS && t4 == _AL      && return 3
        t3 == _OB  && t4 == _AL      && return 5
        t3 == _OB  && t4 == _HB      && return 6

    elseif t2 == _AL      # inner T = Al
        t1 == _OAS && t3 == _OAS && t4 in _SI_TYPES && return 7   # Oas-Al-Oas-Si_a
        t1 == _OAS && t3 == _OB  && t4 in _SI_TYPES && return 8   # Oas-Al-Ob-Si_b
        t1 == _OAS && t3 == _OB  && t4 == _HB       && return 9   # Oas-Al-Ob-Hb
        t1 == _OB  && t3 == _OAS && t4 in _SI_TYPES && return 10  # Ob-Al-Oas-Si_a
        # catch-all for Al
        t3 == _OAS && t4 in _T_TYPES && return 7
        t3 == _OB  && t4 in _SI_TYPES && return 8
        t3 == _OB  && t4 == _HB       && return 9
    end

    @warn "Unknown dihedral: ($t1, $t2, $t3, $t4)"; return 1
end

const _DIHEDRAL_NAMES = Dict(
    1=>"Oss-Si*-Oss-Si*",  2=>"Oas-Si_a-Oss-Si*",
    3=>"Oss-Si_a-Oas-Al",  4=>"Ob-Si_b-Oss-Si*",
    5=>"Oss-Si_b-Ob-Al",   6=>"Oss-Si_b-Ob-Hb",
    7=>"Oas-Al-Oas-Si_a",  8=>"Oas-Al-Ob-Si_b",
    9=>"Oas-Al-Ob-Hb",     10=>"Ob-Al-Oas-Si_a")

# ════════════════════════════════════════════════════════════════
#  IMPROPERS — bridging hydroxyl planarity
# ════════════════════════════════════════════════════════════════
function _add_impropers!(data)
    type_of = data.atom_labels
    id2idx = Dict(data.atom_ids[j] => j for j in eachindex(data.atom_ids))

    ob_indices = [i for i in eachindex(type_of) if type_of[i] == _OB]
    imps = Tuple{Int,Int,Int,Int}[]

    for ob in ob_indices
        ob_id = data.atom_ids[ob]
        si_b = 0; al = 0; hb = 0
        for k in 1:size(data.bonds, 1)
            a1, a2 = data.bonds[k,1], data.bonds[k,2]
            other_id = a1 == ob_id ? a2 : (a2 == ob_id ? a1 : 0)
            other_id == 0 && continue
            haskey(id2idx, other_id) || continue
            t = type_of[id2idx[other_id]]
            t in _SI_TYPES && (si_b = other_id)
            t == _AL       && (al = other_id)
            t == _HB       && (hb = other_id)
        end
        si_b > 0 && al > 0 && hb > 0 && push!(imps, (si_b, ob_id, hb, al))
    end

    n = length(imps)
    if n > 0
        data.impropers = zeros(Int, n, 4)
        data.improper_labels = ones(Int, n)
        for (k, (a,b,c,d)) in enumerate(imps)
            data.impropers[k,:] .= [a,b,c,d]
        end
        data.nimproper_types = 1
    end
    return n
end

# ════════════════════════════════════════════════════════════════
#  MAIN ENTRY POINT
# ════════════════════════════════════════════════════════════════
"""
    refine_topology_types!(data; verbose=true)

Refine bond/angle/dihedral labels to match the .ff numbering.
Must be called AFTER add_zeolite_topology!() AND after restoring
the original 8-type atom labels.

Type map matches hillsauer_alumsil_empty.ff exactly.
"""
function refine_topology_types!(data; verbose::Bool=true)
    type_of = data.atom_labels
    id2idx  = Dict(data.atom_ids[j] => j for j in eachindex(data.atom_ids))

    # ── Bonds ──
    for k in 1:size(data.bonds,1)
        t1 = type_of[id2idx[data.bonds[k,1]]]
        t2 = type_of[id2idx[data.bonds[k,2]]]
        data.bond_labels[k] = _ref_bond(t1, t2)
    end
    data.nbond_types = isempty(data.bond_labels) ? 0 : maximum(data.bond_labels)

    # ── Angles ──
    for k in 1:size(data.angles,1)
        t1 = type_of[id2idx[data.angles[k,1]]]
        tc = type_of[id2idx[data.angles[k,2]]]
        t2 = type_of[id2idx[data.angles[k,3]]]
        data.angle_labels[k] = _ref_angle(t1, tc, t2)
    end
    data.nangle_types = isempty(data.angle_labels) ? 0 : maximum(data.angle_labels)

    # ── Dihedrals ──
    for k in 1:size(data.dihedrals,1)
        t1 = type_of[id2idx[data.dihedrals[k,1]]]
        t2 = type_of[id2idx[data.dihedrals[k,2]]]
        t3 = type_of[id2idx[data.dihedrals[k,3]]]
        t4 = type_of[id2idx[data.dihedrals[k,4]]]
        data.dihedral_labels[k] = _ref_dihedral(t1, t2, t3, t4)
    end
    data.ndihedral_types = isempty(data.dihedral_labels) ? 0 : maximum(data.dihedral_labels)

    # ── Impropers ──
    ni = _add_impropers!(data)

    if verbose
        println("\n  ═══ Refined topology (v2, .ff-matched) ═══")
        for (name, labels, nd) in [
            ("Bonds",     data.bond_labels,     _BOND_NAMES),
            ("Angles",    data.angle_labels,     _ANGLE_NAMES),
            ("Dihedrals", data.dihedral_labels,  _DIHEDRAL_NAMES)]
            ntypes = maximum(labels; init=0)
            println("  $name: $(length(labels))  ($ntypes types)")
            for t in sort(unique(labels))
                println("    type $t  $(get(nd, t, "?")): $(count(==(t), labels))")
            end
        end
        println("  Impropers: $ni  ($(data.nimproper_types) type)")
    end
    return data
end
