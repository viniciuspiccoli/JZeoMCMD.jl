#!/usr/bin/env julia
# ==============================================================================
#  Löwenstein-Compliant Al Substitution in Zeolite Frameworks
#  with Spatial Zoning and T-site Siting Preference
#
#  References:
#    [1] Pelmenschikov et al. J. Phys. Chem. 1992, 96, 7051  (Löwenstein, H geometry)
#    [2] Le et al. Nat. Catal. 2023, 6, 254          (Si-zoning, diffusion)
#    [3] Shen, Qin et al. Angew. Chem. 2025, 64, e202508909 (Al siting, channels vs intersections)
#
# ==============================================================================
#  SCIENTIFIC BACKGROUND
# ==============================================================================
#
#  SPATIAL ZONING (Le et al. 2023)
#    ZSM-5 crystals can have mesoscale Al gradients: "Si-zoned" crystals have a
#    Si-rich exterior (high Si/Al) protecting an Al-rich interior; "Al-zoned"
#    crystals have the opposite (Al concentrated at the rim). Si-zoned crystals
#    show drastically reduced diffusion limitations and external coking,
#    extending catalyst lifetime up to 4× vs. homogeneous samples.
#
#  AL SITING: CHANNELS vs. INTERSECTIONS (Shen, Qin et al. 2025)
#    In MFI (ZSM-5), there are 12 crystallographically distinct T-sites:
#      · T-sites at channel intersections (T7, T8 in Pnma asymmetric unit):
#        These are preferred by TPA⁺ (tetrapropylammonium) as structure-
#        directing agent. Brønsted acid sites at intersections have a TOF
#        13× higher than those in channels (for ethylene aromatization at 973 K).
#      · T-sites in 10-MR channels (T1–T6, T9–T12):
#        Preferred by Na⁺ as SDA. These sites favor the alkene cycle in
#        MTO, giving propylene/higher olefins. More prone to coke formation.
#      · Mixed (TPA⁺ + Na⁺ synthesis): Al in both environments.
#
#  GEOMETRIC CORE Si/Al CALCULATION (hand derivation):
#    Let x = Al/(Si+Al) = 1/(Si/Al + 1)  (Al mole fraction)
#    Al balance:  x_total = f_outer · x_outer + f_inner · x_inner
#    Solving:     x_inner = (x_total − f_outer · x_outer) / f_inner
#
#    Example from image (f_outer = f_inner = 0.5):
#      x_total = 1/34 ≈ 0.02941   (Si/Al = 33)
#      x_outer = 1/19 ≈ 0.05263   (Si/Al = 18)
#      x_inner = 2·(1/34) − 1/19 = 1/17 − 1/19 = 2/323 ≈ 0.00619
#      Si/Al_inner = (1 − x_inner)/x_inner = 321/2 = 160.5  ← core is silicalite-like
#
#    Finite cell note: with only 3 Al atoms in 96 T-sites (Si/Al=33), the outer
#    zone (48 sites) needs ~2.5 Al and the inner zone ~0.3 Al, so in practice
#    2 Al go outer and 1 Al goes inner (Si/Al_inner ≈ 47) or all 3 go outer
#    (pure-Si core). The analytical value is approached for larger supercells.
#
# ==============================================================================
#  ALGORITHM WORKFLOW
# ==============================================================================
#
#  STEP 1  Parse CIF, extract cell parameters, atoms, charges.
#  STEP 2  Expand symmetry to full P1 unit cell; track which ASU (asymmetric unit)
#          site each expanded atom came from — this is what allows T-site
#          classification by label (Si7 → intersection, others → channel).
#  STEP 3  Identify Si T-sites and O atoms. Build Si–O–Si neighbor graph.
#  STEP 4  [OPTIONAL] Classify T-sites by SPATIAL ZONE:
#          Sort by Cartesian distance from cell centroid; inner fraction = inner
#          zone, remaining = outer zone. Zone targets:
#            N_Al_zone = round(N_T_zone / (Si/Al_zone + 1))
#          where Si/Al for the inner zone is derived from the Al balance formula.
#  STEP 5  [OPTIONAL] Classify T-sites by TYPE:
#          Intersection sites (asu_src ∈ {Si7, Si8}): TPA⁺-preferred sites.
#          Channel sites (all others): Na⁺-preferred sites.
#  STEP 6  Greedy Löwenstein-compliant Al placement with priority ordering:
#          Within each zone, the shuffled list is ordered by type preference
#          (intersection-first, channel-first, or homogeneous). Al is placed
#          zone by zone (inner first, then outer) with a shared blocked set,
#          guaranteeing global Löwenstein compliance.
#  STEP 7  Apply Si→Al substitution; adjust charges (Δq = −0.30 e).
#  STEP 8  Place Brønsted H on the LONGEST Al-O bond in each Al-O-Si bridge.
#          Physical basis: Al substitution stretches all Al-O bonds (Al is larger
#          than Si), but one bond elongates significantly more (~1.89 Å vs ~1.71 Å
#          for the other three). This longest bond points toward the pore channel
#          and is the true Brønsted acid site — confirmed by the reference
#          H-ZSM5-3-toluene-T7T8.cif, where all 3 Al sites follow this pattern.
#          H direction: external bisector of Al-O-Si angle (Pelmenschikov 1992),
#          R(O-H) = 0.97 Å.
#  STEP 9  Write P1 CIF.
#
# ==============================================================================
#  USAGE
# ==============================================================================
#
#   julia lowenstein_substitution.jl INPUT.cif SIAL [SEED] [OPTIONS]
#
#  OPTIONS:
#   --zone=OUTER_SIAL[:FRAC]   Si-zoned structure; outer zone Si/Al = OUTER_SIAL,
#                               outer volume fraction FRAC (default 0.5).
#                               Inner zone Si/Al is computed from Al balance.
#   --pref=channel              Prefer Al at channel T-sites (Na⁺-like synthesis)
#   --pref=intersection         Prefer Al at intersection T-sites (TPA⁺-like)
#   --pref=energy               Prefer Al at lowest-energy T-sites in MFI:
#                               T14 > T17 > T20 > T1, then others (random).
#                               Based on Grau-Crespo et al. PCCP 2000 and
#                               Ruiz-Salvador et al. JSC 2013.
#   --pref=dempsey              Maximise minimum Al-Al distance (Dempsey's rule).
#                               Among all Löwenstein-compliant placements, keep
#                               the one where the closest Al-Al pair is furthest.
#                               Used by Woodward et al. Catal. Commun. 2022.
#   --isites=Si7,Si8            Override intersection T-site labels (default Si7,Si8)
#   --esites=Si14,Si17,Si20,Si1 Override energy-ordered T-site labels (MFI default)
#   --supercell=NxMxP           Replicate the unit cell N×M×P times before
#                               substitution. E.g. --supercell=2x2x4 gives the
#                               ~5000-atom supercell used by Woodward et al. 2022.
#   --compute-core              Print zone analysis only; skip substitution
#
#  EXAMPLES:
#   julia lowenstein_substitution.jl MFI_SI.cif 33
#   julia lowenstein_substitution.jl MFI_SI.cif 19 42 --pref=intersection
#   julia lowenstein_substitution.jl MFI_SI.cif 47 42 --pref=energy
#   julia lowenstein_substitution.jl MFI_SI.cif 15 42 --pref=dempsey --supercell=2x2x4
#   julia lowenstein_substitution.jl MFI_SI.cif 33 42 --zone=18 --pref=channel
#   julia lowenstein_substitution.jl MFI_SI.cif 33 --zone=18 --compute-core
#
# ==============================================================================

using Printf, LinearAlgebra, Random

# ==============================================================================
#  DATA STRUCTURE
# ==============================================================================

struct CIFData
    a::Float64; b::Float64; c::Float64
    alpha::Float64; beta::Float64; gamma::Float64
    labels::Vector{String}; elements::Vector{String}
    frac_x::Vector{Float64}; frac_y::Vector{Float64}; frac_z::Vector{Float64}
    charges::Vector{Float64}
end

# ==============================================================================
#  STEP 1 — CIF PARSER
# ==============================================================================

function parse_cif(filename::String)
    lines = readlines(filename)
    a = b = c = alpha = beta = gamma = 0.0
    labels = String[]; elements = String[]
    fx = Float64[]; fy = Float64[]; fz = Float64[]; charges = Float64[]

    in_loop = false; loop_keys = String[]
    col_label = col_symbol = col_x = col_y = col_z = col_charge = -1

    for line in lines
        s = strip(line)
        isempty(s) && continue; startswith(s, "#") && continue

        if startswith(s, "_cell_length_a");     a     = parse(Float64, split(s)[2])
        elseif startswith(s, "_cell_length_b"); b     = parse(Float64, split(s)[2])
        elseif startswith(s, "_cell_length_c"); c     = parse(Float64, split(s)[2])
        elseif startswith(s, "_cell_angle_alpha"); alpha = parse(Float64, split(s)[2])
        elseif startswith(s, "_cell_angle_beta");  beta  = parse(Float64, split(s)[2])
        elseif startswith(s, "_cell_angle_gamma"); gamma = parse(Float64, split(s)[2])
        elseif s == "loop_"
            in_loop = true; loop_keys = String[]
            col_label = col_symbol = col_x = col_y = col_z = col_charge = -1
        elseif in_loop && startswith(s, "_")
            push!(loop_keys, s); idx = length(loop_keys)
            if occursin("_atom_site_label", s) && !occursin("number", s); col_label = idx
            elseif occursin("_atom_site_type_symbol", s); col_symbol = idx
            elseif s == "_atom_site_fract_x"; col_x = idx
            elseif s == "_atom_site_fract_y"; col_y = idx
            elseif s == "_atom_site_fract_z"; col_z = idx
            elseif occursin("_atom_site_charge", s); col_charge = idx
            end
        elseif in_loop && !startswith(s, "_") && !isempty(loop_keys)
            if col_x > 0 && col_y > 0 && col_z > 0
                p = split(s)
                length(p) >= maximum(filter(x->x>0,
                    [col_label,col_symbol,col_x,col_y,col_z])) || continue
                lbl = col_label  > 0 ? p[col_label]  : "X"
                sym = col_symbol > 0 ? p[col_symbol] : lbl
                push!(labels, lbl); push!(elements, sym)
                push!(fx, parse(Float64,p[col_x]))
                push!(fy, parse(Float64,p[col_y]))
                push!(fz, parse(Float64,p[col_z]))
                push!(charges, col_charge > 0 && length(p)>=col_charge ?
                               parse(Float64,p[col_charge]) : 0.0)
            elseif !startswith(s, "_")
                !any(c->c>0, [col_x,col_y,col_z]) && (in_loop = false)
            end
        else
            in_loop = false
        end
    end
    return CIFData(a,b,c,alpha,beta,gamma,labels,elements,fx,fy,fz,charges)
end

# ==============================================================================
#  LATTICE MATH
# ==============================================================================

function frac_to_cart_matrix(a,b,c,al_d,be_d,ga_d)
    al=deg2rad(al_d); be=deg2rad(be_d); ga=deg2rad(ga_d)
    cal=cos(al); cbe=cos(be); cga=cos(ga); sga=sin(ga)
    M=zeros(3,3)
    M[1,1]=a; M[1,2]=b*cga; M[1,3]=c*cbe
    M[2,2]=b*sga; M[2,3]=c*(cal-cbe*cga)/sga
    M[3,3]=c*sqrt(max(0.0,1-cal^2-cbe^2-cga^2+2*cal*cbe*cga))/sga
    return M
end

frac_to_cart(M,fx,fy,fz) = M * [fx,fy,fz]

function min_image_distance(M,f1,f2)
    df=f2.-f1; df=df.-round.(df); return norm(M*df)
end

function min_image_vector(M,f1,f2)
    df=f2.-f1; df=df.-round.(df); return M*df
end

mean_vec(x) = sum(x)/length(x)

# ==============================================================================
#  STEP 2 — SYMMETRY EXPANSION  (returns expanded CIF + asu_src tracking)
# ==============================================================================

function parse_symop_expr(expr::String)
    s=replace(strip(expr)," "=>""); rx=ry=rz=t=0.0
    for (v,idx) in [("x",1),("y",2),("z",3)]
        m=match(Regex("([+-]?[0-9]*(?:\\.[0-9]+)?)"*v),s)
        if m!==nothing
            cs=m.captures[1]
            val=(cs==""||cs=="+") ? 1.0 : cs == "-" ? - 1.0 : parse(Float64,cs)
            idx==1&&(rx=val); idx==2&&(ry=val); idx==3&&(rz=val)
        end
    end
    s2=replace(s,r"[+-]?[0-9]*(?:\.[0-9]+)?[xyz]"=>"")
    if !isempty(s2)
        mf=match(r"([+-]?\d+)/(\d+)",s2)
        if mf!==nothing; t=parse(Float64,mf.captures[1])/parse(Float64,mf.captures[2])
        else mn=match(r"([+-]?[\d.]+)",s2); mn!==nothing&&(t=parse(Float64,mn.captures[1])); end
    end
    return rx,ry,rz,t
end

function parse_symops(filename::String)
    lines=readlines(filename)
    symops=Tuple{Matrix{Float64},Vector{Float64}}[]
    in_s=false
    for line in lines
        s=strip(line); isempty(s)&&continue; startswith(s,"#")&&continue
        if s=="loop_"; in_s=false
        elseif s=="_symmetry_equiv_pos_as_xyz"; in_s=true
        elseif in_s
            if startswith(s,"_")||s=="loop_"; in_s=false
            else
                op=replace(s,r"['\"']"=>""); isempty(op)&&continue
                p=split(op,","); length(p)==3||continue
                R=zeros(3,3); t=zeros(3)
                for (i,pp) in enumerate(p)
                    rx,ry,rz,ti=parse_symop_expr(string(pp))
                    R[i,1]=rx; R[i,2]=ry; R[i,3]=rz; t[i]=ti
                end
                push!(symops,(R,t))
            end
        end
    end
    isempty(symops)&&push!(symops,(Float64[1 0 0;0 1 0;0 0 1],zeros(3)))
    return symops
end

"""
Expand the asymmetric unit to the full P1 cell.
Returns (expanded_cif, asu_src) where asu_src[i] is the original ASU label
for atom i (e.g. "Si7"), enabling T-site type classification.
"""
function expand_symmetry(cif::CIFData, symops; tol=0.01)
    new_labels=String[]; new_elements=String[]
    new_fx=Float64[]; new_fy=Float64[]; new_fz=Float64[]
    new_charges=Float64[]
    asu_src=String[]   # ← NEW: track origin ASU label

    for (R,tvec) in symops
        for i in 1:length(cif.labels)
            f=[cif.frac_x[i],cif.frac_y[i],cif.frac_z[i]]
            fn=R*f.+tvec; fn=fn.-floor.(fn)
            dup=false
            for j in 1:length(new_fx)
                df=[new_fx[j]-fn[1],new_fy[j]-fn[2],new_fz[j]-fn[3]]
                df.-=round.(df)
                if all(abs.(df).<tol); dup=true; break; end
            end
            dup&&continue
            push!(new_elements,cif.elements[i]); push!(new_charges,cif.charges[i])
            push!(new_fx,fn[1]); push!(new_fy,fn[2]); push!(new_fz,fn[3])
            push!(new_labels,"TMP")
            push!(asu_src, cif.labels[i])   # ← store original ASU label
        end
    end

    counts=Dict{String,Int}()
    for i in 1:length(new_labels)
        raw=replace(strip(new_elements[i]),r"[^A-Za-z]+"=>"")
        elem=uppercasefirst(lowercase(raw))
        counts[elem]=get(counts,elem,0)+1
        new_labels[i]="$(elem)$(counts[elem])"
    end
    println("  Symmetry expansion: $(length(cif.labels)) → $(length(new_labels)) atoms")
    exp=CIFData(cif.a,cif.b,cif.c,cif.alpha,cif.beta,cif.gamma,
                new_labels,new_elements,new_fx,new_fy,new_fz,new_charges)
    return exp, asu_src
end

# ==============================================================================
#  SUPERCELL GENERATION  (Woodward et al. Catal. Commun. 2022 use 2×2×4)
# ==============================================================================

"""
Replicate the P1 unit cell na×nb×nc times.
Fractional coordinates are tiled and wrapped; labels are renumbered.
Cell lengths scale accordingly; angles are unchanged.
"""
function make_supercell(cif::CIFData, na::Int, nb::Int, nc::Int)
    (na == 1 && nb == 1 && nc == 1) && return cif

    new_labels   = String[]
    new_elements = String[]
    new_fx = Float64[]; new_fy = Float64[]; new_fz = Float64[]
    new_charges  = Float64[]

    N = length(cif.labels)
    for ia in 0:na-1, ib in 0:nb-1, ic in 0:nc-1
        for i in 1:N
            push!(new_elements, cif.elements[i])
            push!(new_charges,  cif.charges[i])
            push!(new_fx, (cif.frac_x[i] + ia) / na)
            push!(new_fy, (cif.frac_y[i] + ib) / nb)
            push!(new_fz, (cif.frac_z[i] + ic) / nc)
            push!(new_labels, "TMP")
        end
    end

    # Renumber labels sequentially per element
    counts = Dict{String,Int}()
    for i in 1:length(new_labels)
        raw  = replace(strip(new_elements[i]), r"[^A-Za-z]+" => "")
        elem = uppercasefirst(lowercase(raw))
        counts[elem] = get(counts, elem, 0) + 1
        new_labels[i] = "$(elem)$(counts[elem])"
    end

    n_atoms = length(new_labels)
    println("  Supercell $(na)×$(nb)×$(nc): $(N) → $(n_atoms) atoms")

    return CIFData(cif.a * na, cif.b * nb, cif.c * nc,
                   cif.alpha, cif.beta, cif.gamma,
                   new_labels, new_elements,
                   new_fx, new_fy, new_fz, new_charges)
end

# ==============================================================================
#  STEP 3 — T-SITE CONNECTIVITY
# ==============================================================================

function build_tsite_connectivity(cif::CIFData; t_o_cutoff=2.1)
    M=frac_to_cart_matrix(cif.a,cif.b,cif.c,cif.alpha,cif.beta,cif.gamma)
    N=length(cif.labels)
    t_indices=Int[]; o_indices=Int[]
    for i in 1:N
        e=replace(uppercase(strip(cif.elements[i])),r"[^A-Z]+"=>"")
        e=="SI"&&push!(t_indices,i); e=="O"&&push!(o_indices,i)
    end
    println("  Found $(length(t_indices)) Si (T-sites) and $(length(o_indices)) O atoms")

    neighbors=Dict{Int,Set{Int}}()
    for ti in t_indices; neighbors[ti]=Set{Int}(); end
    o_bridges=Dict{Int,Vector{Int}}()

    bridging=0
    for oi in o_indices
        fo=[cif.frac_x[oi],cif.frac_y[oi],cif.frac_z[oi]]
        bonded=Int[]
        for ti in t_indices
            ft=[cif.frac_x[ti],cif.frac_y[ti],cif.frac_z[ti]]
            min_image_distance(M,fo,ft)<t_o_cutoff&&push!(bonded,ti)
        end
        if length(bonded)==2
            push!(neighbors[bonded[1]],bonded[2])
            push!(neighbors[bonded[2]],bonded[1])
            o_bridges[oi]=bonded; bridging+=1
        end
    end
    println("  Found $bridging bridging oxygens")
    nn=[length(neighbors[ti]) for ti in t_indices]
    println("  T-site coordination: min=$(minimum(nn)), max=$(maximum(nn)), " *
            "mean=$(@sprintf("%.1f",mean_vec(nn)))")
    return t_indices,o_indices,neighbors,o_bridges
end

# ==============================================================================
#  GEOMETRIC ZONE ANALYSIS
# ==============================================================================

"""
Compute the Si/Al ratio of the inner (core) zone from the Al balance:

  x = Al/(Si+Al) = 1/(Si/Al + 1)

  Al balance: x_total = f_outer·x_outer + f_inner·x_inner
  Solving:    x_inner = (x_total − f_outer·x_outer) / f_inner
  Si/Al_inner = (1 − x_inner) / x_inner

Example (image): Si/Al_total=33, Si/Al_outer=18, f_outer=0.5
  x_inner = 2/34 − 1/19 = 1/17 − 1/19 = 2/323 ≈ 0.00619
  Si/Al_inner = 321/2 = 160.5  (essentially silicalite-like core)
"""
function compute_core_sial(sial_total::Float64, sial_outer::Float64;
                            f_outer::Float64=0.5)
    f_inner = 1.0 - f_outer
    f_inner < 1e-10 && return Inf
    x_t = 1.0/(sial_total+1.0)
    x_o = 1.0/(sial_outer+1.0)
    x_i = (x_t - f_outer*x_o)/f_inner
    x_i <= 1e-10 && return Inf
    return (1.0 - x_i)/x_i
end

function print_zone_analysis(sial_total, sial_outer, n_t; f_outer=0.5)
    f_inner   = 1.0 - f_outer
    sial_inner = compute_core_sial(sial_total, sial_outer; f_outer=f_outer)
    n_al_tot  = n_t/(sial_total+1)
    n_al_out  = f_outer*n_t/(sial_outer+1)
    n_al_in   = max(0.0, n_al_tot - n_al_out)

    println("\n  ┌─ Geometric zone analysis ────────────────────────┐")
    println("  │ Al balance: x_inner = (x_total − f_outer·x_outer)/f_inner")
    println("  │")
    println("  │ Overall  Si/Al = $(@sprintf("%6.1f",sial_total))" *
            "   →  N_Al ≈ $(@sprintf("%.1f",n_al_tot)) / $n_t T-sites")
    println("  │ Outer ($(@sprintf("%.0f",f_outer*100))%)  Si/Al = $(@sprintf("%6.1f",sial_outer))" *
            "   →  N_Al ≈ $(@sprintf("%.1f",n_al_out))")
    if isinf(sial_inner)
        println("  │ Inner ($(@sprintf("%.0f",f_inner*100))%)  Si/Al = ∞  (pure-Si core)")
    else
        println("  │ Inner ($(@sprintf("%.0f",f_inner*100))%)  Si/Al = $(@sprintf("%6.1f",sial_inner))" *
                "   →  N_Al ≈ $(@sprintf("%.1f",n_al_in))")
    end
    println("  └──────────────────────────────────────────────────┘\n")
    return sial_inner
end

# ==============================================================================
#  T-SITE SPATIAL ZONE AND TYPE CLASSIFIERS
# ==============================================================================

"""
Classify T-sites as inner (closer to cell centroid) or outer by distance.
inner_fraction controls the boundary (default 0.5 = inner 50% by count).

Rationale: the centroid of the T-site cloud approximates the crystal core.
Sorting by Cartesian distance and taking the inner_fraction closest sites as
the "core" is shape-agnostic and works for any crystal morphology.
"""
function classify_tsite_zones(t_indices, cif::CIFData, M; inner_fraction=0.5)
    isempty(t_indices) && return Set{Int}(), Set{Int}()
    n_inner = max(1, round(Int, inner_fraction * length(t_indices)))
    cx = mean_vec([cif.frac_x[ti] for ti in t_indices])
    cy = mean_vec([cif.frac_y[ti] for ti in t_indices])
    cz = mean_vec([cif.frac_z[ti] for ti in t_indices])
    centroid = [cx,cy,cz]
    dists = [norm(min_image_vector(M, centroid,
             [cif.frac_x[ti],cif.frac_y[ti],cif.frac_z[ti]])) for ti in t_indices]
    perm = sortperm(dists)
    inner = Set(t_indices[perm[1:n_inner]])
    outer = Set(t_indices[perm[n_inner+1:end]])
    println("  Spatial zones: $(length(inner)) inner T-sites, " *
            "$(length(outer)) outer T-sites")
    return inner, outer
end

"""
Classify T-sites as channel-type or intersection-type using the original
ASU label tracked during symmetry expansion.

MFI channel intersections (T7, T8 in Pnma asymmetric unit):
  · Preferred by TPA⁺ SDA → favors aromatic cycle, higher TOF per site
    (Shen, Qin et al. Angew. Chem. 2025; Le et al. Nat. Catal. 2023)
  · Al at T8 (straight channel near intersection) used in MD simulations
    of diffusion by Ghorbanpour et al. 2014

MFI channels (T1–T6, T9–T12):
  · Preferred by Na⁺ SDA → favors alkene cycle, propylene selectivity
  · More susceptible to coke formation (Shen, Qin et al. 2025)
"""
function classify_tsite_types(t_indices, asu_src::Vector{String};
                               intersection_source::Vector{String}=["Si7","Si8"])
    channel_sites = Int[]; intersect_sites = Int[]
    for ti in t_indices
        src = ti <= length(asu_src) ? asu_src[ti] : ""
        src in intersection_source ? push!(intersect_sites,ti) : push!(channel_sites,ti)
    end
    n_total = length(t_indices)
    println("  T-site types: $(length(channel_sites)) channel " *
            "($(@sprintf("%.0f",100*length(channel_sites)/n_total))%), " *
            "$(length(intersect_sites)) intersection " *
            "($(@sprintf("%.0f",100*length(intersect_sites)/n_total))%)")
    println("  Intersection source labels: " * join(intersection_source, ", "))
    if isempty(intersect_sites) && intersection_source != [""]
        println("  WARNING: No intersection T-sites found. Check --isites labels " *
                "against ASU labels in your CIF (e.g. 'Si7', 'Si8').")
    end
    return channel_sites, intersect_sites
end

"""
Build a priority-ordered shuffled list of T-sites from a zone,
applying type preference within that zone.
"""
function priority_list(zone_sites::Set{Int}, channel_sites, intersect_sites,
                       pref::Symbol)
    ch = [t for t in channel_sites  if t in zone_sites]
    is = [t for t in intersect_sites if t in zone_sites]
    if pref == :intersection
        return [shuffle(is)..., shuffle(ch)...]
    elseif pref == :channel
        return [shuffle(ch)..., shuffle(is)...]
    else  # :homogeneous
        return shuffle(collect(zone_sites))
    end
end

# ==============================================================================
#  STEPS 5–7 — LÖWENSTEIN SUBSTITUTION WITH ZONE + TYPE PREFERENCE
# ==============================================================================

"""
Greedy Löwenstein-compliant Al substitution supporting:
  · Homogeneous: single Si/Al target across all T-sites
  · Zoned: separate Si/Al targets for inner and outer zones
  · T-site preference: channel-first, intersection-first, or random (homogeneous)

Löwenstein compliance is guaranteed by construction (blocking all T-site
neighbors of any newly placed Al). A final violation count is printed as a
sanity check and should always be 0.
"""
function lowenstein_substitute!(cif::CIFData, target_si_al::Float64;
                                 # Zoning
                                 zoned::Bool=false,
                                 outer_si_al::Float64=Inf,
                                 inner_fraction::Float64=0.5,
                                 # T-site preference
                                 al_preference::Symbol=:homogeneous,
                                 asu_src::Vector{String}=String[],
                                 intersection_source::Vector{String}=["Si7","Si8"],
                                 energy_source::Vector{String}=["Si14","Si17","Si20","Si1"],
                                 # Algorithm
                                 max_attempts::Int=5000, seed::Int=42)
    Random.seed!(seed)
    t_indices,o_indices,neighbors,o_bridges = build_tsite_connectivity(cif)
    N_T = length(t_indices)
    M   = frac_to_cart_matrix(cif.a,cif.b,cif.c,cif.alpha,cif.beta,cif.gamma)

    # ── Zone setup ─────────────────────────────────────────────────────────────
    local inner_sites::Set{Int}, outer_sites::Set{Int}
    local n_al_inner::Int, n_al_outer::Int

    if zoned && !isinf(outer_si_al)
        inner_sites, outer_sites = classify_tsite_zones(t_indices, cif, M;
                                                         inner_fraction=inner_fraction)
        f_outer = 1.0 - inner_fraction
        sial_inner = compute_core_sial(target_si_al, outer_si_al; f_outer=f_outer)
        n_al_outer = max(0, round(Int, length(outer_sites)/(outer_si_al+1)))
        n_al_inner = if isinf(sial_inner)
            0
        else
            max(0, round(Int, length(inner_sites)/(sial_inner+1)))
        end
        println("\n  Zone targets:")
        println("    Inner ($(length(inner_sites)) T-sites, Si/Al=$(@sprintf("%.1f",sial_inner))): N_Al = $n_al_inner")
        println("    Outer ($(length(outer_sites)) T-sites, Si/Al=$(@sprintf("%.1f",outer_si_al))): N_Al = $n_al_outer")
        println("    Total N_Al = $(n_al_inner + n_al_outer)  (overall Si/Al ≈ $(@sprintf("%.1f",target_si_al)))")
    else
        inner_sites = Set(t_indices)
        outer_sites = Set{Int}()
        n_al_inner  = round(Int, N_T/(target_si_al+1))
        n_al_outer  = 0
        println("\n  Homogeneous: N_T=$N_T, N_Al target=$n_al_inner, " *
                "Si/Al target=$(@sprintf("%.1f",target_si_al))")
    end

    # ── T-site type classification ─────────────────────────────────────────────
    channel_sites, intersect_sites = if !isempty(asu_src) &&
                                        al_preference in (:channel, :intersection)
        classify_tsite_types(t_indices, asu_src; intersection_source=intersection_source)
    else
        t_indices, Int[]
    end

    # Energy-ordered sites for MFI (Grau-Crespo 2000; Ruiz-Salvador 2013):
    # T14 > T17 > T20 > T1 are the lowest-energy Al positions.
    # Sites in energy_source are listed in priority order; remaining sites follow randomly.
    energy_priority = if al_preference == :energy && !isempty(asu_src)
        prio  = [ti for es in energy_source
                    for ti in t_indices if ti <= length(asu_src) && asu_src[ti] == es]
        rest  = [ti for ti in t_indices if !(ti in Set(prio))]
        (prio, rest)
    else
        (Int[], t_indices)
    end

    println("  Al preference: $(al_preference)")
    al_preference == :energy &&
        println("  Energy priority labels: " * join(energy_source, " > "))
    al_preference == :dempsey &&
        println("  Dempsey mode: maximise minimum Al-Al distance over $max_attempts shuffles")
    n_al_total = n_al_inner + n_al_outer
    n_al_total == 0 && error("N_Al target = 0: Si/Al ratio too high for this cell size.")
    println("")

    # ── Greedy placement ───────────────────────────────────────────────────────
    best_al_set   = Set{Int}()
    best_n_inner  = 0; best_n_outer = 0
    best_min_al_al = 0.0   # for Dempsey: track best minimum Al-Al distance

    for attempt in 1:max_attempts
        al_set  = Set{Int}()
        blocked = Set{Int}()

        function place_zone!(zone_sites, n_target)
            n_placed = 0
            ordered = if al_preference == :energy
                prio_in_zone = [t for t in energy_priority[1] if t in zone_sites]
                rest_in_zone = shuffle([t for t in energy_priority[2] if t in zone_sites])
                [prio_in_zone..., rest_in_zone...]
            else
                priority_list(zone_sites, channel_sites, intersect_sites, al_preference)
            end
            for ti in ordered
                n_placed >= n_target && break
                ti in blocked && continue
                push!(al_set, ti)
                for nb in neighbors[ti]; push!(blocked, nb); end
                n_placed += 1
            end
            return n_placed
        end

        ni = place_zone!(inner_sites, n_al_inner)
        no = place_zone!(outer_sites, n_al_outer)
        n_placed = ni + no

        if al_preference == :dempsey && n_placed == n_al_total && n_placed > 1
            # Compute minimum Al-Al Cartesian distance (PBC-aware)
            al_list = collect(al_set)
            min_d = Inf
            for i in 1:length(al_list)-1, j in i+1:length(al_list)
                fi = [cif.frac_x[al_list[i]], cif.frac_y[al_list[i]], cif.frac_z[al_list[i]]]
                fj = [cif.frac_x[al_list[j]], cif.frac_y[al_list[j]], cif.frac_z[al_list[j]]]
                d  = min_image_distance(M, fi, fj)
                d < min_d && (min_d = d)
            end
            if min_d > best_min_al_al
                best_min_al_al = min_d
                best_al_set = copy(al_set)
                best_n_inner = ni; best_n_outer = no
            end
        else
            if n_placed > length(best_al_set)
                best_al_set = copy(al_set)
                best_n_inner = ni; best_n_outer = no
            end
            (ni >= n_al_inner && no >= n_al_outer &&
             al_preference != :dempsey) && break
        end
    end

    al_preference == :dempsey && best_min_al_al > 0 &&
        println("  Dempsey best min Al-Al distance: $(@sprintf("%.2f",best_min_al_al)) Å")

    # ── Results ────────────────────────────────────────────────────────────────
    N_Al = length(best_al_set)
    actual_ratio = (N_T - N_Al) / N_Al

    println("  Placed $N_Al Al atoms total:")
    zoned && println("    Inner zone: $best_n_inner Al" *
                     (best_n_inner < n_al_inner ? " (WARNING: target was $n_al_inner)" : ""))
    zoned && println("    Outer zone: $best_n_outer Al" *
                     (best_n_outer < n_al_outer ? " (WARNING: target was $n_al_outer)" : ""))
    println("  Actual Si/Al = $(@sprintf("%.1f", actual_ratio))")

    # Sanity check: count Löwenstein violations (must be 0)
    violations = 0
    for ai in best_al_set
        for nb in neighbors[ai]; nb in best_al_set && (violations += 1); end
    end
    violations ÷= 2
    println("  Löwenstein violations: $violations")
    violations > 0 && println("  WARNING: violations detected — check neighbor graph!")

    # Report T-site type distribution of placed Al
    if !isempty(intersect_sites)
        n_al_is = count(t -> t in Set(intersect_sites), best_al_set)
        n_al_ch = N_Al - n_al_is
        println("  Al type distribution: $n_al_ch channel, $n_al_is intersection")
    end

    # Apply substitution (Step 7)
    for ti in best_al_set
        cif.elements[ti] = "Al"
        old = cif.labels[ti]
        new = replace(old, r"Si" => "Al")
        new == old && (new = "Al" * replace(old, r"[A-Za-z]+" => ""))
        cif.labels[ti] = new
        any(q -> q != 0.0, cif.charges) && (cif.charges[ti] -= 0.30)
    end

    return best_al_set, actual_ratio, o_bridges
end

# ==============================================================================
#  STEP 8 — BRØNSTED ACID H PLACEMENT  (Pelmenschikov et al. 1992)
# ==============================================================================

"""
Place one H per Al on its Al–O–Si bridge oxygen, using the external bisector
of the Al–O–Si angle (Pelmenschikov et al. J. Phys. Chem. 1992, 96, 7051):

  direction_H = −(û_Al + û_Si) / |û_Al + û_Si|

where û_Al and û_Si are PBC-correct unit vectors O→Al and O→Si.
R(O–H) = 0.97 Å. Charge: +0.40 e.
"""
function add_bronsted_H!(cif::CIFData, al_set::Set{Int}, o_bridges::Dict;
                         oh_distance=0.97)
    M     = frac_to_cart_matrix(cif.a,cif.b,cif.c,cif.alpha,cif.beta,cif.gamma)
    M_inv = inv(M)
    h_count   = 0
    skip_count = 0
    bronsted_o_map  = Dict{Int,Int}()   # al_idx → Brønsted O index (Ob)
    bronsted_si_map = Dict{Int,Int}()   # al_idx → Si partner of Ob (Si_b)

    # Build Al → [(o_idx, si_idx, d_al_o)] directly from o_bridges.
    # This eliminates any cutoff mismatch: we use the SAME bridge data that was
    # used to build the neighbor graph, so every Al is guaranteed to find its O.
    #
    # The Brønsted O is the Al–O bond with the LARGEST d(Al–O) among Al–O–Si bridges.
    # Physical basis: the elongated Al–O bond (~1.89 Å vs ~1.71 Å for the other three)
    # points into the pore channel and is the true Brønsted acid site.
    # Verified against H-ZSM5-3-toluene-T7T8.cif: deviation from external bisector < 8°.

    for al_idx in al_set
        fal = [cif.frac_x[al_idx], cif.frac_y[al_idx], cif.frac_z[al_idx]]

        # Collect all (o_idx, si_idx, d_al_o) from o_bridges where al_idx is one endpoint
        candidates = Tuple{Int,Int,Float64}[]   # (o_idx, si_idx, d_al_o)
        for (oi, t_pair) in o_bridges
            al_idx in t_pair || continue         # this O must be bonded to this Al
            # Find the Si partner (the other T-site in the bridge)
            for tj in t_pair
                tj == al_idx && continue
                te = replace(uppercase(strip(cif.elements[tj])), r"[^A-Z]+" => "")
                te == "SI" || continue            # skip if partner became Al (Löwenstein)
                fo  = [cif.frac_x[oi], cif.frac_y[oi], cif.frac_z[oi]]
                d   = min_image_distance(M, fal, fo)
                push!(candidates, (oi, tj, d))
            end
        end

        if isempty(candidates)
            println("  WARNING: Al at index $al_idx has no Al–O–Si bridge — skipping H placement.")
            skip_count += 1
            continue
        end

        # Pick the candidate with the LONGEST Al–O distance → pore-facing bond
        best_oi, best_si, _ = candidates[argmax(last.(candidates))]
        bronsted_o_map[al_idx]  = best_oi   # This O becomes Ob
        bronsted_si_map[al_idx] = best_si   # This Si becomes Si_b

        fo  = [cif.frac_x[best_oi], cif.frac_y[best_oi], cif.frac_z[best_oi]]
        fsi = [cif.frac_x[best_si], cif.frac_y[best_si], cif.frac_z[best_si]]

        # External bisector of Al–O–Si angle (Pelmenschikov et al. 1992)
        vec_o_al = min_image_vector(M, fo, fal)
        vec_o_si = min_image_vector(M, fo, fsi)
        u_al = vec_o_al / norm(vec_o_al)
        u_si = vec_o_si / norm(vec_o_si)
        internal  = u_al .+ u_si
        n_int     = norm(internal)
        direction = n_int > 1e-6 ? -(internal ./ n_int) : -u_al

        cart_o  = frac_to_cart(M, fo...)
        cart_h  = cart_o .+ oh_distance .* direction
        frac_h  = M_inv * cart_h
        frac_h  = frac_h .- floor.(frac_h)

        h_count += 1
        push!(cif.labels,   "H_BAS$h_count")
        push!(cif.elements, "Hb")  # Must match force_field.json pseudo-atom name
        push!(cif.frac_x,   frac_h[1])
        push!(cif.frac_y,   frac_h[2])
        push!(cif.frac_z,   frac_h[3])
        push!(cif.charges,  0.40)
    end

    println("  Added $h_count Brønsted H atoms (R(O–H) = $(oh_distance) Å)")
    skip_count > 0 &&
        println("  WARNING: $skip_count Al site(s) received no H — check framework connectivity.")
    h_count != length(al_set) &&
        println("  WARNING: H count ($h_count) ≠ Al count ($(length(al_set))) — mismatch!")
    return h_count, bronsted_o_map, bronsted_si_map
end

# ==============================================================================
#  STEP 8b — HILL-SAUER TYPE CLASSIFICATION
# ==============================================================================

"""
Classify all framework atoms into proper Hill-Sauer pseudo-atom types and
assign self-consistent charges from the H-S 1995 bond increment model.

Bond increments δ (Hill & Sauer 1995, Table 8):
  Si–Oss: 0.1309   Si–Oas: 0.1265   Si–Ob: 0.1392
  Al–Oas: 0.1694   Al–Ob:  0.0284   Hb–Ob: 0.0839

Resulting atom types and charges:
  Si      4 Oss bonds (far from Al)          q = +0.5236
  Si_a    3 Oss + 1 Oas (next to Al via Oas) q = +0.5192
  Si_b    3 Oss + 1 Ob (next to Brønsted O)  q = +0.5319
  Al      3 Oas + 1 Ob                       q = +0.5366
  Oss     Si–O–Si bridge                     q = −0.2618
  Oas     Si–O–Al bridge (no H)              q = −0.2959
  Ob      Si–O(H)–Al (Brønsted bridge)       q = −0.2515
  Hb      Brønsted proton                    q = +0.0839

Charge neutrality is guaranteed by the bond increment model.
"""
function classify_framework_types!(cif::CIFData, al_set::Set{Int},
                                    o_bridges::Dict, bronsted_o_map::Dict,
                                    bronsted_si_map::Dict)
    # Hill-Sauer 1995 charges (from bond increments, Table 8)
    q_Si   =  0.5236   # 4 × δ(Si-Oss) = 4 × 0.1309
    q_Si_a =  0.5192   # 3 × δ(Si-Oss) + δ(Si-Oas) = 3(0.1309) + 0.1265
    q_Si_b =  0.5319   # 3 × δ(Si-Oss) + δ(Si-Ob)  = 3(0.1309) + 0.1392
    q_Al   =  0.5366   # 3 × δ(Al-Oas) + δ(Al-Ob)  = 3(0.1694) + 0.0284
    q_O    = -0.2618   # Oss: -2 × δ(Si-Oss) = -2 × 0.1309
    q_Oas  = -0.2959   # -(δ(Si-Oas) + δ(Al-Oas)) = -(0.1265 + 0.1694)
    q_Ob   = -0.2515   # -(δ(Si-Ob) + δ(Al-Ob) + δ(Hb-Ob)) = -(0.1392 + 0.0284 + 0.0839)
    q_Hb   =  0.0839   # δ(Hb-Ob) = 0.0839

    ob_set  = Set(values(bronsted_o_map))   # O indices that are Ob
    sib_set = Set(values(bronsted_si_map))  # Si indices that are Si_b

    # Collect all O bonded to Al (but not Ob) → these are Oas
    oas_set = Set{Int}()
    for (oi, t_pair) in o_bridges
        oi in ob_set && continue             # skip Ob, already classified
        any(ti -> ti in al_set, t_pair) || continue  # must be bonded to Al
        push!(oas_set, oi)
    end

    # Collect all Si bonded to Al through Oas → these are Si_a
    # (Si on the other side of each Oas bridge from Al)
    sia_set = Set{Int}()
    for (oi, t_pair) in o_bridges
        oi in oas_set || continue            # must be an Oas oxygen
        for ti in t_pair
            ti in al_set && continue         # skip the Al itself
            # This T-site is Si bonded to Al via Oas → Si_a
            ti in sib_set && continue        # skip if already Si_b
            push!(sia_set, ti)
        end
    end

    # Now classify and assign charges
    n_classified = Dict("Al"=>0, "Si_a"=>0, "Si_b"=>0, "Ob"=>0, "Oas"=>0,
                        "Hb"=>0, "Si"=>0, "O"=>0)

    for i in 1:length(cif.labels)
        elem_upper = uppercase(replace(strip(cif.elements[i]), r"[^A-Za-z]+" => ""))

        if i in al_set
            cif.elements[i] = "Al";   cif.charges[i] = q_Al;   n_classified["Al"] += 1
        elseif i in sib_set
            cif.elements[i] = "Si_b"; cif.charges[i] = q_Si_b; n_classified["Si_b"] += 1
        elseif i in sia_set
            cif.elements[i] = "Si_a"; cif.charges[i] = q_Si_a; n_classified["Si_a"] += 1
        elseif i in ob_set
            cif.elements[i] = "Ob";   cif.charges[i] = q_Ob;   n_classified["Ob"] += 1
        elseif i in oas_set
            cif.elements[i] = "Oas";  cif.charges[i] = q_Oas;  n_classified["Oas"] += 1
        elseif elem_upper == "HB"
            cif.elements[i] = "Hb";   cif.charges[i] = q_Hb;   n_classified["Hb"] += 1
        elseif elem_upper == "SI"
            cif.elements[i] = "Si";   cif.charges[i] = q_Si;   n_classified["Si"] += 1
        elseif elem_upper == "O"
            cif.elements[i] = "O";    cif.charges[i] = q_O;    n_classified["O"] += 1
        end
    end

    println("  Hill-Sauer 1995 type classification:")
    for (t,n) in sort(collect(n_classified))
        n > 0 && println("    $t: $n")
    end

    # Verify charge neutrality
    total_q = sum(cif.charges)
    println("  Net charge: $(@sprintf("%.8f", total_q)) e")
    abs(total_q) > 0.01 && println("  ⚠️  WARNING: not neutral — check classification!")

    # Warn if any Si_b is bonded to >1 Ob (unusual, needs different charge)
    for si_idx in sib_set
        n_ob = count(oi -> oi in ob_set, [oi for (oi,tp) in o_bridges if si_idx in tp])
        n_ob > 1 && println("  ⚠️  Si_b at index $si_idx bonded to $n_ob Ob — " *
                             "charge may need manual adjustment")
    end

    return n_classified
end

# ==============================================================================
#  STEP 9 — WRITE P1 CIF
# ==============================================================================

function write_cif(filename::String, cif::CIFData)
    open(filename,"w") do f
        println(f,"data_zeolite_Al_substituted")
        println(f,"")
        for (k,v) in [("_cell_length_a",cif.a),("_cell_length_b",cif.b),
                      ("_cell_length_c",cif.c),("_cell_angle_alpha",cif.alpha),
                      ("_cell_angle_beta",cif.beta),("_cell_angle_gamma",cif.gamma)]
            println(f,"$k    $(@sprintf("%.6f",v))")
        end
        println(f,"")
        println(f,"_symmetry_space_group_name_H-M  'P 1'")
        println(f,"_symmetry_Int_Tables_number      1")
        println(f,"")
        println(f,"loop_")
        println(f,"_atom_site_label")
        println(f,"_atom_site_type_symbol")
        println(f,"_atom_site_fract_x")
        println(f,"_atom_site_fract_y")
        println(f,"_atom_site_fract_z")
        println(f,"_atom_site_charge")
        for i in 1:length(cif.labels)
            # Write type_symbol as-is: after classify_framework_types! these are
            # proper RASPA3 pseudo-atom names (Si, Si_b, Al, O, Oas, Ob, Hb).
            # Only fall back to element-cleaning for unclassified atoms.
            raw_elem = strip(cif.elements[i])
            if occursin("_", raw_elem) || raw_elem in ("Si","Si_a","Si_b","Al","O","Oss","Oas","Ob","Hb","He")
                type_sym = raw_elem
            else
                type_sym = uppercasefirst(lowercase(
                    replace(raw_elem, r"[^A-Za-z]+" => "")))
            end
            @printf(f,"%-10s  %-4s  %12.8f  %12.8f  %12.8f  %8.4f\n",
                    cif.labels[i], type_sym,
                    cif.frac_x[i], cif.frac_y[i], cif.frac_z[i], cif.charges[i])
        end
    end
    println("\n  Written: $filename")
end

function print_composition(cif::CIFData)
    counts=Dict{String,Int}()
    for e in cif.elements
        elem=replace(uppercase(strip(e)),r"[^A-Z]+"=>"")
        counts[elem]=get(counts,elem,0)+1
    end
    println("\n  Composition:")
    for (e,n) in sort(collect(counts)); println("    $e: $n"); end
    n_si=get(counts,"SI",0); n_al=get(counts,"AL",0)
    n_al>0&&println("    Si/Al = $(@sprintf("%.1f",n_si/n_al))")
end

# ==============================================================================
#  ARGUMENT PARSING
# ==============================================================================

function parse_arguments(args)
    positional = filter(a -> !startswith(a,"--"), args)
    flags      = filter(a ->  startswith(a,"--"), args)

    length(positional) < 2 &&
        error("Usage: julia script.jl INPUT.cif SIAL [SEED] [OPTIONS]")

    input_cif    = positional[1]
    target_si_al = parse(Float64, positional[2])
    rseed        = length(positional) >= 3 ? parse(Int, positional[3]) : 42

    zoned        = false
    outer_si_al  = Inf
    f_outer      = 0.5
    al_pref      = :homogeneous
    isites       = ["Si7","Si8"]
    esites       = ["Si14","Si17","Si20","Si1"]
    compute_only = false
    supercell    = (1,1,1)

    for flag in flags
        if startswith(flag,"--zone=")
            zoned = true
            parts = split(flag[8:end],":")
            outer_si_al = parse(Float64, parts[1])
            length(parts) >= 2 && (f_outer = parse(Float64, parts[2]))
        elseif startswith(flag,"--pref=")
            s = flag[8:end]
            al_pref = s == "channel"      ? :channel      :
                      s == "intersection" ? :intersection  :
                      s == "energy"       ? :energy        :
                      s == "dempsey"      ? :dempsey       : :homogeneous
        elseif startswith(flag,"--isites=")
            isites = String.(split(flag[10:end],","))
        elseif startswith(flag,"--esites=")
            esites = String.(split(flag[10:end],","))
        elseif startswith(flag,"--supercell=")
            parts = split(flag[13:end],"x")
            length(parts) == 3 || error("--supercell format: NxMxP  e.g. 2x2x4")
            supercell = (parse(Int,parts[1]), parse(Int,parts[2]), parse(Int,parts[3]))
        elseif flag == "--compute-core"
            compute_only = true
        end
    end

    return input_cif, target_si_al, rseed, zoned, outer_si_al, f_outer,
           al_pref, isites, esites, compute_only, supercell
end

# ==============================================================================
#  MAIN
# ==============================================================================

function main()
    if length(ARGS) < 2
        println("""
Usage:  julia lowenstein_substitution.jl INPUT.cif SIAL [SEED] [OPTIONS]

Options:
  --zone=OUTER_SIAL[:FRAC]   Si/Al of the outer zone (inner computed from balance)
                              FRAC = outer volume fraction, default 0.5
  --pref=channel              Prefer Al at channel T-sites  (Na⁺ synthesis)
  --pref=intersection         Prefer Al at intersection T-sites  (TPA⁺ synthesis)
  --pref=energy               Prefer lowest-energy MFI T-sites: T14>T17>T20>T1
  --pref=dempsey              Maximise min Al-Al distance (Dempsey\'s rule)
  --isites=Si7,Si8            Intersection T-site ASU labels (default Si7,Si8)
  --esites=Si14,Si17,Si20,Si1 Energy-priority T-site labels (MFI default)
  --supercell=NxMxP           Build N×M×P supercell (e.g. 2x2x4 → ~5000 atoms)
  --compute-core              Print zone geometry only; skip substitution

Examples:
  julia lowenstein_substitution.jl MFI_SI.cif 33
  julia lowenstein_substitution.jl MFI_SI.cif 47 42 --pref=energy
  julia lowenstein_substitution.jl MFI_SI.cif 15 42 --pref=dempsey --supercell=2x2x4
  julia lowenstein_substitution.jl MFI_SI.cif 19 42 --pref=intersection
  julia lowenstein_substitution.jl MFI_SI.cif 33 42 --zone=18 --pref=channel
  julia lowenstein_substitution.jl MFI_SI.cif 33 --zone=18 --compute-core
""")
        return
    end

    input_cif, target_si_al, rseed, zoned, outer_si_al, f_outer,
        al_pref, isites, esites, compute_only, supercell = parse_arguments(ARGS)

    na, nb, nc = supercell
    basename_noext = replace(basename(input_cif), r"\.[^.]+$" => "")
    sc_tag   = (na,nb,nc) == (1,1,1) ? "" : "_$(na)x$(nb)x$(nc)"
    pref_tag = al_pref == :homogeneous ? "" : "_$(al_pref)"
    zone_tag = zoned ? "_zoned$(Int(round(outer_si_al)))" : ""
    output_cif = "$(basename_noext)_SiAl$(Int(round(target_si_al)))$(sc_tag)$(zone_tag)$(pref_tag).cif"

    println("="^62)
    println("  Löwenstein Al Substitution")
    println("="^62)
    println("  Input:       $input_cif")
    println("  Target Si/Al = $target_si_al")
    println("  Seed:        $rseed")
    zoned && println("  Zone mode:   outer Si/Al = $outer_si_al, " *
                     "f_outer = $(@sprintf("%.2f",f_outer))")
    (na,nb,nc) != (1,1,1) && println("  Supercell:   $(na)×$(nb)×$(nc)")
    al_pref != :homogeneous && println("  Al siting preference: $al_pref")
    al_pref == :intersection && println("  Intersection T-sites: " * join(isites, ", "))
    al_pref == :energy       && println("  Energy priority:      " * join(esites, " > "))
    println("")

    # Step 1: Read CIF
    println("--- Step 1: Reading CIF ---")
    cif = parse_cif(input_cif)
    println("  Cell: a=$(cif.a)  b=$(cif.b)  c=$(cif.c)")
    println("  Asymmetric unit: $(length(cif.labels)) atoms")
    print_composition(cif)

    # Step 2: Expand symmetry to P1
    println("\n--- Step 2: Expanding symmetry to P1 ---")
    symops = parse_symops(input_cif)
    println("  Found $(length(symops)) symmetry operation(s)")
    cif, asu_src = expand_symmetry(cif, symops)
    print_composition(cif)

    # Step 2b: Build supercell (if requested)
    if (na,nb,nc) != (1,1,1)
        println("\n--- Step 2b: Building $(na)×$(nb)×$(nc) supercell ---")
        asu_src = repeat(asu_src, na*nb*nc)
        cif     = make_supercell(cif, na, nb, nc)
        print_composition(cif)
    end

    # Geometric zone analysis (when zoned)
    if zoned
        println("\n--- Geometric zone analysis ---")
        t_count_approx = count(e -> replace(uppercase(strip(e)),r"[^A-Z]+"=>"")=="SI",
                                cif.elements)
        print_zone_analysis(target_si_al, outer_si_al, t_count_approx; f_outer=f_outer)
    end

    compute_only && return

    # Steps 3–7: Löwenstein substitution
    println("--- Steps 3–7: Löwenstein Substitution ---")
    al_set, actual_ratio, o_bridges =
        lowenstein_substitute!(cif, target_si_al;
                                zoned               = zoned,
                                outer_si_al         = outer_si_al,
                                inner_fraction      = 1.0 - f_outer,
                                al_preference       = al_pref,
                                asu_src             = asu_src,
                                intersection_source = isites,
                                energy_source       = esites,
                                max_attempts        = 5000,
                                seed                = rseed)

    # Step 8: Brønsted H
    println("\n--- Step 8: Placing Brønsted H atoms ---")
    h_count, bronsted_o_map, bronsted_si_map = add_bronsted_H!(cif, al_set, o_bridges)
    print_composition(cif)

    # Step 8b: Classify all framework atoms into Hill-Sauer types
    # Sets proper type_symbols (Si, Si_b, Al, O, Oas, Ob, Hb) and charges.
    println("\n--- Step 8b: Hill-Sauer type classification ---")
    classify_framework_types!(cif, al_set, o_bridges, bronsted_o_map, bronsted_si_map)

    # Step 9: Write CIF
    println("\n--- Step 9: Writing output CIF ---")
    write_cif(output_cif, cif)

    println("\n" * "="^62)
    println("  Done!  Output: $output_cif")
    println("  Si/Al  = $(@sprintf("%.1f", actual_ratio))")
    println("  NOTES:")
    println("    · Hill-Sauer types: Si, Si_b, Al, O, Oas, Ob, Hb")
    println("    · Charges self-consistent (Schröder & Sauer, JPCC 1996)")
    println("    · H at external bisector of Al–O–Si (Pelmenschikov 1992)")
    println("    · R(O–H) = 0.97 Å")
    println("    · CIF type_symbols match force_field.json pseudo-atom names")
    println("    · Energy minimization recommended before production runs.")
    println("="^62)
end

main()
