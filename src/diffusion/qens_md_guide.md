# QENS and Classical MD for Ethanol Diffusion in H‑ZSM‑5

A working document covering: (1) what QENS actually measures, (2) which observables you should pull from your MD to compare with QENS, (3) what is wrong with your current LAMMPS setup, (4) the force‑field choices and their trade‑offs, (5) initial‑configuration strategies, (6) a polished workflow, and (7) a step‑by‑step recipe for the QENS calculation from your trajectories.

---

## 1. Quasi‑elastic neutron scattering: what it is, what it measures, and where it lives

### 1.1 Physical principle

A neutron incident on a sample with wavevector **k₀** scatters into wavevector **k₁**, exchanging momentum ℏ**Q** = ℏ(**k₀** − **k₁**) and energy ℏω = E₀ − E₁ with the sample. What the spectrometer records, after subtracting empty‑can and resolution contributions, is the dynamic structure factor S(**Q**, ω).

For a powder, only |Q| matters and you obtain S(Q, ω). Two contributions live inside it:

- **Coherent**: scattering centres interfere, giving inter‑particle (collective) information — structural Bragg peaks plus collective dynamics.
- **Incoherent**: each centre scatters independently, giving the self‑correlation function — single‑particle dynamics.

For zeolite + small organic adsorbate, the dominant signal is incoherent scattering from ¹H of the adsorbate, because σ_inc(¹H) ≈ 80 b is one to two orders of magnitude larger than every other relevant cross‑section. That is the physical reason QENS is uniquely well suited to following methanol or ethanol mobility inside an essentially "transparent" silicate framework — to first approximation you only see the protons of the guest.

The dynamic structure factor decomposes around ω = 0 into three pieces:

1. an **elastic** delta peak, broadened only by instrumental resolution, from atoms that don't move on the timescale being probed;
2. a **quasi‑elastic** broadening centred at ω = 0, from translational and rotational motions whose time scales fall inside the instrument's window;
3. **inelastic** wings (vibrations), at energies typically several meV away from the elastic line.

QENS analysis is the modelling of the *quasi‑elastic* part: its energy width as a function of Q tells you *how fast* and the integrated elastic fraction (the EISF) tells you *what geometry*.

### 1.2 The observables you actually fit

For an ensemble of independently moving protons, the incoherent scattering function in the time domain is the **self‑intermediate scattering function**:

$$F_s(Q, t) = \frac{1}{N}\sum_i \big\langle e^{i \mathbf{Q}\cdot[\mathbf{r}_i(t_0+t) - \mathbf{r}_i(t_0)]} \big\rangle_{t_0}$$

For a powder sample (orientational average over **Q̂**),

$$F_s^{\text{powder}}(Q, t) = \frac{1}{N}\sum_i \Big\langle \frac{\sin(Q\,|\Delta\mathbf{r}_i(t)|)}{Q\,|\Delta\mathbf{r}_i(t)|} \Big\rangle_{t_0}$$

This is what you compute in `qens_fsqt.jl` and the formula is correct. Its temporal Fourier transform is the experimental observable:

$$S_{\text{inc}}(Q, \omega) = \frac{1}{\pi}\int F_s(Q, t)\,e^{-i\omega t}\,dt$$

Two limits of F_s(Q, t) carry the whole physical story:

- **Long‑time plateau** B(Q) = lim_{t→∞} F_s(Q, t). Localized (rotational or confined‑translational) motions never decorrelate completely, so F_s decays to a non‑zero plateau. Plotted vs Q this plateau is the **Elastic Incoherent Structure Factor (EISF)**, A₀(Q). Its Q‑dependence is a fingerprint of geometry: isotropic rotation, methyl 3‑site rotation, jump‑in‑a‑sphere, etc.
- **Initial decay rate** Γ(Q). For unbounded translational diffusion, Γ(Q) is the HWHM of the Lorentzian S_inc(Q, ω) and grows with Q² (Fickian, slope = D_s) or saturates at high Q (jump diffusion: Chudley‑Elliott or Hall‑Ross).

### 1.3 Instruments — how the time/length window is chosen

QENS spectrometers come in two geometries:

**Direct geometry** (e.g. IN5 at ILL, LET at ISIS): a chopper selects the incoming neutron energy; the outgoing energy is measured by time of flight. You set E_i, which fixes the resolution and the (Q, ω) window.

**Indirect geometry / backscattering** (e.g. OSIRIS, IRIS at ISIS; IN16B at ILL): the analyser crystals select a fixed *final* energy, so the incoming wavelength is the variable and you scan ω by varying time‑of‑flight from source to sample. The PG(002) analyser on OSIRIS/IRIS gives ΔE ≈ 25 µeV resolution.

| Spectrometer | ΔE (µeV) | Time window (ps) | Length window (Å) |
|---|---|---|---|
| IRIS (PG002) | ~17 | ~3 – 440 | 0.4 – 4 |
| OSIRIS (PG002) | ~25 | ~1 – 100 | 0.4 – 4 |
| IN5 (ILL) | 10 – 100 (depends on E_i) | up to ~300 | 0.1 – 4 |
| IN16B (ILL) | ~0.7 | ~50 – 5000 | 0.1 – 2 |

Key practical point that the Matam 2023 paper hammers home: **the same molecule in the same zeolite can look like rotation on one instrument and translational diffusion on another**, because each instrument simply doesn't see what is too fast or too slow for its window. You should always ask what timescale you are probing with a given configuration before you read off "the" diffusion mechanism.

### 1.4 Sample requirements

QENS samples are typically:
- a few grams of powder packed in an annular Al can (to keep transmission ≈ 90% and avoid multiple scattering);
- pre‑activated zeolite, then loaded with a known amount of (preferably deuterated‑framework / hydrogenated‑guest) adsorbate so that essentially all the incoherent scattering is from the guest;
- mounted on a closed‑cycle refrigerator or furnace stick (typical T range 4 – 700 K).

A "background" run with empty zeolite and a "resolution" run at low T (~10 K) of the loaded sample are essential — the latter defines what counts as "elastic".

### 1.5 What QENS cannot do

- It cannot resolve macroscopic D_s much below ~10⁻¹² m²/s on backscattering instruments — the Lorentzian becomes narrower than the resolution function and disappears into the elastic peak. Spin‑echo (NSE) extends this by another decade.
- It is not single‑molecule sensitive; you always work with ensemble‑averaged self‑correlations.
- It does not by itself separate the dynamics of two different hydrogenous species (e.g. methanol vs water): all ¹H contribute to the same S_inc(Q, ω). Selective deuteration is the only experimental fix; MD is the natural complement.

---

## 2. From QENS observables to MD analyses

The whole reason QENS pairs so naturally with classical MD is that the **time and length scales the two techniques cover are essentially identical**: 1 ps – 1 ns and 1 – 20 Å. MD gives you all atomic coordinates as a function of time, so any S(Q, ω) you can measure you can also compute.

### 2.1 Three families of MD analyses, and what each one corresponds to

| MD quantity | QENS observable | Physical info |
|---|---|---|
| MSD ⟨\|Δr(t)\|²⟩ → Einstein D_s | Slope of HWHM vs Q² at low Q (Fickian) | Long‑range translational diffusivity |
| Direct F_s(Q, t) from atomic trajectories | F_s(Q, t) (NSE) or its Fourier transform S_inc(Q, ω) (TOF/BS) | Full QENS lineshape, jump diffusion fits |
| Long‑time plateau of F_s(Q, t) vs Q | EISF A₀(Q) | Geometry: isotropic / 3‑site / sphere / etc. |
| Rotational F_s computed from r_H(t) − r_COM(t) | High‑Q broadening, EISF | Rotational dynamics decoupled from translation |

The recipe in the Armstrong 2020 and Matam 2023 reviews is exactly this: compute F_s(Q, t) directly from atomic displacements, fit it as a sum of one or two exponentials, and read off the corresponding Lorentzian widths Γ_i = 1/τ_i. The constant baseline B(Q) is the EISF.

### 2.2 What kind of MD do you need

For ethanol diffusion in H‑ZSM‑5 the right tool is **classical MD with a flexible all‑atom force field**, in the NVT ensemble (Nosé‑Hoover or CSVR thermostat) for equilibration, and either NVT or NVE for production. Reasons:

- **Timescales**. You need ≥ 1 ns of production for translational diffusion to be in the diffusive regime in MFI at 373 K (the Woodward paper used 2 ns; Dunn et al. used 10 ns for esters and 100 ns for slower aldehydes). AIMD cannot reach this for a 2×2×2 ZSM‑5 supercell.
- **Sampling**. You will compute time‑lag averages of F_s(Q, t) over many origins t₀ and over many adsorbate molecules. ~5–10 independent runs of ~5 ns each is a sensible target.
- **Framework treatment**. A flexible framework is the modern standard (Dunn 2024, the recent O'Malley group papers) and is needed for honest QENS comparison — a frozen framework changes both quantitative diffusivities and the local potential that the guest sees, especially near the BAS. The Catlow rigid‑ion potentials let the framework breathe at the right amplitude for THz vibrations and do not require freezing.

### 2.3 The single most important MD setup choice for QENS comparison

Match the instrument's resolution to your dump frequency. If you want to compare to an OSIRIS spectrum (ΔE ≈ 25 µeV ↔ ~50 ps time resolution, energy window ±0.5 meV ↔ ~4 ps minimum time resolution), then:

- dump every Δt_dump ≤ 0.5 ps (ideally 0.1–0.2 ps), so the high‑frequency end of S(Q, ω) is captured;
- run for ≥ 1 ns of production *after* equilibration, so the low‑frequency end ω ~ 25 µeV (~165 ps decay) is resolved;
- have ≥ 50 protons per molecule × number of molecules so the noise on F_s(Q, t) is acceptable.

This is the simulation analogue of the experimentalist's choice between IRIS and IN16B: pick a simulation length and dump frequency that brackets the instrument window you care about.

---

## 3. Diagnosis of your current LAMMPS run

### 3.1 Why the trajectory is "frozen": three coupled problems

**Problem 1 (the lethal one): `fix rmom` is never unfixed.**

```lammps
# in Phase 2
fix rmom ethanol momentum 1000 linear 1 1 1 angular
...
# Phase 3 — no `unfix rmom`
# Phase 4 — no `unfix rmom`
```

`fix momentum` zeroes the linear and angular momentum of the group every 1000 steps **without rescaling velocities**. With a *frozen* framework, every collision an ethanol molecule has with the wall is asymmetric — momentum is transferred into the (force‑zeroed) framework and lost. The drift in the ethanol's COM is then explicitly killed by `fix momentum` every 0.5 ps. Combined, this is a slow, monotonic energy drain. The thermo log shows it directly: T_ethanol falls from 376 K → 18 K over 5 ns, KE from 1008 → 49 kcal/mol. The MSD growth flattens to zero accordingly.

**Problem 2: Phase 4 is pure NVE with no thermostat to compensate.**

Even without `fix rmom`, NVE on a small group inside a frozen wall will drift in temperature because the wall absorbs random momentum. Phase 4 should be either a long NVT (Nosé‑Hoover with a τ ~ 100 fs that is gentle enough not to bias the dynamics, as Dunn et al. show in their thermostat tests) or a flexible‑framework NVE.

**Problem 3 (under the hood): the framework is rigid by construction.**

`fix freeze framework setforce 0.0 0.0 0.0` plus `velocity framework set 0.0 0.0 0.0` removes all framework dynamics. This is not, by itself, the cause of T crash — the energy crash is the rmom + NVE combination. But a rigid framework changes the diffusion physics in ways that matter for your comparison to QENS:
- the BAS Hb cannot reorient or hydrogen‑bond hop;
- the framework does not provide thermal coupling to the guest;
- the 10‑MR pore "windows" do not breathe, so jump rates between intersections are systematically wrong.

The Dunn 2024 paper explicitly argues against frozen frameworks (Section 1: *"it is not possible to know a priori how zeolite flexibility will affect diffusion, so results with a fixed framework cannot be extrapolated to one which is flexible"*).

### 3.2 Other things to fix in the input

- `comm_modify cutoff 20.0` is unnecessarily large for a 10 Å cutoff and inflates ghost‑atom storage. With a flexible framework the default + `neighbor 2.0 bin` is fine; with a rigid framework you don't need a neighbor list rebuild on framework atoms at all.
- `compute Teth ethanol temp` is correct, but `thermo_modify temp Teth` then makes the *thermodynamic pressure* compute use ethanol's temperature, which is wrong (pressure is ill‑defined with a frozen subsystem anyway). The pressure in your log is meaningless — that's why the warning fires. Just don't print pressure if framework is frozen.
- The dump frequency `${dumpfreq} 4000` = 2 ps at 0.5 fs is **far too coarse for QENS analysis**. You will lose all motion faster than ~4 ps (i.e., you cannot resolve anything above ~1 meV). For QENS analysis, dump xu yu zu of the ethanol H atoms every 100–200 fs (every 200–400 steps).

### 3.3 The ethanol initial configuration

80 ethanol molecules in a 2×2×2 ZSM‑5 supercell = **10 ethanol per primitive unit cell**. That is roughly twice the highest loading in any of the papers I've read for MFI (Woodward et al.: 3 and 5 methanol/uc; Dunn et al.: 1 promoter/uc + 2 methanol/promoter; max plausible ethanol loading from saturation isotherms is ~6 per uc). At 10/uc the random insertion is bound to leave bad contacts that energy minimisation alone cannot fully heal. Even if you fix the rmom bug, the initial PE will be enormous and the first ps of dynamics will be a shockwave.

---

## 4. Force field — assessment and recommendations

### 4.1 Your current choice (Woodward‑style framework + OPLS ethanol)

Reading `forcefield_hzsm5_ethanol.in` against the Woodward 2022 SI:

- **Framework charges (Hill–Sauer style)**: Si = +4, Al = +3, O = −2, Ob = −1.426, Hb = +0.426. Per Al substitution the net charge change is (Al−Si) + (Ob_change) + Hb = −1 + 0.574 + 0.426 = 0. ✓ Consistent. The "net charge 0.01" warning is the rounding error from the file's precision; harmless.
- **Buckingham A, ρ, C** for Si–O, Al–O, O–O: match the Woodward SI Table SI1 after eV → kcal/mol conversion. ✓
- **Ob–Hb Morse** with D = 162.6 kcal/mol, α = 2.20, r_e = 0.985 Å: matches the Woodward SI. ✓
- **OPLS bonded for ethanol**: numbers from LigParGen for ethanol look right (the charges sum to 0 to 4 decimal places).
- **Ethanol–framework cross LJ**: you adapted the methanol‑zeolite cross terms from the Woodward SI. The σ_OC ≈ 4.31 Å looks unphysically large compared to Lorentz–Berthelot but this is a deliberate Kiselev‑style "hard‑wall" parameterisation; the Woodward / O'Malley group uses it consistently, so as long as you adopt the *same* approach for ethanol, you are internally consistent.

So as a force field this is internally consistent and a defensible choice for varying Si/Al ratios — its key advantage over the Catlow potentials is precisely that it parameterises Al, Ob and Hb explicitly, which Catlow does not. **Keep it.**

### 4.2 The Catlow rigid‑ion alternative (Dunn 2024 style)

The Catlow potentials are:
- Si–O Buckingham (A = 1283.9 eV, ρ = 0.32, C = 10.66 eV·Å⁶ in the Sanders/Catlow form);
- O–O Buckingham (A = 22764.0 eV, ρ = 0.149, C = 27.88 eV·Å⁶);
- formal charges Si = +4, O = −2;
- harmonic O–Si–O bending (this is what makes the framework well‑behaved at room T).

These are excellent for siliceous frameworks and have been used at scale, but they have **no Al or BAS terms**. To extend them to H‑ZSM‑5 you would need to graft on Al/Ob/Hb interactions consistent with the rest, which is exactly what the Hill–Sauer / Woodward parameterisation already does. So unless you only want to study siliceous MFI, your current Woodward force field is the right one.

### 4.3 If you want to harden the FF later

- Replace the O–H₂O‑style cross LJ with the explicit Boronat/Corma or DREIDING‑like parameters used in the Cnudde / Van Speybroeck papers, if you ever go to ab initio MD comparison.
- Cross‑validate by computing methanol diffusivities at the same loadings and Si/Al ratios as Woodward et al. (Table 1 of their paper). If your D_s values land within a factor of 2 of theirs at 373 K, the FF is calibrated.

---

## 5. Initial configuration: stop using random insertion

Random insertion at 80 ethanol/2×2×2 supercell guarantees bad contacts. Three saner options, in order of increasing fidelity:

### 5.1 Reduce the loading to literature values

For QENS comparison the loading should match an actual experimental loading. Pick one of:
- **3 ethanol/uc** = 24 ethanol in 2×2×2 (Woodward "low loading")
- **5 ethanol/uc** = 40 ethanol in 2×2×2 (Woodward "high loading", most clusters formed)
- **n_BAS × 2** = 2 ethanol per Brønsted site, mirroring a "saturation" loading per site

At Si/Al = 19 in a 2×2×2 ZSM‑5 (96 T sites total, so 5 Al substitutions ≈ Si/Al = 19), 2 ethanol/BAS = 10 ethanol — way more tractable than 80.

### 5.2 GCMC pre‑equilibration (recommended)

A grand‑canonical Monte Carlo equilibration in RASPA (or LAMMPS GCMC) at the desired chemical potential / loading gives you:
- an equilibrated ethanol configuration that respects the framework topology,
- correct adsorption sites at the channel intersections / BAS sites,
- no overlapping atoms by construction.

Workflow:
1. Convert your LAMMPS framework data to RASPA CIF + PDB.
2. Run NVT‑GCMC at 373 K with the same Lennard‑Jones and Coulomb parameters until the loading converges to your target.
3. Export the final ethanol positions and **convert the united‑atom or all‑atom GCMC outputs back to your LAMMPS all‑atom representation**: RASPA can run all‑atom OPLS ethanol natively, so you don't need the UA → AA conversion if you set RASPA up with all‑atom ethanol from the start. If you do start from UA (CH3, CH2, OH pseudo‑atoms), the safe thing is to take the COM of each molecule and the orientation of the C–C–O backbone and rebuild the all‑atom ethanol on top of that with reasonable hydrogen positions, then minimise.

### 5.3 Sequential insertion + relaxation in LAMMPS

If you don't want to add RASPA to your toolchain, the cleanest pure‑LAMMPS option is:

```
fix gcmc ethanol gcmc 100 1 0 0 ${seed} ${T} -1.0 0.0 mol et_template
```

This inserts/removes ethanol molecules in the canonical/grand‑canonical ensemble using LAMMPS's `fix gcmc`. The template `et_template` is a single ethanol molecule from `ethanol.lmp`. You let it run until you have the target N_ethanol, then turn off `fix gcmc` and run normal MD. This is much cleaner than your random‑insertion script because it accepts/rejects on an actual Boltzmann criterion, so you never build an overstrained configuration.

### 5.4 What I'd actually recommend

Given where you are: 

1. **Start over with 5 ethanol per primitive unit cell = 40 ethanol total**, matching the Woodward 5‑MPUC condition. This gives you a direct literature comparison.
2. Either use `fix gcmc` to insert them, or extend your Julia script to use a much larger heavy‑atom cutoff (3.0 Å) and accept fewer molecules per attempt while it relaxes between batches.
3. After insertion, do a long, careful minimization: `min_style cg` then `min_style hftn`, to dmax = 0.05 Å, with the framework free.

---

## 6. Polished workflow

Below is the order of operations I'd run, with the rationale at each step. The corrected LAMMPS input is in §A1 at the end of the document.

### Step 1 — Build a flexible H‑ZSM‑5 supercell

You already have this: `mfi_sial_19.data` is a 2×2×2 supercell with 5 Al substitutions at the lowest‑energy T sites, BAS at Ob positions, and full bonded topology (3112 bonds / 6224 angles / 9456 dihedrals / 40 impropers). The full topology is what makes the framework *flexible* under MD.

Do NOT drop the framework topology in the build script when you go to a flexible run. Your current script (`build_hzsm5_ethanol.jl`) explicitly drops it. Modify the script to keep both framework and ethanol topology and concatenate them with offset bond/angle/dihedral type IDs — this is the bigger change.

### Step 2 — Insert ethanol at literature loading

Use `fix gcmc` (LAMMPS) or RASPA. Target loading: 5 ethanol/uc → 40 molecules in 2×2×2.

### Step 3 — Minimise

Conjugate gradient with `dmax = 0.02`, then `hftn` to tighten. With a flexible framework this is essential to remove insertion artefacts everywhere, not just on the ethanol.

### Step 4 — Heat & equilibrate

- 50 ps of NVT at 50 K with `nve/limit 0.05` to safely thaw without explosions (the `nve/limit` damps any over‑energetic ethanol that was inserted into a bad spot).
- 200 ps of NVT (Nosé‑Hoover, τ_T = 100 fs) ramping from 50 K to 373 K.
- 200 ps of NVT at 373 K, full system thermostatted.

Use **CSVR (`fix temp/csvr`)** rather than Berendsen for production NVT — it has the correct canonical ensemble distribution. Berendsen is OK for the warmup ramp but biases velocities at long times.

### Step 5 — Production

Two sensible choices:

- **NVT with CSVR**, τ_T = 100 fs, 5 ns. This is the cleanest for QENS comparison (canonical, T‑controlled, no energy drift).
- **NVE** *only if* you want to verify the canonical result and your equilibrated config has very stable energy. Always check `etotal` is stationary, not drifting.

Repeat 5× with different velocity seeds for proper error bars on D_s and F_s(Q, t).

### Step 6 — Dumps for analysis

Two separate dumps:

```lammps
# Dump 1 — sparse, for MSD / D_s and visualisation
dump  dprod  ethanol  custom 1000  traj_prod.lammpstrj  id mol type xu yu zu
# 1000 steps × 0.5 fs = 0.5 ps cadence, 5 ns ÷ 0.5 ps = 10000 frames

# Dump 2 — dense, ethanol H atoms only, for F_s(Q, t)
group ethanol_H type 12 13 14 15 16 17
dump  dqens  ethanol_H custom 200 traj_ethanol_H_qens.lammpstrj id mol type xu yu zu
# 200 steps × 0.5 fs = 0.1 ps cadence, 5 ns ÷ 0.1 ps = 50000 frames
```

The 0.1 ps cadence resolves OSIRIS energies up to ~20 meV, which is more than enough.

### Step 7 — Analyse (next section)

---

## 7. Step‑by‑step QENS calculation from your trajectory

Your `qens_fsqt.jl` script is mathematically right. Here is how to use it correctly together with the rest of the analysis pipeline; what your output files mean; and how to fit them.

### 7.1 The MSD and D_s (the easy one)

From `traj_prod.lammpstrj`:

1. Extract centre of mass of each ethanol molecule (use the C atoms or compute mass‑weighted COM):
$$\mathbf{R}_m(t) = \frac{1}{M_{eth}}\sum_{i \in m} m_i \mathbf{r}_i(t)$$
2. Compute the lag‑averaged MSD:
$$\mathrm{MSD}(t) = \frac{1}{N_{mol}\,N_{t_0}}\sum_{m, t_0} |\mathbf{R}_m(t_0+t) - \mathbf{R}_m(t_0)|^2$$
3. Discard the first 100–200 ps (sub‑diffusive regime).
4. Linear fit MSD(t) = 6 D_s t + offset. The slope ÷ 6 is D_s in Å²/ps. Convert to m²/s by × 10⁻⁸.
5. Repeat per Cartesian direction for D_xx, D_yy, D_zz to expose anisotropy (in MFI the b axis = straight pores has higher D, exactly as Dunn et al. show).

A clean way to do this in `kinisi` (the Bayesian package the Dunn paper used) is:

```python
from kinisi.parser import LammpsParser
from kinisi.diffusion import DiffusionAnalyzer
p = LammpsParser('traj_prod.lammpstrj', specie='C')      # or per-molecule COM
da = DiffusionAnalyzer.from_universe(p.universe, time_step=0.5, step_skip=1000)
da.diffusion(start_dt=200)                               # discard first 200 ps
print(da.D, '±', da.D_offset)                            # m²/s
```

`kinisi` gives you posterior distributions on D_s, not just point estimates — which is the right thing to report.

### 7.2 The self‑intermediate scattering function F_s(Q, t)

This is the central QENS observable. Your script computes it correctly. The full pipeline:

**(a) Choose Q values** that match the experimental Q range. OSIRIS/IRIS has Q ≈ 0.2 – 1.85 Å⁻¹. A good grid: Q ∈ {0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 1.9} Å⁻¹.

**(b) Choose lag range and frame cadence**. With a 0.1 ps dump cadence and a 5 ns trajectory, you can compute F_s(Q, t) for lags from 0 to ~1000 ps with good statistics (cap your max lag at ~25% of trajectory length to keep the t₀ average meaningful).

**(c) Run `qens_fsqt.jl`**:

```bash
julia qens_fsqt.jl traj_ethanol_H_qens.lammpstrj \
  --dt-frame-ps 0.1 \
  --discard-ps 200 \
  --max-lag-ps 1000 \
  --q 0.3 0.5 0.7 1.0 1.3 1.6 1.9 \
  -o fsqt.csv
```

A note on the script's complexity: it is currently O(N_lag × N_t0 × N_atom × N_q) which is fine for small lags but explodes if you push max_lag to 1000 ps with 50000 frames. Two options: (i) reduce N_t0 by striding over t₀ (every 10th frame is plenty for a 50000‑frame trajectory); (ii) port to FFT‑based computation, which is O(N log N) per atom per Q, using e.g. `numpy.fft` in Python. For your first runs the direct method is fine.

**(d) Fit F_s(Q, t) to a sum of exponentials plus a baseline**:

$$F_s(Q, t) = A_1(Q) e^{-\Gamma_1(Q) t} + A_2(Q) e^{-\Gamma_2(Q) t} + B(Q)$$

You typically need at most two exponentials in zeolite work — one fast (rotational, ~ps) and one slow (translational, ~10–100 ps). The plateau B(Q) is the EISF.

A simple Python fit:

```python
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

df = pd.read_csv('fsqt.csv')
def biexp(t, A1, G1, A2, G2, B):
    return A1*np.exp(-G1*t) + A2*np.exp(-G2*t) + B

t = df['t_ps'].values
results = []
for col in df.columns[1:]:
    Q = float(col.split('_')[2])
    fs = df[col].values
    popt, _ = curve_fit(biexp, t, fs, p0=[0.3, 1.0, 0.4, 0.05, 0.3],
                        bounds=([0,0,0,0,0],[1,100,1,10,1]))
    results.append({'Q': Q, 'A1': popt[0], 'G1_1/ps': popt[1],
                    'A2': popt[2], 'G2_1/ps': popt[3], 'EISF': popt[4]})
print(pd.DataFrame(results))
```

### 7.3 The EISF and its geometric interpretation

Plot B(Q) vs Q. Compare against analytical models (all in the Matam 2023 paper):

- **Free isotropic rotation** of the methyl group around its axis: A₀(Q) = (1/3)[1 + 2 j₀(Q a √3)] with a = C–H bond projected length ~1.02 Å.
- **Confined diffusion in a sphere** of radius R: A₀(Q) = [3 j₁(QR) / (QR)]².
- **Long‑range translation**: B(Q) → 0 at all Q (no plateau).
- **Mixed model**: p_mobile × A₀_localised(Q) + (1 − p_mobile) × 1, where p_mobile is the mobile fraction.

In your case ethanol has **two rotational modes** to think about: the methyl group (CH₃) rotation and the OH rotation. Each contributes a different A₀(Q) shape; you may need a sum.

### 7.4 The Lorentzian widths and jump diffusion

Plot Γ_slow(Q) vs Q². Three regimes you should be able to distinguish:

- **Fickian translation** (low Q): Γ ∝ Q² with slope D_s. Use the same D_s here as you got from the MSD — they should agree to ~20%.
- **Jump diffusion** (intermediate Q): Γ saturates following Chudley‑Elliott Γ(Q) = (1/τ)[1 − sin(Qd)/(Qd)], where d is the jump length and τ is the residence time. Fit d, τ; D_s is recovered as d²/(6τ).
- **Confined diffusion** (low Q, in sphere of radius R): Γ is Q‑independent at low Q and equals D_s/R² (Volino–Dianoux model).

The Q‑dependence pattern of Γ_slow is the diagnostic for which regime ethanol is in at your loading and Si/Al — exactly the analysis Matam 2023 walk through for methanol.

### 7.5 Sanity checks before publishing anything

- **Energy conservation**: in the production run, total energy should drift by less than ~0.1% of KE per ns. If it drifts more, your timestep is too long or your LJ cutoff is too tight.
- **Temperature stationarity**: T_inst should fluctuate around T_target with σ_T ≈ T × √(2/(3N)). For 40 ethanol × 9 atoms = 360 atoms, that's ~16 K at 373 K. If the running average drifts, your thermostat is wrong.
- **MSD ↔ F_s(Q, t) consistency**: at small Q, the slope of −d log F_s/dt at t > 100 ps must equal D_s × Q² where D_s is the MSD‑derived diffusion coefficient. This is the strongest internal consistency check.
- **Resolution convolution**: if you want a *direct visual* match to the experimental S(Q, ω), Fourier‑transform F_s(Q, t) (or take the analytic Lorentzian transform of your fit), then convolve with a Gaussian of FWHM equal to the instrument resolution (25 µeV for OSIRIS). MDANSE does this automatically; in Python a 4‑line scipy.signal.fftconvolve does the same.

---

## A1. A corrected LAMMPS input (drop‑in replacement)

This keeps the spirit of your file but fixes the rmom bug, switches to a flexible framework with proper thermostat, sets sensible dump cadences for QENS, and removes the meaningless pressure printing.

```lammps
# H-ZSM-5 + ethanol — flexible-framework production run
units           real
atom_style      full
boundary        p p p
newton          on

bond_style      harmonic
angle_style     harmonic
dihedral_style  opls
improper_style  cvff

variable T          equal 373.0
variable dt         equal 0.5
variable seed       equal 20260424

# read a data file that retains BOTH framework and ethanol topology
read_data       data.hzsm5_sial19_ethanol_full
include         forcefield_hzsm5_ethanol_flexible.in
special_bonds   lj/coul 0.0 0.0 0.5

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

timestep        ${dt}

group           framework type 1 2 3 4 5 6 7 8
group           ethanol   type 9 10 11 12 13 14 15 16 17
group           ethanol_H type 12 13 14 15 16 17

# --- Minimisation ---
min_style       cg
min_modify      dmax 0.02
minimize        1.0e-4 1.0e-6 20000 200000
write_data      data.minimised

# --- Stage 1: thaw at 50 K with displacement limit ---
velocity        all create 50.0 ${seed} mom yes rot yes dist gaussian
fix             nvt1 all nve/limit 0.05
fix             tstat1 all temp/csvr 50.0 50.0 100.0 ${seed}
run             100000          # 50 ps
unfix           nvt1
unfix           tstat1

# --- Stage 2: heat 50 K -> 373 K over 200 ps ---
fix             nvt2 all nvt temp 50.0 ${T} 100.0
run             400000          # 200 ps
unfix           nvt2

# --- Stage 3: equilibrate at 373 K for 200 ps ---
fix             nvt3 all nvt temp ${T} ${T} 100.0
run             400000
unfix           nvt3

# --- Production: NVT/CSVR, 5 ns ---
reset_timestep  0
fix             prod all nvt temp ${T} ${T} 100.0
thermo          1000
thermo_style    custom step temp pe ke etotal vol
thermo_modify   flush yes

# Sparse dump for MSD / visualisation
dump            d1 ethanol custom 1000 traj_prod.lammpstrj id mol type xu yu zu
dump_modify     d1 sort id

# Dense dump for QENS F_s(Q, t)
dump            d2 ethanol_H custom 200 traj_ethanol_H_qens.lammpstrj id mol type xu yu zu
dump_modify     d2 sort id

# MSD on ethanol C atoms
group           ethanol_C type 9 10
compute         msd ethanol_C msd com yes
fix             msdout all ave/time 200 1 200 c_msd[1] c_msd[2] c_msd[3] c_msd[4] file msd.dat mode scalar

run             10000000         # 5 ns at 0.5 fs
write_restart   restart.production.${T}K
write_data      data.final.${T}K
```

Two key things to notice vs your original:

1. **No `fix momentum`**. In NVT with a fully flexible system there is no need to remove COM motion because the thermostat couples to all atoms.
2. **All atoms are integrated, framework included**. The Hill–Sauer / Woodward force field with framework bonds, angles, BAS Morse term, and Buckingham non‑bonded gives you a *physical* flexible framework at room temperature.

If you really want to keep a frozen framework for a smoke‑test, the only changes from this template are:

- replace `fix prod all nvt ...` with `fix prod ethanol nvt ...` and add `fix freeze framework setforce 0.0 0.0 0.0`;
- explicitly `unfix` any `fix momentum` you add (or, better, never add one);
- accept that quantitative comparison to QENS is approximate.

---

## A2. What to compare against in the literature

When you have your first F_s(Q, t) and D_s for ethanol/H‑ZSM‑5, sanity‑check against:

| System | T | D_s (m²/s) | Source |
|---|---|---|---|
| Methanol / siliceous MFI, 5 mpuc | 373 K | 2.3 × 10⁻¹⁰ | Woodward 2022 |
| Methanol / H‑ZSM‑5 Si/Al = 15, 5 mpuc | 373 K | 0.65 × 10⁻¹⁰ | Woodward 2022 |
| Methanol / H‑ZSM‑5 Si/Al = 95, 5 mpuc | 373 K | 3.06 × 10⁻¹⁰ | Woodward 2022 |
| Methanol / H‑ZSM‑5 Si/Al = 140 | 373 K | 1.6 × 10⁻⁹ | Matam 2023 (QENS) |
| Methanol / H‑ZSM‑5 Si/Al = 30 | 300 K | ~10⁻¹¹ | Jobic 1986 (QENS, IN10) |
| Methyl ester promoters / siliceous MFI | 423 K | 0.5–6 × 10⁻⁵ cm²/s ~ 10⁻⁹ m²/s | Dunn 2024 |

Your ethanol values should fall slightly below methanol (ethanol is bigger), so expect D_s ~ 10⁻¹⁰ – 10⁻⁹ m²/s in the same loading range, with the same Si/Al trend.
