# JZeoMCMD

[![Build Status](https://github.com/viniciuspiccoli/JZeoMCMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/viniciuspiccoli/JZeoMCMD.jl/actions/workflows/CI.yml?query=branch%3Amain)


### check parameters and conversions - add original papers for the force fields
### check analysis and meaning of them
### check scripts that handle the transformations

### add different Si/Al structure with 2x2x2 supercell and their data file from ovito
### need to add aluminossilicate to reload_adsorbate.jl



Add baseline unit tests. (OK)
Add ExternalInputFiles. (OK)
Add validate_inputs. ( )
Add generic file staging. ( )
Add a new workflow overload using external files. ( )
Run it with the current silicalite-1 inputs. ( )
Run it with the current Si/Al = 19 inputs. ( )
Move ethanol and MFI metadata out of the generic layer. ( )
Replace restart parsing. ( )
Introduce time-series and block statistics. ( )
Replace convergence criteria. ( )
Deprecate internal scientific input directories. ( )





┌─────────────────────────────────────────────────┐
│  Cycle 1 (initial)                              │
│                                                 │
│  test.data (Ovito)                              │
│       ↓  build_loaded_zeolite.jl                │
│  loaded_zeolite.lmp + run_loaded.in             │
│       ↓  LAMMPS NPT-MD                          │
│  loaded_npt_final.lmp                           │
│       ↓  write_cif.jl                           │
│  relaxed_MFI.cif                                │
└────────────────┬────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────┐
│  Cycle N (iterative)                            │
│                                                 │
│  relaxed_MFI.cif → RASPA3 GCMC → restart.json   │
│       ↓                              ↓          │
│       │    reload_adsorbate.jl ←─────┘          │
│       │    (reads loaded_npt_final.lmp          │
│       │     + new restart.json)                 │
│       ↓                                         │
│  cycleN_loaded.lmp                              │
│       ↓  LAMMPS NPT-MD (same run_loaded.in)     │
│  cycleN_npt_final.lmp                           │
│       ↓  write_cif.jl                           │
│  cycleN_relaxed.cif → RASPA3 GCMC → ...         │
└─────────────────────────────────────────────────┘




# what to do:

check: force_field.json
check charges and lennard jones (check why using lj parameters using Bai.). Why Si for silicate and Si for aluminossilicate are different?
check params.toml and params_loader.jl
adjust the write_lammps_data functions in: reload_adsorbate.jl, build_loaded_zeolite.jl, add_zeolite_topology.jl
check the paramters in generate_nb_tables.jl

add this: nohup ./your_script.sh > my_output.log 2>&1 &




base_dir/                         ← user provides these:
├── MFI_SI.data                   ← supercell data (any name)
├── MFI_SI.cif                    ← unit cell CIF (any name)
├── hillsauer_nb.table            ← NB table
├── run_npt.in                    ← LAMMPS input
└── cycle_01/...                  ← workflow creates these

the package provides (from raspa_inputs/):
force_field.json                  ← H-S charges + Bai LJ + TraPPE
ethanol.json                      ← TraPPE-UA molecule definition
simulation.json.template          ← cycle 1 (2×2×2 unit cells)
simulation.json.template_next     ← cycle 2+ (1×1×1 supercell)


wp = WorkflowParams(...)
wp.initial_cif  = "MFI_SI.cif"      # ← in base_dir
wp.initial_data = "MFI_SI.data"      # ← in base_dir
run_gcmc_md_workflow(wp)


remove json files from the inside also!!!!


src/
├── JZeoMCMD.jl                    ← updated (new includes + exports)
├── read_lammps_data.jl            ← unchanged
├── add_zeolite_topology.jl        ← fixed minmax bug
├── refine_topology_types.jl       ← NEW
├── alumino_support.jl             ← NEW (replaces old version)
├── build_loaded_zeolite.jl        ← updated (configurable eth topology defs)
├── reload_adsorbate.jl            ← updated (silica/alumino dispatch)
├── workflow.jl                    ← updated (is_alumino flag)
├── write_cif.jl                   ← unchanged
├── analysis.jl                    ← unchanged
├── params_loader.jl               ← unchanged
├── generate_nb_tables.jl          ← unchanged
├── setup_pressure_sweep.jl        ← unchanged
└── ff/
    ├── hillsauer_alumsil_empty.ff   ← NEW
    └── hillsauer_alumsil_loaded.ff  ← NEW



    # JZeoMCMD

**A Julia package for coupled GCMC/NPT-MD simulations of molecular adsorption in flexible zeolites.**

JZeoMCMD automates the iterative grand canonical Monte Carlo (GCMC) / isothermal–isobaric molecular dynamics (NPT-MD) workflow for studying adsorption-induced deformation in zeolite frameworks. It supports both all-silica and aluminosilicate (H-form) zeolites with ethanol as the adsorbate, using the Hill–Sauer / mHSFF force field for framework flexibility.

---

## Table of Contents

- [Physical Context](#physical-context)
- [Systems Under Study](#systems-under-study)
- [Force Field Description](#force-field-description)
- [Coupled GCMC/NPT-MD Workflow](#coupled-gcmcnpt-md-workflow)
- [Simulation Details](#simulation-details)
- [Convergence and Statistics](#convergence-and-statistics)
- [Code Architecture](#code-architecture)
- [Installation](#installation)
- [Input Files](#input-files)
- [Usage Examples](#usage-examples)
- [References](#references)

---

## Physical Context

Zeolites are microporous aluminosilicate minerals with well-defined channel and cavity systems at the molecular scale. Their frameworks are not rigid: upon adsorption of guest molecules, the unit cell parameters change — a phenomenon known as **adsorption-induced deformation** (AID). This deformation feeds back into the adsorption thermodynamics: a swollen or contracted pore accommodates a different number of molecules than a rigid one.

Capturing this coupling requires simultaneous treatment of adsorption equilibrium (grand canonical ensemble) and framework flexibility (isothermal–isobaric ensemble). Neither GCMC alone (which fixes the framework) nor NPT-MD alone (which fixes the chemical potential) can capture this self-consistently. The solution is an **iterative GCMC/NPT-MD scheme** where the two simulations exchange structural information until convergence (Daou et al. 2021; Emelianova et al. 2023).

The physical observables of interest include adsorption isotherms, cell parameter evolution with loading, volumetric strain, channel-specific occupancy (straight vs. sinusoidal vs. intersections for MFI), and diffusion coefficients.

---

## Systems Under Study

### MFI Topology (ZSM-5 / Silicalite-1)

MFI is an orthorhombic zeolite (space group *Pnma*) with a three-dimensional channel system consisting of straight channels along the *b*-axis (5.3 x 5.6 A), sinusoidal channels along the *a*-axis (5.1 x 5.5 A), and intersections where the two channel systems cross (~9 A cavity).

<!-- Add images of your zeolite structures:

![MFI framework along b-axis](docs/images/mfi_b_axis.png)
*Figure 1: MFI framework viewed along the b-axis showing straight channels.*

For side-by-side images use HTML:

<p align="center">
  <img src="docs/images/mfi_b_axis.png" width="45%" alt="MFI along b">
  <img src="docs/images/mfi_a_axis.png" width="45%" alt="MFI along a">
</p>

Tips for GitHub images:
  - Use relative paths from the repo root
  - PNG format, ~800px wide
  - Place a blank line before and after the image line
  - Alt text improves accessibility
-->

Two variants are supported:

| System | Composition (2x2x2 supercell) | Atom types | Atoms |
|--------|-------------------------------|------------|-------|
| Silicalite-1 (all-silica) | Si768 O1536 | Si, O | 2304 |
| H-ZSM-5 (Si/Al ~ 18) | Si728 Al40 O1536 H40 | Si, Si\_a, Si\_b, Al, O, Oas, Ob, Hb | 2344 |

In H-ZSM-5, each Al substitution creates a Bronsted acid site (Si-OH-Al bridge) with distinct atom types: Si\_a (bonded to Oas), Si\_b (bonded to Ob), Oas (bridging Si\_a-Al), Ob (bridging Si\_b-Al-Hb), and Hb (Bronsted proton).

---

## Force Field Description

The force field combines multiple models, each optimized for its interaction type. This is the standard approach in the Sholl group's flexible-framework adsorption studies (Boulfelfel et al. 2016; Daou et al. 2021).

### Framework Intramolecular (Bonded)

The framework uses the **class II (CFF)** functional form from Hill and Sauer (1994, 1995) with the **mHSFF** modification (Boulfelfel et al. 2016) that corrects the T-O-T equilibrium angle to 150 degrees.

Quartic bond stretch:

    E_bond = K2*(b - b0)^2 + K3*(b - b0)^3 + K4*(b - b0)^4

Quartic angle bend:

    E_angle = H2*(theta - theta0)^2 + H3*(theta - theta0)^3 + H4*(theta - theta0)^4

Three-term cosine torsion:

    E_dihedral = V1*[1 - cos(phi - phi1)] + V2*[1 - cos(2*phi - phi2)] + V3*[1 - cos(3*phi - phi3)]

Cross-coupling terms (bond-bond, bond-angle, angle-angle-torsion) are included. These are essential for reproducing the coupling between stretching and bending modes in Si-O-Si linkages.

An improper torsion E\_imp = K\_chi * chi^2 is applied at the Si\_b-Ob-Hb-Al center in H-ZSM-5 to maintain planarity of the Bronsted acid site.

| Interaction | Source |
|---|---|
| Si-O bonds, O-Si-O/Si-O-Si angles | Hill and Sauer, J. Phys. Chem. 1994, 98, 1238 |
| Si-O-Si theta0 = 150 deg (mHSFF) | Boulfelfel et al., J. Phys. Chem. C 2016, 120, 14140 |
| Al-O, acid-site angles/dihedrals | Hill and Sauer, J. Phys. Chem. 1995, 99, 9536 |

### Framework Nonbonded (Host-Host)

Framework atoms interact nonbondedly through a **purely repulsive** inverse-9 potential:

    E_nb(r) = (A_i * A_j) / r^9

where A\_i are atom-specific repulsion parameters (H-S 1994 Table 7, H-S 1995 Table 7):

| Atom | A\_i (kcal^1/2 mol^-1/2 A^9/2) |
|------|------|
| Si | 432.33 |
| O | 239.61 |
| Al | 143.18 |
| Hb | 7.77 |

There is no attractive term — the framework is held together by covalent bonded interactions. The A/r^9 only prevents nonbonded framework atoms (4+ bonds apart) from overlapping. Since LAMMPS has no native pair\_style for this form, these are **tabulated** via `pair_style table spline 2000`.

### Adsorbate (Guest-Guest)

Ethanol uses the **TraPPE-UA** force field (Siepmann et al.):

| Site | sigma (A) | epsilon/kB (K) | epsilon (kcal/mol) | Charge (e) |
|------|-----------|-----------------|---------------------|------------|
| CH3 | 3.750 | 98.0 | 0.19475 | 0.000 |
| CH2 | 3.950 | 46.0 | 0.09141 | +0.265 |
| O\_eth | 3.020 | 93.0 | 0.18481 | -0.700 |
| H\_eth | 0.0 | 0.0 | 0.0 | +0.435 |

Guest-guest nonbonded interactions use standard **LJ 12-6**:

    E_LJ(r) = 4*eps_ij * [(sigma_ij/r)^12 - (sigma_ij/r)^6]

### Host-Guest Cross-Interactions

The A/r^9 was never parametrized for adsorbate interactions and lacks the attractive well needed for adsorption. Host-guest interactions therefore use **separate LJ 12-6 parameters** from Bai et al. (2013), cross-mixed with TraPPE-UA via Lorentz-Berthelot rules:

    sigma_ij = (sigma_i + sigma_j) / 2
    eps_ij   = sqrt(eps_i * eps_j)

Bai framework LJ parameters:

| Atom | sigma (A) | epsilon/kB (K) | Source |
|------|-----------|-----------------|--------|
| O\_fw | 3.30 | 53.0 | Bai et al., J. Phys. Chem. C 2013, 117, 24375 |
| Si\_fw | 2.30 | 22.0 | Bai et al. 2013 |

Al uses Si parameters; Hb has no LJ (Coulomb only).

### Summary: LAMMPS `pair_style hybrid/overlay`

| Sub-style | Interaction | Potential | Source |
|---|---|---|---|
| `table spline 2000` | fw-fw | A/r^9 (repulsive) | Hill-Sauer 1994/1995 |
| `lj/cut 12.0` | fw-eth, eth-eth | LJ 12-6 | Bai 2013 x TraPPE |
| `coul/long 12.0` | all-all | Ewald summation | - |

### Charge Model

Framework charges use the Hill-Sauer bond increment method (Eq. 11 of H-S 1995): q\_i = sum of delta\_ij over bonded neighbors.

| Atom | Charge (e) | Derivation |
|------|-----------|------------|
| Si | +0.5236 | 4 x delta(Si-O) |
| Si\_a | +0.5192 | 3 x delta(Si-O) + delta(Si-Oas) |
| Si\_b | +0.5319 | 3 x delta(Si-O) + delta(Si-Ob) |
| Al | +0.5366 | 3 x delta(Al-Oas) + delta(Al-Ob) |
| O (Oss) | -0.2618 | 2 x delta(O-Si) |
| Oas | -0.2959 | delta(O-Si\_a) + delta(O-Al) |
| Ob | -0.2515 | delta(O-Si\_b) + delta(O-Al) + delta(O-Hb) |
| Hb | +0.0839 | delta(H-Ob) |

Ethanol charges follow TraPPE-UA. Total system is charge-neutral.

---

## Coupled GCMC/NPT-MD Workflow

### Overview

```
           Initial CIF (experimental)
                    |
          +---------v---------+
    +---->|   RASPA3 GCMC     |<-----------+
    |     |   (rigid fw)      |            |
    |     +---------+---------+            |
    |               | N molecules          |
    |     +---------v---------+            |
    |     |  Merge fw + ads   |            |
    |     |  (Julia)          |            |
    |     +---------+---------+            |
    |               | loaded.lmp           |
    |     +---------v---------+            |
    |     |  LAMMPS NPT-MD   |            |
    |     |  (flexible fw)   |            |
    |     +---------+---------+            |
    |               | deformed cell        |
    |     +---------v---------+            |
    |     |  Extract CIF      |------------+
    |     |  (Julia)          | distorted.cif
    |     +---------+---------+
    |               |
    |     +---------v---------+
    +--NO-| Converged?        |
          +---------+---------+
                    | YES
              Final result
```

Each iteration consists of four steps:

1. **GCMC** (RASPA3) — determine equilibrium loading in the current framework.
2. **Build** (Julia) — merge RASPA3 adsorbate positions into the framework.
3. **NPT-MD** (LAMMPS) — relax the loaded framework under constant T, P.
4. **Extract CIF** (Julia) — write distorted framework for the next GCMC cycle.

---

## Simulation Details

### GCMC Stage (RASPA3)

| Parameter | Value |
|---|---|
| Ensemble | Grand canonical (muVT), rigid framework |
| Moves | Translation, rotation, insertion, deletion, reinsertion |
| Adsorbate | TraPPE-UA ethanol (4-site united atom) |
| Electrostatics | Ewald summation (automatic precision) |
| Cutoff | 12 A (LJ and real-space Coulomb) |
| Typical cycles | 5000 init + 5000 equil + 20000 prod |

### NPT-MD Stage (LAMMPS)

The protocol follows Daou et al. (2021):

| Phase | Ensemble | Duration | Purpose |
|-------|----------|----------|---------|
| 0 | Soft potential | ~1 ps | Remove overlaps (skipped if 0 adsorbate) |
| 1 | Minimization | CG, tol 1e-6 | Local minimum |
| 2 | NVT ramp | 10 -> T K, staged | Thermalize adsorbate then framework |
| 3 | NPT equilibration | 200 ps | Cell relaxation |
| 4 | NPT production | 2 ns | Data collection |

Technical parameters:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Timestep | 0.1 fs (inner) | Required for stiff class2 Si-O bonds |
| r-RESPA | 2:1 outer:inner | Bonded every step, nonbonded every 2 |
| Thermostat | Nose-Hoover, tchain=6 | Daou et al. protocol |
| Barostat | Anisotropic, pchain=6 | Independent a, b, c relaxation |
| tdamp | 10 fs | Thermostat relaxation |
| pdamp | 1000 fs | Barostat relaxation |
| Cutoff | 12 A | Both vdW and Coulomb |
| Ewald precision | 1e-6 | Long-range electrostatics |
| special\_bonds | 0.0 0.0 1.0 | Full nonbonded for 1-4+ pairs |

The barostat is `aniso` (not `iso`) because MFI is orthorhombic — channels run in different directions, so the mechanical response is anisotropic.

---

## Convergence and Statistics

### Convergence Criteria

The workflow monitors two quantities over a sliding window of the last W iterations (default W = 5):

- **Loading**: coefficient of variation CV(N\_ads) over the window.
- **Volume**: coefficient of variation CV(V) over the window.

Convergence is declared when both CV values fall below a tolerance (default 0.5%). Maximum iterations are capped at `max_iterations` (default 10-15).

In practice, volume converges within 3-5 cycles. Loading converges more slowly due to GCMC stochasticity.

### Statistical Validity

- **GCMC loading fluctuations**: With `raspa_n_prod = 5000`, individual cycles show +/-20-30% variance. This is inherent to the grand canonical ensemble at finite sampling. Increasing to 50000+ reduces CV to ~3-5%.
- **NPT cell fluctuations**: Volume varies ~0.1-0.3% between cycles — expected thermodynamic fluctuation for a 2x2x2 supercell.
- **Reporting**: Average over the last W converged iterations. Error bars = standard deviation across these iterations.
- **Zero-loading regime**: At very low pressures, RASPA3 may return 0 molecules. The code handles this (empty framework MD), but these points have large relative uncertainty — multiple independent runs are recommended.

---

## Code Architecture

```
JZeoMCMD/src/
+-- JZeoMCMD.jl                    # Module definition, includes, exports
+-- read_lammps_data.jl            # LAMMPS data file reader/writer
+-- add_zeolite_topology.jl        # Distance-based bond/angle/dihedral detection
+-- refine_topology_types.jl       # Refine bonded types for aluminosilicate
+-- alumino_support.jl             # Aluminosilicate configs and helpers
+-- build_loaded_zeolite.jl        # Merge framework + adsorbate into LAMMPS data
+-- reload_adsorbate.jl            # Strip old adsorbate, insert new from RASPA3
+-- write_cif.jl                   # Write CIF from LAMMPS data
+-- workflow.jl                    # Master GCMC/NPT-MD iterative loop
+-- setup_pressure_sweep.jl        # Directories for isotherm pressure points
+-- analysis.jl                    # Log parsing, channel occupancy, strain
+-- generate_nb_tables.jl          # Generate A/r^9 pair tables
+-- params_loader.jl               # TOML configuration loader
+-- ff/
    +-- hillsauer_silica_loaded.ff       # All-silica + ethanol
    +-- hillsauer_alumsil_empty.ff       # H-ZSM-5 framework only
    +-- hillsauer_alumsil_loaded.ff      # H-ZSM-5 + ethanol
```

Key design decisions:

- **Topology-only data files**: `.lmp` files contain atoms and topology. All coefficients live in `.ff` include files. This decouples structure from parameters.
- **Zero-adsorbate safety**: When 0 molecules are loaded, all adsorbate types and bonded type counts are still registered in the data file header, ensuring `.ff` coefficients have valid indices. LAMMPS scripts use `if "${nads} == 0"` jumps to skip adsorbate phases.
- **Aluminosilicate type refinement**: Topology builder collapses Si variants and O variants for distance-based detection, then restores original 8 types and calls `refine_topology_types!` for correct Hill-Sauer bonded indices.

---

## Installation

```bash
git clone https://github.com/<your-username>/JZeoMCMD.jl.git
julia -e 'using Pkg; Pkg.develop(path="JZeoMCMD.jl")'
```

### Dependencies

- **Julia** >= 1.9
- **LAMMPS** with CLASS2, KSPACE, MOLECULE, EXTRA-PAIR packages
- **RASPA3** compiled with ethanol TraPPE-UA support
- Julia packages: Printf, LinearAlgebra, Statistics, JSON, TOML, Dates

---

## Input Files

### All-silica MFI

| File | Description | Source |
|------|-------------|--------|
| `MFI_SI.data` | LAMMPS data (OVITO export) | Export from OVITO |
| `MFI_SI.cif` | Initial CIF | IZA database |
| `hillsauer_nb.table` | A/r^9 tables (3 pairs) | `julia generate_nb_tables.jl` |
| `hillsauer_silica_loaded.ff` | Force field | Provided in repo |
| `run_npt.in` | LAMMPS input | Provided in repo |
| `raspa_inputs/` | RASPA3 templates | Provided in repo |

### Aluminosilicate H-ZSM-5

| File | Description | Source |
|------|-------------|--------|
| `corrected.data` | 8-type LAMMPS data (OVITO) | Export with Al/Si\_a/Si\_b types |
| `corrected.cif` | CIF with pseudo-atom names | Must match force\_field.json |
| `corrected_with_topology.data` | Refined topology | `julia run_update_data.jl` |
| `hillsauer_alumsil_nb.table` | A/r^9 tables (36 pairs) | `julia generate_nb_tables.jl --rmin 0.5` |
| `hillsauer_alumsil_loaded.ff` | Force field | Provided in repo |
| `run_npt_alumsil_loaded.in` | LAMMPS input | Provided in repo |

**Important**: CIF `_atom_site_type_symbol` must use RASPA3 pseudo-atom names (Si, Si\_a, Si\_b, Al, O, Oas, Ob, Hb), not element symbols. These must match `"name"` fields in `force_field.json`.

---

## Usage Examples

### Running the Workflow

**All-silica:**

```julia
using JZeoMCMD

wp = WorkflowParams(
    raspa_n_init      = 5000,
    raspa_n_equil     = 5000,
    raspa_n_prod      = 20000,
    lammps_npt_steps  = 2000000,
    lammps_timestep   = 0.25,
    temperature       = 303.0,
    pressure          = 1e5,
    max_iterations    = 15,
)
wp.initial_cif  = "MFI_SI.cif"
wp.initial_data = "MFI_SI.data"
wp.lammps_exe   = "mpirun -np 4 lmp_mpi"
wp.lammps_input = "run_npt.in"
wp.ff_include   = "hillsauer_silica_loaded.ff"

run_gcmc_md_workflow(wp)
```

**Aluminosilicate:**

```julia
using JZeoMCMD

wp = WorkflowParams(
    raspa_n_init      = 5000,
    raspa_n_equil     = 5000,
    raspa_n_prod      = 20000,
    lammps_npt_steps  = 2000000,
    lammps_timestep   = 0.25,
    temperature       = 373.0,
    pressure          = 1e5,
    is_alumino        = true,
    nfw_atoms         = 2344,
    table_file        = "hillsauer_alumsil_nb.table",
    max_iterations    = 15,
)
wp.initial_cif  = "corrected.cif"
wp.initial_data = "corrected_with_topology.data"
wp.lammps_exe   = "mpirun -np 4 lmp_mpi"
wp.lammps_input = "run_npt_alumsil_loaded.in"
wp.ff_include   = "hillsauer_alumsil_loaded.ff"

run_gcmc_md_workflow(wp)
```

### Standalone Scripts

**Build aluminosilicate topology** (one-time):

```bash
julia run_update_data.jl
# Output: corrected_with_topology.data
# Expected: 3112 bonds (6 types), 6224 angles (10 types),
#           9456 dihedrals (10 types), 40 impropers (1 type)
```

**Generate nonbonded tables:**

```bash
# All-silica
julia generate_nb_tables.jl

# Aluminosilicate (rmin=0.5 for O-H bonds at 0.954 A)
julia generate_nb_tables.jl --rmin 0.5
```

**Empty framework validation** (no ethanol):

```bash
julia run_update_data.jl
mpirun -np 4 lmp_mpi -in run_npt_alumsil_empty.in
awk '{print $1, $2, $3, $7}' cell_params.dat | tail -10
# Expected: a ~ 41.1, b ~ 40.7, c ~ 27.3 (H-S prediction, ~2.5% > experiment)
```

**Pressure sweep:**

```julia
using JZeoMCMD

pressures = [100, 1000, 5000, 10000, 50000, 100000]  # Pa
setup_pressure_sweep("isotherm_303K", pressures;
    template_dir=".", temperature=303.0)
```

---

## Output Structure

```
cycle_01/
+-- raspa/                    # RASPA3 GCMC results
+-- loaded.lmp                # Framework + adsorbate data file
+-- lammps/
|   +-- loaded.lmp            # Copy for LAMMPS
|   +-- hillsauer_*.ff        # Force field include
|   +-- *_nb.table            # Nonbonded tables
|   +-- run_npt*.in           # LAMMPS input
|   +-- cell_params.dat       # Cell parameters per timestep
|   +-- traj.lammpstrj        # Trajectory
|   +-- loaded_npt_final.lmp  # Final structure
+-- distorted.cif             # CIF for next cycle
convergence.csv               # Running convergence record
```

---

## References

1. Hill, J.-R.; Sauer, J. *J. Phys. Chem.* **1994**, 98, 1238-1244. (Silica force field)
2. Hill, J.-R.; Sauer, J. *J. Phys. Chem.* **1995**, 99, 9536-9550. (Aluminosilicate extension)
3. Boulfelfel, S. E. et al. *J. Phys. Chem. C* **2016**, 120, 14140-14148. (mHSFF)
4. Daou, A. S. S. et al. *J. Phys. Chem. C* **2021**, 125, 5296-5305. (Flexible framework protocol)
5. Emelianova, A. et al. **2023**. (AID in zeolites 4A and 13X)
6. Bai, P. et al. *J. Phys. Chem. C* **2013**, 117, 24375-24387. (TraPPE-zeo)
7. Martin, M. G.; Siepmann, J. I. *J. Phys. Chem. B* **1998**, 102, 2569-2577. (TraPPE-UA)
8. Fang, H. et al. *J. Phys. Chem. C* **2018**, 122, 12880-12891. (CH4 in zeolites)

---

## License

<!-- MIT / Apache-2.0 / GPL-3.0 — choose yours -->

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

*README last updated: April 2026*