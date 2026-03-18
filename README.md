# GitHub scripts — README

## Overview

This repository contains MATLAB code organized into two main themes:

1. **Markov kinetics model**
   Tools for simulating multibond detachment, converting lifetimes into survival curves, estimating effective dissociation behavior, and evaluating Poisson occupancy / geometric coverage.

2. **Polymers and diffusion sim**
   Tools for 2D multivalent diffusion simulations on a binder landscape, 3D multi-tether Brownian simulations, trajectory visualization, and downstream motion analysis.

The codebase includes both **standalone functions** and **workspace-driven scripts**. The functions are reusable building blocks, while the scripts are typically intended for parameter exploration, batch simulation, or post-processing.

---

## Folder structure

### `Markov kinetics model/`

#### `Calculate_koff_multivalent.m`
Main fitting function for constructing multivalent dissociation curves and fitting them with an effective single-exponential model with offset.

**Outputs:**
- `k` — fitted effective off-rate
- `A` — fitted amplitude
- `C` — fitted offset
- `data` — generated time-series data used for the fit

**Local helper:**
- `fit_single_exp_offset`

#### `Calculate_koff_multivalent_with_error.m`
Script that propagates uncertainty in dissociation times by sampling lower and upper bounds and recomputing effective behavior. Useful for sensitivity analysis and uncertainty visualization.

#### `circle_coverage.m`
Script for Poisson occupancy and coverage calculations in circular interaction areas. It estimates probabilities such as the number of binders in a capture region and can generate simulated distance/coverage maps.

#### `Parameters_for_simulation.m`
Example setup script defining parameter sets for multibond lifetime simulations. It runs lifetime simulations for several parameter choices and converts them to survival curves.

#### `simulate_poisson_multibond_lifetimes.m`
Core function for simulating particle detachment lifetimes when:
- the number of bonds follows a Poisson distribution,
- zero-bond particles are excluded,
- bond count is capped at 5,
- the off-rate depends on bond number.

The output lifetimes are quantized to a sampling interval `dt`.

#### `survival_from_lifetimes.m`
Converts a vector of lifetimes into a Kaplan–Meier-style survival curve, with censoring above a maximum observation time.

---

### `Polymers and diffusion sim/`

#### `analyse_simulation_diffusion.m`
Post-processing script for 2D diffusion simulations. It assumes a precomputed `track` cell array and extracts quantities such as:
- total path length,
- maximum displacement,
- diffusion coefficient from MSD,
- average number of bonds.

#### `analyze3DDistribution.m`
Utility function that converts 3D positions into a 2D histogram of the XY distribution for visualization.

#### `basicTransientDiffusionAnalysisv1.m`
Large transient-diffusion / MSS analysis function for classifying motion segments in trajectories. It accepts tracks in matrix or structure form and returns segment-level diffusion classification metrics.

This is a self-contained analysis routine with several embedded helper functions for segmentation, classification, localization precision estimation, confinement analysis, and asymmetric diffusion analysis.

#### `brownianTether3DMulti.m`
Core 3D Brownian dynamics simulation of a sphere tethered to a surface through multiple anchors.

Features include:
- overdamped translation and rotation,
- multiple tethers,
- steric wall constraint,
- nonlinear tether-force models,
- optional plotting of simulated trajectories.

**Outputs:**
- `tArray` — time vector
- `posArray` — sphere-center positions
- `attachmentArray` — attachment points on the particle over time
- `anchorPos` — fixed anchor positions on the surface

**Local helpers:**
- `projectuv`
- `rodrigues`
- `sphere2cart`
- `forceCalc`

#### `calculate_average_motion.m`
Utility function that computes a weighted average motion class from the output of transient diffusion / MSS analysis.

#### `exploring_landscape_diffusion.m`
Batch simulation script for sweeping:
- `k_on`
- binder concentration

For each parameter combination it runs repeated calls to the 2D diffusion simulator and stores the outputs.

#### `func_MSS_Tracks_simulation.m`
Wrapper function that converts the simulated tracks into the format expected by the MSS pipeline and runs transient diffusion analysis.

**Local helper:**
- `subfunc_MSS_convertData`

#### `generate_statistics_3D_binder.m`
Batch script for running `brownianTether3DMulti` over a range of binder counts. It then centers each trajectory and computes radial XY displacements for statistical analysis.

#### `multivalent_diffusion_sim.m`
Core 2D multivalent binder-field simulation.

The particle moves on a 2D plane by forming and breaking bonds to nearby binders. Its position is determined by the center of mass of currently attached bonds, and the simulation ends when no bonds remain.

**Outputs:**
- `track` — simulation results, including positions and bond history
- `binders_loc` — binder locations used in the simulation

#### `multivalent_diffusion_sim_fix_r.m`
Variant of the 2D multivalent diffusion simulation in which the interaction radius is treated differently / kept fixed according to the implementation in the file. This is best viewed as an alternative simulator for comparison against `multivalent_diffusion_sim.m`.

#### `plot_selected_trajectories.m`
Plotting utility for visualizing a selected subset of trajectories from a simulation output array.

**Local helper:**
- `plot_subplot`

#### `showParticleVideo3DMulti.m`
Visualization utility for generating a 3D animation of a multi-tether simulation, including sphere projections onto coordinate planes and tether lines.

---

## Typical workflows

### 1. Markov / multibond lifetime workflow
A typical sequence is:

1. Define `mu` and bond-specific `koff`
2. Run `simulate_poisson_multibond_lifetimes`
3. Convert lifetimes with `survival_from_lifetimes`
4. Explore effective behavior with `Calculate_koff_multivalent`
5. Optionally assess uncertainty with `Calculate_koff_multivalent_with_error`

### 2. 2D multivalent diffusion workflow
A typical sequence is:

1. Define `binders`, `kinetics`, and `time`
2. Run `multivalent_diffusion_sim` or `multivalent_diffusion_sim_fix_r`
3. Explore parameter sweeps with `exploring_landscape_diffusion`
4. Analyze output using `analyse_simulation_diffusion`
5. Visualize trajectories with `plot_selected_trajectories`
6. Convert tracks for MSS analysis using `func_MSS_Tracks_simulation`
7. Summarize motion classes with `calculate_average_motion`

### 3. 3D tether / polymer workflow
A typical sequence is:

1. Run `brownianTether3DMulti`
2. Batch-run over binder counts using `generate_statistics_3D_binder`
3. Inspect XY distributions with `analyze3DDistribution`
4. Create trajectory movies with `showParticleVideo3DMulti`

---

## Scripts vs functions

The repository contains a mix of:

### Functions
Reusable code with explicit inputs and outputs, for example:
- `simulate_poisson_multibond_lifetimes`
- `survival_from_lifetimes`
- `Calculate_koff_multivalent`
- `multivalent_diffusion_sim`
- `multivalent_diffusion_sim_fix_r`
- `brownianTether3DMulti`
- `basicTransientDiffusionAnalysisv1`

### Scripts
Files that operate on variables in the current workspace and are generally intended for interactive use or batch runs, for example:
- `Calculate_koff_multivalent_with_error`
- `circle_coverage`
- `Parameters_for_simulation`
- `analyse_simulation_diffusion`
- `exploring_landscape_diffusion`
- `generate_statistics_3D_binder`

For portability and reproducibility, the functions are easier to reuse programmatically, while the scripts are convenient for exploratory work.

---

## MATLAB requirements

Several files appear to rely on core MATLAB plus common toolboxes.

### Likely required toolbox functionality
- **Statistics and Machine Learning Toolbox**
  - `poisspdf`
  - `poisscdf`
  - `poissrnd`
  - `histcounts2`
  - `pdist`

- **Optimization Toolbox**
  - `lsqcurvefit`
  - `optimoptions`

Some visualization functions also rely on standard MATLAB graphics and video-writing functionality.

---

## Notes for use

- Units are not identical across all files. In particular, some diffusion simulations use normalized units, while the 3D tether simulation is written in SI-scale units.
- Several scripts assume variables such as `track`, `time`, `binders`, or `kinetics` already exist in the workspace.
- Some files are exploratory or analysis-oriented rather than packaged as general-purpose functions.
- The repository contains both simulation code and analysis code, so it is useful to keep output `.mat` files and intermediate variables organized per workflow.

---

## Short summary

This repository provides a compact MATLAB toolkit for:
- multibond lifetime simulation,
- survival and effective off-rate analysis,
- 2D multivalent diffusion on heterogeneous binder fields,
- 3D multi-tether Brownian dynamics,
- visualization and MSS-based motion classification.

It is best used as a research workflow repository rather than as a single packaged MATLAB toolbox.
