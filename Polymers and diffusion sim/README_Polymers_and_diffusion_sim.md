# README — Polymers and diffusion sim

## Overview
This folder contains MATLAB scripts and functions for 2D multivalent diffusion simulations, 3D multi-tether Brownian simulations, visualization of simulated trajectories, and downstream trajectory analysis.

## Contents

### `multivalent_diffusion_sim.m`
Core 2D multivalent diffusion simulator. The particle moves on a 2D binder landscape by forming and breaking bonds to nearby binders.

**How to call**
```matlab
[track, binders_loc] = multivalent_diffusion_sim(binders, kinetics, time)
```

### `multivalent_diffusion_sim_fix_r.m`
Variant of the 2D multivalent diffusion simulation with a fixed interaction-radius implementation.

**How to call**
```matlab
[track, binders_loc] = multivalent_diffusion_sim_fix_r(binders, kinetics, time)
```

### `brownianTether3DMulti.m`
Core 3D Brownian dynamics simulation of a sphere tethered to a surface through multiple anchors.

**How to call**
```matlab
[tArray, posArray, attachmentArray, anchorPos] = brownianTether3DMulti(N_t)
```

### `showParticleVideo3DMulti.m`
Visualization utility for generating an animation of a multi-tether 3D simulation.

**How to call**
```matlab
showParticleVideo3DMulti(tArray, posArray, attachmentArray, anchorPos)
```

### `analyze3DDistribution.m`
Function that converts 3D positions into a 2D XY histogram for visualization.

**How to call**
```matlab
out = analyze3DDistribution(posArray, ...)
```

### `analyse_simulation_diffusion.m`
Post-processing script for 2D diffusion simulations. Extracts trajectory-level quantities from simulated tracks.

**How to call**
```matlab
analyse_simulation_diffusion
```

### `exploring_landscape_diffusion.m`
Batch simulation script that sweeps kinetic and binder-concentration parameters and stores repeated 2D diffusion simulations.

**How to call**
```matlab
exploring_landscape_diffusion
```

### `generate_statistics_3D_binder.m`
Batch script for running `brownianTether3DMulti` over a range of binder counts and extracting radial XY displacements.

**How to call**
```matlab
generate_statistics_3D_binder
```

### `func_MSS_Tracks_simulation.m`
Wrapper function that converts simulated tracks into the format expected by the transient-diffusion / MSS analysis pipeline.

**How to call**
```matlab
out = func_MSS_Tracks_simulation(track, ...)
```

### `basicTransientDiffusionAnalysisv1.m`
Transient-diffusion / MSS analysis function for classifying motion segments in trajectories.

**How to call**
```matlab
results = basicTransientDiffusionAnalysisv1(data, ...)
```

### `calculate_average_motion.m`
Utility function that computes an average motion class from transient-diffusion analysis output.

**How to call**
```matlab
avg_motion = calculate_average_motion(results)
```

### `plot_selected_trajectories.m`
Plotting utility for visualizing selected trajectories from a simulation output array.

**How to call**
```matlab
plot_selected_trajectories(track, ...)
```

## Minimal example workflows

### 2D diffusion workflow
```matlab
[track, binders_loc] = multivalent_diffusion_sim(binders, kinetics, time);
analyse_simulation_diffusion
plot_selected_trajectories(track)
```

### 3D tether workflow
```matlab
[tArray, posArray, attachmentArray, anchorPos] = brownianTether3DMulti(N_t);
showParticleVideo3DMulti(tArray, posArray, attachmentArray, anchorPos)
```
