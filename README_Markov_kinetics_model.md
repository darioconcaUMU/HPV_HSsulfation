# README — Markov kinetics model

## Overview
This folder contains MATLAB scripts and functions for simulating multibond detachment, converting lifetimes into survival curves, estimating effective dissociation behavior, and evaluating Poisson occupancy or coverage effects.

## Contents

### `Calculate_koff_multivalent.m`
Main function for generating multivalent dissociation behavior and fitting an effective single-exponential model with offset.

**How to call**
```matlab
[k, A, C, data] = Calculate_koff_multivalent(...)
```

### `Calculate_koff_multivalent_with_error.m`
Script for propagating uncertainty in dissociation times by sampling lower and upper bounds and recalculating effective dissociation behavior.

**How to call**
```matlab
Calculate_koff_multivalent_with_error
```

### `Parameters_for_simulation.m`
Parameter setup script for multibond lifetime simulations. Defines parameter sets and runs example simulations.

**How to call**
```matlab
Parameters_for_simulation
```

### `circle_coverage.m`
Script for Poisson occupancy and geometric coverage calculations in circular interaction areas.

**How to call**
```matlab
circle_coverage
```

### `simulate_poisson_multibond_lifetimes.m`
Function for simulating detachment lifetimes when the number of bonds follows a Poisson distribution, excluding zero-bond particles and allowing bond-number-dependent off-rates.

**How to call**
```matlab
lifetimes = simulate_poisson_multibond_lifetimes(...)
```

### `survival_from_lifetimes.m`
Function that converts a vector of lifetimes into a survival curve with optional censoring at a maximum observation time.

**How to call**
```matlab
[t, S] = survival_from_lifetimes(lifetimes, ...)
```

## Minimal example workflow
```matlab
Parameters_for_simulation
lifetimes = simulate_poisson_multibond_lifetimes(...);
[t, S] = survival_from_lifetimes(lifetimes, ...);
[k, A, C, data] = Calculate_koff_multivalent(...);
```
