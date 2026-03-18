# README --- Simulation SNR Scripts

## Overview

This repository contains MATLAB scripts for simulating signal-to-noise
ratio (SNR) conditions and evaluating their impact on detection,
tracking, and analysis. The scripts generate synthetic data, apply noise
models, and quantify resulting signal properties.

------------------------------------------------------------------------

## Main Scripts

### `simulate_SNR.m`

Primary simulation script.

**What it does:** - Generates synthetic signal data - Adds noise with
controlled SNR levels - Produces simulated datasets for analysis

**How to call:**

``` matlab
simulate_SNR
```

------------------------------------------------------------------------

### `analyze_SNR.m`

Analysis of simulated SNR datasets.

**What it does:** - Computes signal statistics - Evaluates detection
performance as a function of SNR - Outputs summary metrics

**How to call:**

``` matlab
analyze_SNR
```

------------------------------------------------------------------------

## Supporting Functions

### `generate_signal.m`

Creates synthetic signal traces.

``` matlab
signal = generate_signal(params)
```

### `add_noise.m`

Applies noise to signals.

``` matlab
noisy_signal = add_noise(signal, snr)
```

### `compute_SNR_metrics.m`

Computes SNR-related metrics.

``` matlab
metrics = compute_SNR_metrics(signal)
```

### `plot_SNR_results.m`

Visualizes simulation outputs.

``` matlab
plot_SNR_results(results)
```

------------------------------------------------------------------------

## Minimal Example Workflow

``` matlab
% Generate simulated data
simulate_SNR

% Analyze results
analyze_SNR
```

------------------------------------------------------------------------

## Notes

-   Scripts are designed to run independently with internally defined
    parameters.
-   Outputs are typically stored in workspace variables or saved to MAT
    files.
-   All required user-defined functions are included in this repository.
