# README --- TIRF / Tracking Analysis Scripts

## Overview

This repository contains MATLAB scripts and apps for analyzing particle
trajectories, intensity signals, and surface interactions, primarily
from imaging data (e.g., TIRF). The workflow is organized around main
entry-point scripts, analysis functions, and supporting utilities.

------------------------------------------------------------------------

## Main Entry Points (Most Recent)

### 1. `Main_v3_1.m`

Primary script to run the full analysis pipeline.

**What it does:** - Loads input data - Runs tracking and preprocessing -
Applies filtering and corrections - Outputs processed trajectories and
intermediate results

**How to call:**

``` matlab
Main_v3_1
```

------------------------------------------------------------------------

### 2. `Main_quantify_4v.m`

Latest version of the quantification pipeline.

**What it does:** - Processes tracked data - Computes quantitative
metrics - Aggregates results across datasets

**How to call:**

``` matlab
Main_quantify_4v
```

------------------------------------------------------------------------

### 3. `tirf_analyzer_3v1.m`

Analysis script focused on TIRF datasets.

**What it does:** - Handles image-based analysis - Extracts intensity
and spatial information - Interfaces with downstream quantification
functions

**How to call:**

``` matlab
tirf_analyzer_3v1
```

------------------------------------------------------------------------

### 4. `Main_UI_v3_1.mlapp`

Graphical user interface.

**What it does:** - Provides UI for parameter selection - Allows
interactive execution

**How to run:**

``` matlab
Main_UI_v3_1
```

------------------------------------------------------------------------

### 5. `select_threshold.mlapp`

Interactive threshold selection tool.

**What it does:** - Allows manual threshold selection

**How to run:**

``` matlab
select_threshold
```

------------------------------------------------------------------------

## Core Functional Modules

### Pipeline / Control

-   `func_main_v3_1.m`
-   `func_reanalysis_v3_0.m`
-   `repeat_analysis.m`

### Tracking & Trajectory Processing

-   `track_analysis.m`
-   `separate_tracks.m`
-   `remove_duplicates_2.m`
-   `filter_duration.m`
-   `filter_duration_short.m`

### Intensity & Signal Analysis

-   `analysis_fract_intensity.m`
-   `intensity_histogram.m`
-   `add_fraction_analysis.m`

### Cumulative / Kinetic Analysis

-   `cumulative_analysis.m`
-   `cumulative_analysis_batch_2v0.m`
-   `cumulative_analysis_batch_quantify.m`
-   `extract_fit.m`
-   `extract_fit_cumulative.m`
-   `extract_IR_cumulative.m`
-   `doubleExpFit.m`

### Image Processing / Detection

-   `find_peaks_2D.m`
-   `get_size_image.m`
-   `threshold_selection.m`
-   `threshold_sweep.m`

### Surface / Spatial Analysis

-   `surface_coverage.m`
-   `surface_coverage_stuck.m`

### Correction & Preprocessing

-   `bleach_correction_step.m`
-   `blink_correction_4.m`

### Output / Visualization

-   `create_final_videoNEW.m`
-   `create_final_video.m`

### Utility

-   `calculate_properties.m`

------------------------------------------------------------------------

## Minimal Example Workflow

``` matlab
Main_v3_1
Main_quantify_4v
```

------------------------------------------------------------------------

## Legacy / Older Scripts (Footnote)

-   `Main_quantify.m`
-   `Main_quantify_3v1.m`
-   `create_final_video.m`
-   `cumulative_analysis_batch_2v0.m`
-   `func_reanalysis_v3_0.m`
-   `extract_fit.asv`
-   `playground.m`
