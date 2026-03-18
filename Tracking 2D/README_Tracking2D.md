# README --- 2D Tracking Analysis Scripts

## Overview

This repository contains MATLAB functions for processing, analyzing, and
visualizing 2D particle tracking data. The workflow includes trajectory
extraction, cleaning, segmentation, diffusion analysis, and
visualization.

------------------------------------------------------------------------

## Main Functions

### `analysis_full_tracks.m`

Main analysis routine for full trajectories.

**What it does:** - Processes full particle tracks - Computes
trajectory-based metrics - Calls internal helper function `analyseTrack`

**How to call:**

``` matlab
analysis_full_tracks(data)
```

------------------------------------------------------------------------

### `compile_subtracks.m`

Aggregates and summarizes segmented trajectories.

**What it does:** - Combines subtracks into a unified dataset - Computes
summary statistics - Uses internal helper `summerise_data`

**How to call:**

``` matlab
compile_subtracks(subtracks)
```

------------------------------------------------------------------------

### `classify_confined_free_by_fit.m`

Classifies motion type based on fitting.

**What it does:** - Fits models to trajectories - Classifies tracks into
confined or free diffusion

**How to call:**

``` matlab
classify_confined_free_by_fit(tracks)
```

------------------------------------------------------------------------

### `histogram_D_from_MSD.m`

Generates diffusion coefficient distributions.

**What it does:** - Computes diffusion coefficients from MSD - Produces
histograms of diffusion values

**How to call:**

``` matlab
histogram_D_from_MSD(msd_data)
```

------------------------------------------------------------------------

## Trajectory Processing

### `extract_track_position.m`

Extracts position data from track structures.

``` matlab
positions = extract_track_position(track)
```

### `remove_double_tracks.m`

Removes duplicate or overlapping tracks.

``` matlab
clean_tracks = remove_double_tracks(tracks)
```

### `clean_outliers.m`

Filters out anomalous trajectory points.

``` matlab
clean_data = clean_outliers(data)
```

### `compile_subtracks.m`

Segments and compiles trajectory subsets.

``` matlab
compile_subtracks(subtracks)
```

------------------------------------------------------------------------

## Motion and Diffusion Analysis

### `plot_MSD_mean.m`

Plots mean squared displacement.

``` matlab
plot_MSD_mean(msd_data)
```

### `plot_avg_track.m`

Plots average trajectory behavior.

``` matlab
plot_avg_track(tracks)
```

### `motion_correction_and_visualization.m`

Applies motion correction and visualizes results.

``` matlab
motion_correction_and_visualization(data)
```

------------------------------------------------------------------------

## Geometry / Surface Interaction

### `calculate_hep_stretch.m`

Computes stretch or extension metrics.

``` matlab
calculate_hep_stretch(data)
```

### `numerical_hep_stertch.m`

Numerical implementation of stretch calculation.

``` matlab
numerical_hep_stertch(data)
```

------------------------------------------------------------------------

## Input / Data Handling

### `unpack_TrackMate_XML.m`

Reads and converts TrackMate XML files.

``` matlab
tracks = unpack_TrackMate_XML(filename)
```

------------------------------------------------------------------------

## Minimal Example Workflow

``` matlab
% Load tracking data
tracks = unpack_TrackMate_XML('file.xml');

% Clean data
tracks = remove_double_tracks(tracks);
tracks = clean_outliers(tracks);

% Analyze trajectories
analysis_full_tracks(tracks);

% Compute diffusion
histogram_D_from_MSD(tracks);
```

------------------------------------------------------------------------

## Dependencies

The scripts rely on standard MATLAB functions and may require:

-   Statistics and Machine Learning Toolbox:
    -   `nlinfit`
    -   `kmeans`
    -   `pdist`
-   MATLAB XML support:
    -   `xmlread`

------------------------------------------------------------------------

## Notes

-   Helper functions such as `analyseTrack` and `summerise_data` are
    defined within their respective parent files.
-   All required user-defined functions are included in this repository.
