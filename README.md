MATLAB Code and Data for Mixed-Substrate Fermentation Manuscript
======================================================

This repository contains the MATLAB source code and precomputed data used to generate the results and figures in the manuscript:

**“Establishing Quantitative Guidelines for Microbial Mixed Substrate Fermentation”**

## Requirements

The code was tested using `MATLAB R2024b` on `Windows 11`.

The `Symbolic Math Toolbox` is required for scripts that use `solve`. No nonstandard hardware is required.

No installation is needed beyond downloading the repository and opening it in MATLAB. Setup should take less than five minutes.

## Files

Scripts are organized by figure number:

- `FigX_*_run.m` recomputes the data used for Figure X. These calculations are optional when using the supplied precomputed data and may take longer to run.
- `FigX_*_plot.m` generates Figure X using the corresponding saved data.

Precomputed outputs are provided in the `data` folder, for example:

```text
data/data_FigX.mat
```

## Reproducing the figures

Set the repository folder as the current folder in MATLAB and run the corresponding plotting script. For example:

```matlab
Fig1_2_plot
```

or:

```matlab
Fig4_1_plot
```

The plotting script loads the supplied data and generates the corresponding manuscript figure. Most plotting scripts should finish within a few seconds.

To recompute the underlying data, first run the corresponding `FigX_*_run.m` script and then run the plotting script. Computation time varies by analysis and may range from a few seconds to several hours.

## Using the code with modified parameters

Model parameters and analysis settings can be changed in the corresponding `FigX_*_run.m` scripts. Use a different output filename when testing modified settings to avoid overwriting the supplied manuscript data.
