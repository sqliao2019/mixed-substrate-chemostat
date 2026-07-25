MATLAB Code and Data for Mixed-Substrate Fermentation Manuscript
======================================================

This repository contains the MATLAB source code and precomputed data used to generate the results and figures in the manuscript:

**“Establishing Quantitative Guidelines for Microbial Mixed Substrate Fermentation”**

## Requirements

No installation is needed beyond downloading the repository and opening it in MATLAB. Setup should take less than five minutes.

The `Symbolic Math Toolbox` is required for scripts that use `solve`. No nonstandard hardware is required.

The code was tested using `MATLAB R2024b` on `Windows 11`.

## Repository contents

Scripts are organized by figure number:

- `FigX_*_run.m` recomputes the data used for Figure X.
- `FigX_*_plot.m` generates Figure X using the corresponding saved data.

Precomputed outputs are provided in the `data` folder, for example:

```text
data/data_FigX.mat
```

## Running the example workflows

The supplied code and data reproduce the representative cases shown in the manuscript and serve as examples for running the analysis.

Set the repository folder as the current folder in MATLAB and run the corresponding plotting script. For example:

```matlab
Fig1_2_plot
```

or:

```matlab
Fig4_1_plot
```

Each plotting script loads the corresponding supplied data, when required, and generates the corresponding manuscript figure. The expected output is the figure associated with that script, and most plotting scripts should finish within a few seconds.

## Recomputing data and modifying parameters

To recompute the underlying data, first run the corresponding `FigX_*_run.m` script and then run the plotting script. Computation time varies by analysis and may range from a few seconds to several hours.

Users can explore other conditions by modifying the model parameters and analysis settings in the corresponding `FigX_*_run.m` scripts. When generating new results, use a different output filename to avoid overwriting the supplied manuscript data.

