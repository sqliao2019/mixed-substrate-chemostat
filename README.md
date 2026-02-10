MATLAB Code for Mixed-Substrate Fermentation Manuscript
======================================================

This folder contains MATLAB scripts and precomputed data used to generate the
figures in the manuscript:

"Establishing Quantitative Guidelines for Microbial Mixed Substrate Fermentation."

Scripts are organized by figure number:

  `FigX_*_run.m`   : (optional) recomputes data for Figure X (may be time-consuming)
  
  `FigX_*_plot.m`  : generates Figure X using saved outputs in the data/ folder

Precomputed outputs are provided in:

  data/data_FigX.mat

To reproduce a figure, run the corresponding plotting script in MATLAB, e.g.:

  Fig1_2_plot
  
  Fig4_1_plot

Requirements: MATLAB + Symbolic Math Toolbox (for scripts using solve).
