# Spatial-Temporal GPR

MATLAB code for spatial-temporal Gaussian process regression, including marginal-likelihood, generalized cross-validation (GCV), and Stein's unbiased risk estimate (SURE) hyperparameter selection.

## Usage

Start MATLAB in this directory. Review the dataset selection in Load.m and the settings in the driver before running:

- Main4GPR.m: compares marginal-likelihood, GCV, and SURE tuning.
- MainGPR_ML.m: marginal-likelihood tuning.
- Load.m: loads or constructs the selected dataset; currently selects weather data from TAVG20.mat.

MATLAB and Optimization Toolbox (fmincon) are required. Additional requirements may depend on the selected kernel and plotting functions. The supplied MAT files include datasets and saved experimental quantities. Numerical experiments have not been rerun for this publication, and no minimum MATLAB release is asserted.

## Citation and attribution

If you use this code in research, please acknowledge this repository and cite the underlying methods you use. The original implementation credits these papers:

1. A. Carron, M. Todescato, R. Carli, L. Schenato, and G. Pillonetto, "Machine Learning meets Kalman Filtering," 55th Conference on Decision and Control, 2016, pp. 4594–4599.
2. M. Todescato, A. Carron, R. Carli, G. Pillonetto, and L. Schenato, "Efficient Spatio-Temporal Gaussian Regression via Kalman Filtering," arXiv:1705.01485.

The original README is preserved in [UPSTREAM_README.md](UPSTREAM_README.md). Its reference to main.m and its directory layout describe the upstream project, rather than this folder's current layout.

## Licensing

This folder includes third-party code: plotResults.m carries a University of Padova copyright notice (2017), credits Andrea Carron and Marco Todescato, and specifies GPL-3.0-or-later. These notices are preserved, and the GPL version 3 text is provided in [COPYING](COPYING).

No blanket MIT license is asserted over this folder. The provenance and license scope of the remaining files have not been independently established; inclusion here does not replace any applicable upstream terms.
