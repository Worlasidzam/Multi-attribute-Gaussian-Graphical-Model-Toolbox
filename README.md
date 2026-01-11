# Multi-attribute Gaussian Graphical Model Toolbox

This repository provides MATLAB implementations for sparse estimation and model selection in multi-attribute Gaussian graphical models (MA-GGMs). The code accompanies a Master’s thesis and follows the methodology developed in:

J. K. Tugnait,
“Multi-Attribute Graph Estimation With Sparse-Group Non-Convex Penalties,”
IEEE Access, 2025.

The focus is on high-dimensional settings where the number of variables is large relative to the sample size and regularization is required for stable, interpretable graph recovery.

---

## Overview

In a multi-attribute graphical model, each node in a graph is associated with a vector of attributes rather than a single scalar variable. Conditional dependencies between nodes are encoded through the block structure of the precision matrix.

This toolbox supports:
- Sparse precision matrix estimation
- Conditional independence graph recovery
- Model selection under small-sample, high-dimensional regimes
- Evaluation on synthetic graph models and real data

---

## Mathematical Setting

Let

x = [z₁ᵀ z₂ᵀ … z_pᵀ]ᵀ ∈ ℝ^(mp)

be a zero-mean Gaussian random vector with covariance Σ and precision matrix Ω = Σ⁻¹.

Edges in the graph correspond to nonzero off-diagonal blocks of Ω, encoding conditional dependence between nodes.

---

## Methods Implemented

### Graph Estimation
- Penalized Gaussian log-likelihood formulation
- ADMM-based optimization
- Element-wise and block-wise sparsity control

### Regularization
- Lasso (ℓ₁ penalty)
- Log-Sum Penalty (LSP)
- Sparse-group penalties with mixing parameter α

### Model Selection
- Stability Selection (with and without pruning)
- Bayesian Information Criterion (BIC)
- Cross-Validation (CV)

---

## Graph Models (Synthetic Data)

The toolbox supports controlled generation of synthetic graphs:
- Erdős–Rényi (ER) graphs
- Barabási–Albert (BA) graphs
- Chain graphs

Ground-truth precision matrices are constructed to ensure positive definiteness, controlled sparsity, and realistic intra-node attribute correlation (AR-type decay).

---

## Real Data Experiments

The framework is designed to support real multivariate datasets, including financial time series with multiple attributes per node. All model selection methods implemented here can be applied consistently to both synthetic and real data.

---
## Repository Structure

```text
.
├── scripts/
│ ├── run_stability_bic.m # Stability Selection and BIC experiment
│ └── run_cv_only.m # Cross-validation experiments
├── src/
│ └── Core estimation, ADMM, penalties, graph generation
├── .gitignore
├── LICENSE
└── README.md
```



## Requirements

- MATLAB R2020b or later
- Standard MATLAB toolboxes only

---

## Usage

After cloning the repository, add it to the MATLAB path:

addpath(genpath(pwd))

Main experiment scripts are located in the scripts directory:
- run_stability_bic.m – Stability Selection and BIC
- run_cv_only.m – Cross-Validation based model selection

These scripts reproduce results reported in the associated thesis.

---

## Status

This repository is under active development. Code structure and documentation are being refined as part of ongoing research.

---

## Citation

If you use this code, please cite:

J. K. Tugnait,
“Multi-Attribute Graph Estimation With Sparse-Group Non-Convex Penalties,”
IEEE Access, 2025.

and the associated Master’s thesis:

Model Selection for Multi-Attribute Gaussian Graphical Models
Auburn University, 2026

---

## Author

Worlasi Kofi Dzam
Master’s Student, Electrical Engineering
Auburn University

---

## License

This project is released under the MIT License. See the LICENSE file for details.

---

## Contact

Email: wkd0014@auburn.edu

For questions or collaboration inquiries, please open an issue or contact the author directly.
