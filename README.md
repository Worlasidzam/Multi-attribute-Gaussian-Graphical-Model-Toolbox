# Multi-attribute Gaussian Graphical Model Estimation.
Multi-attribute Gaussian graphical model estimation and model selection for high dimensional learning and Inference.

This repository provides MATLAB implementations for sparse estimation and model selection in multi-attribute Gaussian graphical models (MA-GGMs). The code accompanies a Master’s thesis and follows the methodology developed in:

J. K. Tugnait, "Multi-Attribute Graph Estimation With Sparse-Group Non-Convex Penalties," IEEE Access, 2025.

The focus is on high-dimensional settings where the number of variables is large relative to the sample size and regularization is required for stable and interpretable graph recovery.

---

## Overview

In a multi-attribute graphical model, each node in a graph is associated with a vector of attributes rather than a single scalar variable. Conditional dependencies between nodes are encoded through the block structure of the precision matrix.

This library supports:
- Sparse precision matrix estimation
- Conditional independence graph recovery
- Model selection using multiple criteria
- Evaluation on synthetic and real-world datasets

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
- Stability Selection (pruned and unpruned)
- Bayesian Information Criterion (BIC)
- Cross-Validation (CV)

---

## Graph Models (Synthetic Data)

The library supports synthetic graph generation for controlled evaluation:
- Erdős–Rényi (ER) graphs
- Barabási–Albert (BA) graphs
- Chain graphs 

Ground-truth precision matrices are constructed to ensure positive definiteness, controlled sparsity, and realistic intra-node attribute correlation.

---

## Real Data Experiments

The framework is designed to support real multivariate datasets, including financial time series and other multi-attribute measurements. All model selection methods can be applied consistently across synthetic and real data.

---

## Repository Structure (Planned)

.
├── src/                 Core estimation and optimization routines  
├── model_selection/     BIC, CV, Stability Selection  
├── graph_generation/    ER, BA, Chain graph generators  
├── experiments/         Scripts reproducing thesis results  
├── examples/            Minimal working examples  
└── README.md  

The codebase is currently being cleaned and modularized for public release.

---

## Requirements

- MATLAB R2020b or later
- Standard MATLAB toolboxes only

---

## Usage

After cloning the repository, add it to the MATLAB path:

addpath(genpath(pwd))

Example experiment scripts reproduce results reported in the associated thesis. Additional documentation and examples will be expanded in future revisions.

---

## Status

This repository is under active development. Code organization, documentation, and examples are being refined.

---

## Citation

If you use this code, please cite:

J. K. Tugnait,  
"Multi-Attribute Graph Estimation With Sparse-Group Non-Convex Penalties,"  
IEEE Access, 2025.

and the associated Master's thesis:
"Model Selection for Multi-attribute Gaussian Graphical Models"
Auburn University, 2026

---

## Author

Worlasi Kofi Dzam  
Master's Student, Electrical Engineering  
Auburn University

---

## License

A permissive open-source license will be added in a future release.

---

## Contact
email: wkd0014@auburn.edu

For questions or collaboration inquiries, please open an issue or contact the author directly.
