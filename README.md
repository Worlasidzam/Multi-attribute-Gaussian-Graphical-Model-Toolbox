# Multi-Attribute Gaussian Graphical Model: Regularization Selection Toolbox

This repository provides MATLAB implementations for sparse estimation and model selection in multi-attribute Gaussian graphical models (MA-GGMs). The code accompanies the following publications:

**Thesis:**
> W. K. Dzam, "Model Selection for Multi-Attribute Gaussian Graphical Models," M.S. Thesis, Dept. of Electrical and Computer Engineering, Auburn University, Auburn, AL, 2026.

**Conference Paper:**
> W. K. Dzam and J. K. Tugnait, "Regularization Selection for High-Dimensional Multi-Attribute Graph Learning," submitted to *Asilomar Conference on Signals, Systems & Computers*, Pacific Grove, CA, Oct. 2026.

**Estimation Framework:**
> J. K. Tugnait, "Multi-Attribute Graph Estimation with Sparse-Group Non-Convex Penalties," *IEEE Access*, vol. 13, pp. 80174–80190, 2025.

---

## Overview

In a multi-attribute graphical model, each node in a graph is associated with a vector of *m* attributes rather than a single scalar variable. Conditional independence between nodes *i* and *j* corresponds to the entire *m × m* off-diagonal block Ω^(ij) of the precision matrix being zero.

This toolbox implements and compares four regularization selection methods for estimating the graph structure:

1. **BIC** — Bayesian Information Criterion
2. **Cross-Validation (CV)** — K-fold cross-validation using log-likelihood loss
3. **StARS** — Stability Approach to Regularization Selection
4. **StARS + BIC** — Two-stage method: StARS for initial selection, BIC for refinement

Each method is evaluated with two penalty types:
- **Lasso** (ℓ₁) — convex group-sparse penalty
- **LSP** (Log-Sum Penalty) — non-convex sparse-group penalty solved via LLA + ADMM

---

## Key Findings

- **LSP consistently outperforms Lasso** across all selection methods and graph types
- **StARS + BIC (LSP)** achieves the best overall performance, combining stability-based screening with information-theoretic refinement
- **BIC (LSP)** is the best single-stage method
- **CV requires significantly larger sample sizes** (N ≥ 1200) to perform well
- Results are consistent across Erdős–Rényi, Barabási–Albert, and chain graph structures

---

## Installation

```bash
git clone https://github.com/Worlasidzam/Multi-attribute-Gaussian-Graphical-Model-Toolbox.git
cd Multi-attribute-Gaussian-Graphical-Model-Toolbox
```

In MATLAB:
```matlab
addpath(genpath(pwd))
```

---

## Usage

### Synthetic Experiments

Main experiment scripts are located in the `scripts/` directory:

- `run_stability_bic.m` — Stability Selection and BIC experiments
- `run_cv_only.m` — Cross-Validation based model selection

These scripts reproduce results reported in the associated thesis.

### Real Data Experiments

The framework is designed to support real multivariate datasets, including financial time series with multiple attributes per node (e.g., S&P 100 stock data with open, high, low, close prices). All model selection methods implemented here can be applied consistently to both synthetic and real data.

---

## Repository Structure

```text
.
├── scripts/
│   ├── run_stability_bic.m       # Stability Selection and BIC experiment
│   └── run_cv_only.m             # Cross-validation experiments
├── src/
│   └── Core estimation, ADMM, penalties, graph generation
├── data/
│   └── Financial dataset (S&P 100, source: Yahoo Finance)
├── .gitignore
├── LICENSE
└── README.md
```

---

## Requirements

- MATLAB R2020b or later
- Standard MATLAB toolboxes only

---

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `p` | Number of nodes | 100 |
| `m` | Attributes per node | 4 |
| `alpha` | Element-wise vs group penalty weight | 0.05 |
| `epsilon` | LSP non-convexity parameter | 10⁻⁴ |
| `rho` | Initial ADMM penalty parameter | 2 |
| `beta` | StARS instability threshold | 0.05 |
| `p_thr` | StARS pruning threshold | 0.60 |
| `K` | Number of CV folds | 5 |
| `N_sub` | Number of StARS subsamples | 20 |

---

## Performance Metrics

- **F₁-score:** Harmonic mean of precision and recall (higher is better, ideal = 1)
- **Hamming distance:** Total number of false positive + false negative edge decisions (lower is better, ideal = 0)

Results are averaged over 50 Monte Carlo trials with one standard deviation reported.

---

## Citation

If you use this code, please cite:

```bibtex
@mastersthesis{dzam2026thesis,
  author  = {Dzam, Worlasi Kofi},
  title   = {Model Selection for Multi-Attribute Gaussian Graphical Models},
  school  = {Auburn University},
  year    = {2026},
  address = {Auburn, AL},
  type    = {M.S. Thesis}
}

@inproceedings{dzam2026asilomar,
  author    = {Dzam, Worlasi K. and Tugnait, Jitendra K.},
  title     = {Regularization Selection for High-Dimensional Multi-Attribute Graph Learning},
  booktitle = {Proc. Asilomar Conference on Signals, Systems \& Computers},
  year      = {2026},
  address   = {Pacific Grove, CA},
  note      = {Submitted}
}

@article{tugnait2025,
  author  = {Tugnait, Jitendra K.},
  title   = {Multi-Attribute Graph Estimation with Sparse-Group Non-Convex Penalties},
  journal = {IEEE Access},
  volume  = {13},
  pages   = {80174--80190},
  year    = {2025}
}
```

---

## Author

**Worlasi Kofi Dzam**
M.S. in Electrical Engineering, Auburn University, 2026
Email: wkd0014@auburn.edu

---

## Acknowledgments

This work was conducted under the supervision of Prof. Jitendra K. Tugnait. The estimation framework (ADMM solver, LLA scheme) follows [Tugnait, IEEE Access, 2025]. This work was supported by the National Science Foundation under Grant CCF-2308473.

---

## License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.
