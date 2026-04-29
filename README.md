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

This toolbox implements and compares three regularization selection methods for estimating the graph structure:

1. **Stability Selection (SS)** — Stability Approach to Regularization Selection (StARS), with optional pruning via selection probabilities
2. **BIC** — Bayesian Information Criterion
3. **Cross-Validation (CV)** — K-fold cross-validation using log-likelihood loss

Each method is evaluated with two penalty types:
- **Lasso** (ℓ₁) — convex group-sparse penalty
- **LSP** (Log-Sum Penalty) — non-convex sparse-group penalty solved via LLA + ADMM

---

## Key Findings

- **LSP consistently outperforms Lasso** across all selection methods and graph types
- **BIC (LSP)** provides the best overall balance between accuracy and sparsity at moderate-to-large sample sizes
- **SS (LSP) with pruning** is highly effective in sparse or small-sample regimes
- **CV requires significantly larger sample sizes** (N ≥ 1200) to perform well, but achieves near-perfect recovery with LSP at large N
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

Scripts in `scripts/`:

- `run_stability_bic.m` — Stability Selection and BIC experiments
- `run_cv_only.m` — Cross-Validation based model selection

### Real Data Experiments

Scripts in `scripts/`:

- `multiAttr_real_ss_all_new.m` — Stability Selection on S&P 100 data
- `multiAttr_real_bic_all_new.m` — BIC on S&P 100 data
- `multiAttr_real_cv_all_new.m` — Cross-Validation on S&P 100 data

The S&P 100 dataset consists of 97 stocks with four daily attributes (high, low, close, volume) over 1259 trading days (April 2015 – April 2020), sourced from Yahoo Finance.

---

## Repository Structure

```text
.
├── scripts/
│   ├── run_stability_bic.m             # SS and BIC synthetic experiments
│   ├── run_cv_only.m                   # CV synthetic experiments
│   ├── multiAttr_real_ss_all_new.m     # SS on S&P 100 data
│   ├── multiAttr_real_bic_all_new.m    # BIC on S&P 100 data
│   └── multiAttr_real_cv_all_new.m     # CV on S&P 100 data
├── src/
│   ├── BAModel_mod.m                   # Barabási–Albert graph generator
│   ├── GenGraphPrec.m                  # Graph and precision matrix generation
│   ├── bic_selec.m                     # BIC selection
│   ├── bisection_uplim.m              # Lambda grid upper limit via bisection
│   ├── dataGen.m                       # Gaussian data generation
│   ├── opt_admm_ma.m                   # ADMM solver
│   ├── opt_admm_ma_adap.m              # ADMM solver (LSP/adaptive)
│   ├── optimize_admm_ma.m              # ADMM optimizer wrapper
│   ├── performance.m                   # F1-score and Hamming distance
│   └── stab_selec_modv4.m              # Stability selection
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
