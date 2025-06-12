# Proximal-MCMC-Quantile-TrendFiltering

This project implements a fully Bayesian framework for **quantile trend filtering** using Proximal MCMC methods. It replaces classical ADMM point estimation with posterior sampling under a nonsmooth prior.

---

## 📌 Motivation

Traditional trend filtering fails under:
- Non-Gaussian noise
- Sensor drift
- Spiky, asymmetric residuals

We target baseline drift estimation (e.g., in air quality data) using **quantile regression** with:
- Smoothed quantile loss ($\tilde{\rho}_\tau$)
- $\ell_1$ penalty on $(k+1)$-order differences (TV prior)
- Full Bayesian inference via **Proximal ULA and Proximal MALA**

---

## 🧮 Model

We consider the model:

$$
\pi(\theta \mid y) \propto \exp\left(-\sum_{i=1}^n \tilde{\rho}_\tau(y_i - \theta_i)\right) \cdot \exp\left(-\lambda \|D^{k+1} \theta\|_1\right)
$$

where:
- $\tilde{\rho}_\tau$ is a smooth approximation to the quantile loss
- $D^{k+1}$ is the $(k+1)$-order difference matrix (e.g., total variation)
- $\lambda$ controls regularization

---

## ⚙️ Samplers

We implement:

- **P-ULA**: Proximal Unadjusted Langevin Algorithm  
- **P-MALA**: Proximal Metropolis-Adjusted Langevin Algorithm

Each step involves a **proximal update**:

$$
\text{prox}_{\lambda g}(z) = \arg\min_\theta \left\{ \lambda \|D^{k+1} \theta\|_1 + \frac{1}{2\lambda} \|\theta - z\|^2 \right\}
$$

which we solve via ADMM or generalized lasso.

---

## 📈 Results

- Posterior median tracks true baseline
- Posterior credible bands widen near volatility
- Better MSE and uncertainty quantification than deterministic ADMM
- P-MALA outperforms P-ULA when smoothness is low

---

## 🔍 Diagnostics

- Trace plots and autocorrelations verify mixing
- Posterior histograms symmetric and concentrated
- Effective Sample Sizes (ESS) $> 300$ for all parameters

---

## 📁 Files

| File                  | Description                              |
|-----------------------|------------------------------------------|
| `PULA_sampler.R`      | Proximal ULA implementation              |
| `PMALA_sampler.R`     | Proximal MALA implementation             |
| `prox_tv.R`           | Proximal operator for TV penalty         |
| `generate_data.R`     | Synthetic data generation (sinusoidal + noise) |
| `plot_results.R`      | Posterior median, bands, and diagnostics |
| `qtf_admm.R`          | ADMM benchmark implementation            |

---

## 🧠 References

- Pereyra et al. (2016), *Proximal MCMC Algorithms*
- Durmus & Moulines (2018), *High-dimensional Bayesian inference via Proximal MALA*
- Heng et al. (2023), *Bayesian Quantile Trend Filtering*

---

## 👤 Author

**Felipe Bedoya**  
PhD Candidate in Statistics, Rice University  
📧 fb31@rice.edu  
🌐 [GitHub Profile](https://github.com/FbStatsQuant)
