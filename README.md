# Scenario-Robust-Beam-Design-for-Multi-Aperture-FSO-Links-Under-Pointing-Jitter-and-Blockage
This work presents a scenario-robust beam design for multi-aperture FSO links under anisotropic pointing jitter and partial blockage. The method jointly optimizes beam-center offset and divergence to maximize the worst-state 5% lower-tail collected power, improving robustness while quantifying the trade-off with mean optical collection.

## Software Requirement

numpy==2.3.5
scipy==1.17.0
pandas==2.2.3
matplotlib==3.10.7

  "test_samples_per_state": 60000,
  "repeat_trials": 10,
  "repeat_samples_per_state": 20000,
  "base_seed": 20260919

  green_experiments.py` drives both new studies. `green_independent_solver.py` is isolated from the original model's sampling, coupling, and tail routines. It implements:

- direct two-dimensional aperture-intensity integration, with the squared-radius change of variable given in Eq. (13);
- 32 radial Gauss–Legendre nodes and 128 periodic azimuths, checked against 16 × 64 and 64 × 256 rules;
- a PCHIP coupling table on 0–550 mm displacement, at 0.2 mm spacing, with extrapolation prohibited;
- symmetric eigenvalue-based covariance roots rather than the original Cholesky path;
- independently scrambled 30-dimensional Sobol sampling for the complete channel;
- full sample sorting for each tail mean, rather than the original partition operation.

The independent route uses the same Gaussian envelope, pointing distribution, correlated lognormal factors, retained-power blockage multipliers, geometry, and fixed parameter vectors. Its agreement therefore verifies numerical implementation, not those physical assumptions.

The independent comparison uses four retained footprints on all seven S6 states, with 4,096, 16,384, and 65,536 samples per state and eight independent replications at every budget. Base seeds are 202610240600 for ordinary Monte Carlo and 202610240700 for scrambled Sobol; replication increments are 1,000,000 and state increments 1009. Budgets reuse the same scramble prefix. A separate paired-input comparison uses 8,192 samples per state and seeds 202610240500 + 1009 × state index. There are 276 off-grid quadrature checks for each footprint.

At 65,536 samples per state, the S6 joint-tail means are 0.08791184 for analytic coupling/Monte Carlo and 0.08781993 for independent integration/Sobol (rounded values; consult CSVs for full precision). The relative difference is −0.105%. The maximum absolute relative mean difference across the four fixed footprints is 0.504%. Every independent mean falls within its corresponding ordinary-Monte-Carlo 95% interval at this budget. The same footprint ordering holds for both methods at all three budgets. These observations do not establish distributional equivalence or a universal quasi-Monte-Carlo convergence rate for the tail/minimum functional.

The parameter sweep holds the legacy S3 joint, S6 mean-power, and S6 joint-tail footprints fixed. Each of eleven one-at-a-time settings has six independent 20,000-sample/state S6 trials. Base seeds are 202610240900 + 1,000,000 × (trial − 1), with state increments 1009. Underlying draws are shared across cases and designs within each trial/state. Tested settings are:

| Parameter | Values |
| --- | --- |
| Log-irradiance standard deviation σln | 0.15, 0.38, 0.65 |
| Latent correlation length rc | 20, 55, 120 mm |
| Selected pair's retained-power upper endpoint | 0, 0.12, 0.30 |
| Micro-event probability per aperture/realization | 0%, 1%, 5% |
| Micro-event retained-power upper endpoint | 0.1, 0.4, 0.8 |

Only one parameter changes in each non-nominal case. These ranges are controlled design stressors, not measured uncertainty intervals. The S6 joint-tail mean spans 0.06934–0.09707 and its ratio-of-means gain over the retained mean-power design spans 48.3%–161.8%. No simultaneous-uncertainty guarantee is inferred.


