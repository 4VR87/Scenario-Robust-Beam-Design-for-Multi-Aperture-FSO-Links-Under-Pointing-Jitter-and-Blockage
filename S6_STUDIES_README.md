## Run the numerical studies

Use Python 3.12 and the versions in `requirements_numerical.txt`. From the extracted package directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python revision_experiments.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python channel_diagnostics.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python pink_experiments.py --workers 3
python verify_reported_results.py
```

## Scenario and probability definitions

* Aperture 0 is central; apertures 1–6 are at ring azimuths 0, 60, 120, 180, 240, and 300 degrees.
* S3 contains no pair blockage plus the original pairs (3,4), (1,2), and (5,6). S6 contains no pair blockage plus all six adjacent outer pairs (1,2), (2,3), (3,4), (4,5), (5,6), and (6,1). These are four and seven states, respectively.
* Strong blockage is a retained-power multiplier Uniform[0,0.12] drawn independently for each member of the selected pair in each realization. Every aperture also independently experiences a micro-event with probability 0.01 = 1% per realization conditional on a state. On a micro-event, another independent retained-power multiplier Uniform[0,0.4] applies. Simultaneous factors multiply. This is not a state occurrence probability or a temporal blockage process.
* Scintillation factors are mean-normalized lognormal variables with latent correlation Cij = exp[-(Dij/0.055 m)^(5/3)]. The original source adds 1e-12 times the identity; its negligible effect on analytic mean is disclosed. Positive definiteness and a 200,000-draw mean check are in `model_checks.json`.

## New S6 protocol

Optimization uses 5,000 samples **per state**, not 60,000. DE uses bounds x0,y0 in [-0.1,0.1] m and divergence half-angle in [60,350] microradians, maxiter 65, popsize 11, tol 0.0002, atol 0, best1bin, mutation (0.5,1), recombination 0.7, immediate updates, one worker, and final polishing. Three-variable population size is 33; center-only size is 22. Cached pointing and fading/blockage products retain common random numbers across candidate footprints.

The five joint 5% searches use SAA base seeds 20261024400 + 100000*j and DE seeds 701+j for j=0,...,4. The first search, alias `S6_joint_alpha05`, is retained by a predeclared rule. No independent test score chooses the retained restart. The center-only search holds divergence at 100 microradians. Alpha 0.01 and 0.10, average-SNR, average-outage, and mean-power searches use the same S6 SAA base pool. Every seed and final status is recorded in `optimizations.csv` and individual optimization JSON files. `history_*.csv` records every generation's best score before final polishing. Final nfev includes polishing and excludes callback score computations; elapsed time includes callbacks.

Independent validation has ten trials of 20,000 samples per state, base seed 20261024100 + 100000*(trial-1), and state increments 1009. Candidate designs share streams within each trial/state, enabling a paired joint-minus-center interval. A separate 60,000-sample/state test uses base seed 20261024200. Randomized-sector mixture seeds are 20261024300 (channel) and 20261024301 (uniform pair draw). Mixture lower tails are not minima of conditional lower tails; this experiment introduces per-draw orientation randomness, not unseen shadow shapes.

## Link model and interpretation

The illustrative receiver sums seven equal-weight photodetector currents without coherent optical-field combining, branch selection, or maximal-ratio weighting. All branches contribute independent electronic noise even when blocked; fading remains spatially correlated. For equiprobable OOK, average launch power Pavg and ON power 2*Pavg, fixed transmission Tf, responsivity R, per-branch current-noise density iN, and electrical noise-equivalent bandwidth B:

* sigma_sum^2 = 7*iN^2*B;
* kappa = R*Pavg*Tf/sigma_sum;
* decision-distance SNR gamma_s = (kappa*eta_s)^2;
* conditional BER Pb(eta_s) = Q(kappa*eta_s), using a perfect adaptive midpoint threshold;
* eta_threshold = Q_inverse(0.001)/kappa;
* state outage = Pr(eta_s < eta_threshold), state mean BER = E[Q(kappa*eta_s)].

Nominal parameters: Pavg -15 dBm, Tf 0.10, R 0.9 A/W, iN 5 pA/sqrt(Hz), B 100 MHz. This gives eta_threshold 0.1436374. Geometric loss is already in eta and is not counted again as a free-space path-loss factor. Constant Gaussian receiver noise, uncoded OOK, and threshold knowledge are assumptions, not a measured link budget. Signal-dependent noise, coding, tracking dynamics, and hardware calibration require further system definition.

The S6 average-SNR criterion maximizes the equal-state average E[eta_s^2]; the average-outage criterion minimizes the equal-state average Pr(eta_s < eta_threshold) at -15 dBm. These are objective-level adaptations motivated by references [11] and [13], not replications of those papers' MIMO/UAV architectures or algorithms. Their arithmetic weights are comparison conventions, not measured state probabilities. The empirical outage criterion is discontinuous and only heuristically optimized.

The deep-tail design need not minimize outage at a particular threshold. At -15 dBm it has lower worst-state average BER than center-only design but higher worst-state outage. The ranking changes along the fixed-design power sweep. A 119.07% S3 power-tail increase is therefore not reported as a universal communication or availability gain.

## Files and manuscript mapping

`pink_revision_data/`:

* `model_checks.json`: analytic mean, eigenvalues, Cholesky check, empirical normalization and micro-event frequencies.
* `optimizations.csv`, `optimization_*.json`, `history_*.csv`: eleven new S6 searches, including five independent 5% joint searches, other tail levels, center-only, and three comparison objectives.
* `designs.json`: every retained parameter vector in meters and microradians; historical S3 vectors are distinct from new S6 vectors.
* `independent_trials.csv`, `independent_trial_summary.csv`: ten-trial metrics, sample SD and 95% Student-t intervals with nine degrees of freedom. These quantify Monte Carlo estimation for each fixed design, not population-optimizer or deployment uncertainty.
* `s6_joint_center_comparison.json`: paired difference CI using within-trial common random numbers.
* `independent_state_metrics.csv`: independent 60,000-sample conditional-state checks, including all three tail levels, outage and mean BER.
* `randomized_sector_holdout.csv`: fresh randomized adjacent-sector mixture; no new pair geometry is introduced.
* `link_config.json`, `link_power_sweep.csv`: complete illustrative link parameters and fixed-footprint power curves. Curve estimates use the separate large test and need not equal the ten-trial table means.
* `quantized_control_trials.csv`: continuous versus 1 mm/1 microradian rounded controls on the same trial pools.
* `environment.json`: numerical versions, platform, state order, and DE settings.


