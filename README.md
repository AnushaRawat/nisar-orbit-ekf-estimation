# EKF-Based Satellite Orbit Estimation

An orbit-estimation experiment comparing an Extended Kalman Filter (EKF) against noisy six-dimensional satellite-state measurements, using an SGP4-propagated reference trajectory. The filter uses two-body gravity plus Earth's J2 perturbation, fourth-order Runge–Kutta (RK4) propagation, and covariance propagation through the state-transition matrix (STM).

## Key results

**72-hour position RMSE, averaged over 30 independent Monte Carlo realisations per measurement-noise level.** The final comparison uses the same process-noise parameter, `sigma_a = 2e-8 km/s^2`, at every measurement-noise level.

| Position measurement noise, σₚ | Raw measurement RMSE | EKF RMSE | RMSE improvement |
|---:|---:|---:|---:|
| 0.1 km | 0.1722 km | 0.0852 km | 50.5% |
| 0.5 km | 0.8612 km | 0.2883 km | 66.5% |
| 1.0 km | 1.7225 km | 0.4851 km | 71.8% |
| 2.0 km | 3.4450 km | 0.8239 km | 76.1% |

RMSE values are means across seeds; see the full evaluation table below for standard deviations. The percentage improvement is computed per realisation and then averaged, so it may not exactly equal the percentage calculated from the displayed mean RMSE values.

## Objective

Assess whether an EKF can reduce position-estimation error when its measurements contain Gaussian noise, and examine performance as the assumed measurement-noise level changes.

## Method

1. **Reference trajectory:** propagate a satellite orbit using SGP4 and save position and velocity states at 10-minute intervals for 72 hours.
2. **Measurements:** add independent, zero-mean Gaussian perturbations to each position and velocity component. The Monte Carlo experiment uses position-noise standard deviations of 0.1, 0.5, 1.0 and 2.0 km. Velocity-noise standard deviation is set to 1% of the numerical position-noise value, in km/s.
3. **Filter dynamics:** model two-body Earth gravity plus the J2 oblateness perturbation. Atmospheric drag and other perturbations included in SGP4 are not modelled by the EKF.
4. **Numerical propagation:** integrate the six-dimensional state and the state-transition matrix together using RK4, with a 30-second internal integration step.
5. **EKF update:** use the full state measurement `[x, y, z, vx, vy, vz]`, propagate the covariance, and apply the measurement update with a Joseph-form covariance update.
6. **Evaluation:** compare raw measurements and filtered estimates with the SGP4 reference. Report position RMSE over 24-, 48- and 72-hour windows, along with NEES and NIS consistency diagnostics.

All states use kilometres and kilometres per second; time is in seconds internally.

## Process-noise selection

A separate tuning sweep tested the following process-noise parameter values:

`1e-9, 3e-9, 1e-8, 2e-8, 3e-8, 5e-8, 1e-7, 3e-7, 1e-6 km/s²`.

The sweep used 10 tuning seeds per combination. The lowest-RMSE grid value differed across measurement-noise levels: `1e-8 km/s²` was best among the tested values for position noise of 0.1, 0.5 and 1.0 km, while `1e-9 km/s²` was best for 2.0 km. For the final controlled comparison, `2e-8 km/s²` was selected as a common compromise and then held fixed across all four noise levels. It should not be interpreted as a universally optimal value.

## Final evaluation results

Values below are mean ± standard deviation across 30 independent noise realisations. The 24-hour and 48-hour windows are nested subsets of the same 72-hour runs; they are not independent trials. The first 12 time steps are excluded as filter warm-up.

| σₚ (km) | Window | Raw RMSE (km) | EKF RMSE (km) | Improvement (%) |
|---:|---:|---:|---:|---:|
| 0.1 | 24 h | 0.1737 ± 0.0074 | 0.0842 ± 0.0061 | 51.5 ± 3.3 |
| 0.1 | 48 h | 0.1725 ± 0.0049 | 0.0848 ± 0.0046 | 50.8 ± 2.4 |
| 0.1 | 72 h | 0.1722 ± 0.0035 | 0.0852 ± 0.0036 | 50.5 ± 1.8 |
| 0.5 | 24 h | 0.8685 ± 0.0372 | 0.2845 ± 0.0311 | 67.2 ± 3.4 |
| 0.5 | 48 h | 0.8625 ± 0.0244 | 0.2881 ± 0.0240 | 66.6 ± 2.6 |
| 0.5 | 72 h | 0.8612 ± 0.0173 | 0.2883 ± 0.0192 | 66.5 ± 2.1 |
| 1.0 | 24 h | 1.7370 ± 0.0744 | 0.4929 ± 0.0580 | 71.6 ± 3.2 |
| 1.0 | 48 h | 1.7250 ± 0.0489 | 0.4889 ± 0.0486 | 71.7 ± 2.7 |
| 1.0 | 72 h | 1.7225 ± 0.0347 | 0.4851 ± 0.0411 | 71.8 ± 2.3 |
| 2.0 | 24 h | 3.4741 ± 0.1488 | 0.8755 ± 0.1104 | 74.8 ± 3.1 |
| 2.0 | 48 h | 3.4499 ± 0.0977 | 0.8441 ± 0.0972 | 75.5 ± 2.7 |
| 2.0 | 72 h | 3.4450 ± 0.0694 | 0.8239 ± 0.0864 | 76.1 ± 2.4 |

## Consistency diagnostics

For a six-dimensional state, the nominal expected value of NEES and NIS is 6 under the relevant consistency assumptions. With the fixed process-noise parameter, the ensemble-average results were:

| σₚ (km) | Mean NEES | NEES time steps within ensemble bounds | Mean NIS | NIS time steps within ensemble bounds |
|---:|---:|---:|---:|---:|
| 0.1 | 5.15 | 56% | 5.70 | 94% |
| 0.5 | 6.79 | 49% | 5.85 | 95% |
| 1.0 | 7.35 | 41% | 5.90 | 95% |
| 2.0 | 7.02 | 51% | 5.92 | 95% |

The NIS means are close to 6. NEES is also of a similar scale, though the 1.0 km case is slightly above the stated ensemble-average 95% bounds of 4.82–7.30. These are diagnostics, not proof of perfect consistency; covariance calibration and dynamical model mismatch remain relevant limitations.

## Numerical checks

The implementation includes checks for:

- **Analytical acceleration Jacobian:** maximum relative error against central finite differences of `5.66e-10` in the recorded run.
- **RK4 step-size convergence:** position difference of `3.58e-5 km` between 30-second and 5-second internal steps over one measurement interval; `4.10e-7 km` between 10-second and 5-second steps.
- **One-step model mismatch:** mean 7.3 m and maximum 11.3 m between the J2-based EKF propagation and the SGP4 reference over a 600-second interval in the recorded run.

These checks support the numerical implementation, but do not eliminate the physical differences between the filter dynamics and SGP4.

## Limitations

- This is a simulation experiment, not an operational orbit-determination system.
- SGP4 supplies the reference trajectory, while the EKF uses a simplified gravity-plus-J2 model. Their dynamical mismatch can contribute to estimation error.
- Measurements are simulated with independent Gaussian noise; real tracking errors can be correlated, biased, or non-Gaussian.
- The chosen process-noise parameter is a common compromise for this comparison, not a universal optimum.
- Regenerating the orbit using a newer TLE may change the reference trajectory and results.

## Repository contents

- `01_orbit.ipynb` — generates the SGP4 reference orbit and noisy input files.
- `02_ekf_clean.ipynb` — notebook walkthrough of the fixed-process-noise experiment.
- `ekf_orbit.py` — main EKF implementation and experiment runner.
- `data/README.md` — required input files, formats and regeneration instructions.
- `requirements.txt` — Python dependencies.

## Setup and run

Use Python 3.10 or later. From the repository root, create and activate a virtual environment if desired, then install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
python -m pip install -r requirements.txt
```

Run `01_orbit.ipynb` first to generate the required data files in `data/`. The orbit notebook fetches a TLE from CelesTrak, so an internet connection is required for that step. Then run the experiment:

```bash
python ekf_orbit.py
```

The script writes `results/summary.json`, figures in `figures/`, and a filtered state file at `data/ekf_orbit_corrected.npz`. Generated data and outputs are excluded from Git by default.
