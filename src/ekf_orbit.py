"""
EKF-based satellite orbit estimation experiment.

Reference: SGP4 states in data/truth_orbit.npz
Measurements: Gaussian perturbations stored in data/noisy_orbit.npz
Filter dynamics: two-body gravity + J2, integrated with RK4.
Units: km, seconds, km/s, and km/s^2.
"""

import json
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.stats import chi2

# Experiment settings
MU = 398600.4418          
RE = 6378.137             
J2 = 1.08263e-3
INTEGRATOR_STEP = 30.0    # Internal RK4 step (s)
WARMUP_STEPS = 12         # Excluded from RMSE and consistency metrics
VEL_NOISE_RATIO = 0.01    # sigma_v = 0.01 * sigma_p, in corresponding units
NOISE_LEVELS = [0.1, 0.5, 1.0, 2.0]  # Position measurement sigma (km)
SIGMA_A = 2e-8            # Common process-noise parameter (km/s^2)
N_EVAL_SEEDS = 30
EVAL_SEED0 = 2000

TRUTH_PATH = "data/truth_orbit.npz"
NOISY_PATH = "data/noisy_orbit.npz"
OUT_DIR = "results"
FIG_DIR = "figures"

_D = np.array([1.0, 1.0, 3.0])


# ORBIT DYNAMICS AND STATE TRANSITION MATRIX

def acceleration(r):
    """Two-body gravity plus J2, neglected drag."""
    rho = np.linalg.norm(r)
    a = -MU * r / rho**3
    s = (r[2] / rho) ** 2
    c = 1.5 * J2 * MU * RE**2 / rho**5
    return a + c * r * (5.0 * s - _D)       

def acceleration_jacobian(r):
    rho = np.linalg.norm(r)
    I = np.eye(3)
    rrT = np.outer(r, r)
    z = r[2]
    ez = np.array([0.0, 0.0, 1.0])

    d_twobody = -MU * (I / rho**3 - 3.0 * rrT / rho**5)

    k0 = 1.5 * J2 * MU * RE**2
    dT1 = 5.0 * (np.outer(r, 2.0 * z * ez) / rho**7
                 + z * z * I / rho**7
                 - 7.0 * z * z * rrT / rho**9)
    dT2 = np.diag(_D) / rho**5 - 5.0 * np.outer(_D * r, r) / rho**7
    return d_twobody + k0 * (dT1 - dT2)

def f_state(x):
    return np.concatenate([x[3:6], acceleration(x[:3])])

def F_matrix(x):
    F = np.zeros((6, 6))
    F[:3, 3:] = np.eye(3)
    F[3:, :3] = acceleration_jacobian(x[:3])
    return F

def _aug_deriv(y):
    x = y[:6]
    Phi = y[6:].reshape(6, 6)
    return np.concatenate([f_state(x), (F_matrix(x) @ Phi).ravel()])

def propagate(x, dt, h=INTEGRATOR_STEP):
    """RK4 propagation of state and STM together. Returns (x_next, Phi)."""
    n = int(np.ceil(dt / h))
    hh = dt / n
    y = np.concatenate([np.asarray(x, float), np.eye(6).ravel()])
    for _ in range(n):
        k1 = _aug_deriv(y)
        k2 = _aug_deriv(y + 0.5 * hh * k1)
        k3 = _aug_deriv(y + 0.5 * hh * k2)
        k4 = _aug_deriv(y + hh * k3)
        y = y + (hh / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y[:6], y[6:].reshape(6, 6)

def Q_matrix(dt, sigma_a):
    """Discrete random-acceleration process noise (sigma_a in km/s^2)."""
    q = sigma_a**2
    Q_axis = q * np.array([[dt**4 / 4.0, dt**3 / 2.0],
                           [dt**3 / 2.0, dt**2]])
    Q = np.zeros((6, 6))
    for i in range(3):
        Q[np.ix_([i, i + 3], [i, i + 3])] = Q_axis
    return Q

def run_ekf(times, meas, sigma_p, sigma_v, sigma_a, truth=None):
    N = len(times)
    R = np.diag([sigma_p**2] * 3 + [sigma_v**2] * 3)
    I6 = np.eye(6)

    x = meas[0].copy()
    P = R.copy()

    states = np.zeros((N, 6))
    sigma = np.zeros((N, 6))
    pos_rms_sigma = np.zeros(N)
    nees = np.full(N, np.nan)
    nis = np.full(N, np.nan)

    states[0] = x
    sigma[0] = np.sqrt(np.diag(P))
    pos_rms_sigma[0] = np.sqrt(np.trace(P[:3, :3]))

    for k in range(1, N):
        dt = times[k] - times[k - 1]

        # predict
        x_pred, Phi = propagate(x, dt)
        P_pred = Phi @ P @ Phi.T + Q_matrix(dt, sigma_a)
        P_pred = 0.5 * (P_pred + P_pred.T)

        # update (H = I)
        nu = meas[k] - x_pred
        S = P_pred + R
        S = 0.5 * (S + S.T)
        K = solve(S, P_pred.T, assume_a="pos").T          # P_pred S^-1
        x = x_pred + K @ nu

        IK = I6 - K
        P = IK @ P_pred @ IK.T + K @ R @ K.T
        P = 0.5 * (P + P.T)

        states[k] = x
        sigma[k] = np.sqrt(np.clip(np.diag(P), 0.0, None))
        pos_rms_sigma[k] = np.sqrt(max(np.trace(P[:3, :3]), 0.0))
        nis[k] = nu @ solve(S, nu, assume_a="pos")
        if truth is not None:
            e = x - truth[k]
            nees[k] = e @ solve(P, e, assume_a="pos")

    return dict(states=states, sigma=sigma, pos_rms_sigma=pos_rms_sigma,
                nees=nees, nis=nis)

def pos_err(a, b):
    return np.linalg.norm(a[:, :3] - b[:, :3], axis=1)

def windowed_rmse(err, dt, warmup=WARMUP_STEPS):
    """RMSE over [warmup, window_end). NOTE: windows are nested, not independent."""
    out = {}
    for hrs in (24, 48, 72):
        end = min(int(round(hrs * 3600.0 / dt)) + 1, len(err))
        seg = err[warmup:end]
        out[hrs] = float(np.sqrt(np.mean(seg**2))) if len(seg) else np.nan
    return out

def make_measurements(truth, sigma_p, sigma_v, seed):
    rng = np.random.default_rng(seed)
    scale = np.array([sigma_p] * 3 + [sigma_v] * 3)
    return truth + rng.standard_normal(truth.shape) * scale

def simulate_run(truth, times, sigma_p, sigma_v, sigma_a, seed):
    """One Monte Carlo realisation -> lightweight result dict."""
    dt = float(np.median(np.diff(times)))
    meas = make_measurements(truth, sigma_p, sigma_v, seed)
    res = run_ekf(times, meas, sigma_p, sigma_v, sigma_a, truth=truth)
    raw = windowed_rmse(pos_err(meas, truth), dt)
    ekf = windowed_rmse(pos_err(res["states"], truth), dt)
    return dict(raw=raw, ekf=ekf, nees=res["nees"], nis=res["nis"])

def aggregate(runs, warmup=WARMUP_STEPS):
    """Mean/std of RMSE, improvement, and ensemble-averaged NEES/NIS consistency."""
    n = len(runs)
    out = {"n_runs": n}
    for hrs in (24, 48, 72):
        raw = np.array([r["raw"][hrs] for r in runs])
        ekf = np.array([r["ekf"][hrs] for r in runs])
        imp = (raw - ekf) / raw * 100.0
        out[hrs] = dict(raw_mean=raw.mean(), raw_std=raw.std(ddof=1) if n > 1 else 0.0,
                        ekf_mean=ekf.mean(), ekf_std=ekf.std(ddof=1) if n > 1 else 0.0,
                        imp_mean=imp.mean(), imp_std=imp.std(ddof=1) if n > 1 else 0.0)

    lo = chi2.ppf(0.025, 6 * n) / n
    hi = chi2.ppf(0.975, 6 * n) / n
    for name in ("nees", "nis"):
        arr = np.array([r[name] for r in runs])[:, warmup:]    # (n, T)
        ens = np.nanmean(arr, axis=0)                           # ensemble avg per step
        out[name] = dict(mean=float(np.nanmean(ens)),
                         frac_in_bounds=float(np.mean((ens >= lo) & (ens <= hi))),
                         lo=lo, hi=hi)
    return out

def check_jacobian(x):
    """Analytic acceleration Jacobian vs central finite differences."""
    r = x[:3]
    num = np.zeros((3, 3))
    for i in range(3):
        d = np.zeros(3)
        d[i] = 1e-3
        num[:, i] = (acceleration(r + d) - acceleration(r - d)) / 2e-3
    ana = acceleration_jacobian(r)
    err = np.max(np.abs(ana - num)) / np.max(np.abs(num))
    print(f"[check] Jacobian max relative error vs finite diff: {err:.2e}")
    assert err < 1e-5, "Analytic Jacobian does not match finite differences"

def check_step_size(x0, dt):
    xs = {h: propagate(x0, dt, h=h)[0] for h in (30.0, 10.0, 5.0)}
    print("[check] RK4 step-size convergence (position difference, km):")
    print(f"        30 s vs 5 s : {np.linalg.norm(xs[30.0][:3] - xs[5.0][:3]):.2e}")
    print(f"        10 s vs 5 s : {np.linalg.norm(xs[10.0][:3] - xs[5.0][:3]):.2e}")

def check_one_step_model_error(truth, dt):
    """Truth(k) propagated one step vs truth(k+1): the model-mismatch floor."""
    e = [np.linalg.norm(propagate(truth[k], dt)[0][:3] - truth[k + 1, :3])
         for k in range(len(truth) - 1)]
    print(f"[check] One-step model error vs truth: mean={np.mean(e) * 1e3:.1f} m  "
          f"max={np.max(e) * 1e3:.1f} m")
    print("        (difference between the EKF propagation model and the SGP4 reference)")

def make_plots(times, truth, eval_table):
    """Create the plots used to summarise the fixed-process-noise experiment."""
    os.makedirs(FIG_DIR, exist_ok=True)
    hours = times / 3600.0

    # Position error over time for one reproducible measurement realisation.
    sigma_p = NOISE_LEVELS[0]
    sigma_v = sigma_p * VEL_NOISE_RATIO
    measurements = make_measurements(truth, sigma_p, sigma_v, EVAL_SEED0)
    raw_error = pos_err(measurements, truth)
    filtered = run_ekf(times, measurements, sigma_p, sigma_v, SIGMA_A, truth=truth)
    ekf_error = pos_err(filtered["states"], truth)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(hours, raw_error, color="0.65", linewidth=1, label="Raw measurement error")
    ax.plot(hours, ekf_error, linewidth=1.3, label=f"EKF error (sigma_a={SIGMA_A:.0e})")
    ax.axvspan(0, WARMUP_STEPS * (times[1] - times[0]) / 3600,
               color="k", alpha=0.06, label="Warm-up (excluded from metrics)")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("3-D position error (km)")
    ax.set_title(f"Position error over time (position noise = {sigma_p} km)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "error_vs_time.png"), dpi=150)
    plt.close(fig)

    # Improvement over raw measurements at each measurement-noise level.
    means = [eval_table[s]["agg"][72]["imp_mean"] for s in NOISE_LEVELS]
    stds = [eval_table[s]["agg"][72]["imp_std"] for s in NOISE_LEVELS]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(NOISE_LEVELS, means, yerr=stds, fmt="o-", capsize=4)
    ax.set_xscale("log")
    ax.set_xlabel("Position measurement noise sigma_p (km)")
    ax.set_ylabel("Position RMSE improvement (%)")
    ax.set_title(f"EKF improvement over 72 h ({N_EVAL_SEEDS} seeds)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "improvement_vs_noise.png"), dpi=150)
    plt.close(fig)
    print(f"[plots] Saved figures to {FIG_DIR}/")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    print("[EKF] Loading data...")
    with np.load(TRUTH_PATH) as data:
        truth = np.asarray(data["states"], dtype=float)
    with np.load(NOISY_PATH) as data:
        file_measurements = np.asarray(data["states"], dtype=float)
        times = np.asarray(data["times_minutes"], dtype=float) * 60.0

    if truth.shape != file_measurements.shape or len(times) != len(truth):
        raise ValueError("Truth, measurement, and time arrays must have matching lengths.")
    if not np.all(np.diff(times) > 0):
        raise ValueError("Times must be strictly increasing.")

    dt = float(np.median(np.diff(times)))
    print(f"[EKF] {len(times)} epochs, measurement interval = {dt:.0f} s")
    print(f"[EKF] Fixed sigma_a = {SIGMA_A:.0e} km/s^2")
    print(f"[EKF] Evaluation: {N_EVAL_SEEDS} seeds per noise level")

    print("\nNumerical checks")
    check_jacobian(truth[0])
    check_step_size(truth[0], dt)
    check_one_step_model_error(truth, dt)

    evaluation = {}
    start = time.time()
    print("\n" + "=" * 78)
    print("FIXED PROCESS-NOISE EXPERIMENT")
    print("=" * 78)

    for sigma_p in NOISE_LEVELS:
        sigma_v = sigma_p * VEL_NOISE_RATIO
        runs = []
        for seed in range(EVAL_SEED0, EVAL_SEED0 + N_EVAL_SEEDS):
            runs.append(simulate_run(truth, times, sigma_p, sigma_v, SIGMA_A, seed))
        evaluation[sigma_p] = {"sigma_a": SIGMA_A, "agg": aggregate(runs)}

    print(f"Position RMSE (km), mean ± std across seeds; warm-up = {WARMUP_STEPS} steps.")
    print("The 24 h and 48 h windows are nested within the 72 h runs, not independent trials.\n")
    print(f"{'sigma_p':>8} {'window':>7} {'raw RMSE':>16} {'EKF RMSE':>16} {'improvement %':>15}")
    for sigma_p in NOISE_LEVELS:
        agg = evaluation[sigma_p]["agg"]
        for hours in (24, 48, 72):
            row = agg[hours]
            print(f"{sigma_p:>8.1f} {str(hours) + ' h':>7} "
                  f"{row['raw_mean']:>8.4f}±{row['raw_std']:.4f} "
                  f"{row['ekf_mean']:>8.4f}±{row['ekf_std']:.4f} "
                  f"{row['imp_mean']:>7.1f}±{row['imp_std']:.1f}")

    n = N_EVAL_SEEDS
    lo = chi2.ppf(0.025, 6 * n) / n
    hi = chi2.ppf(0.975, 6 * n) / n
    print(f"\nConsistency diagnostics (ensemble-average 95% bounds: [{lo:.2f}, {hi:.2f}], expected value 6)")
    print(f"{'sigma_p':>8} {'mean NEES':>11} {'NEES in bounds':>16} {'mean NIS':>10} {'NIS in bounds':>15}")
    for sigma_p in NOISE_LEVELS:
        stats = evaluation[sigma_p]["agg"]
        print(f"{sigma_p:>8.1f} {stats['nees']['mean']:>11.2f} "
              f"{stats['nees']['frac_in_bounds'] * 100:>15.0f}% "
              f"{stats['nis']['mean']:>10.2f} "
              f"{stats['nis']['frac_in_bounds'] * 100:>14.0f}%")

    make_plots(times, truth, evaluation)

    def to_json(value):
        if isinstance(value, dict):
            return {str(k): to_json(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [to_json(v) for v in value]
        if isinstance(value, (np.floating, np.integer)):
            return value.item()
        return value

    summary = {
        "config": {
            "warmup_steps": WARMUP_STEPS,
            "noise_levels_km": NOISE_LEVELS,
            "sigma_a_km_per_s2": SIGMA_A,
            "evaluation_seeds": N_EVAL_SEEDS,
            "evaluation_seed_start": EVAL_SEED0,
            "velocity_noise_ratio": VEL_NOISE_RATIO,
            "integrator_step_seconds": INTEGRATOR_STEP,
        },
        "evaluation": evaluation,
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(to_json(summary), f, indent=2)

    # Filter the measurements saved by the orbit-generation notebook.
    filtered = run_ekf(
        times, file_measurements, NOISE_LEVELS[0],
        NOISE_LEVELS[0] * VEL_NOISE_RATIO, SIGMA_A, truth=truth
    )
    np.savez(
        os.path.join("data", "ekf_orbit_corrected.npz"),
        states=filtered["states"],
        sigma=filtered["sigma"],
        nees=filtered["nees"],
        nis=filtered["nis"],
        times_seconds=times,
        times_minutes=times / 60.0,
        sigma_a=SIGMA_A,
        sigma_p=NOISE_LEVELS[0],
        sigma_v=NOISE_LEVELS[0] * VEL_NOISE_RATIO,
    )
    print(f"\n[EKF] Saved {OUT_DIR}/summary.json and data/ekf_orbit_corrected.npz")
    print(f"[EKF] Total runtime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
