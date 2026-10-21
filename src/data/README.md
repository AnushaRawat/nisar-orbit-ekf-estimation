# Input data

The orbit-generation notebook (`sgp4_orbit.py`) creates the input files required by `ekf_orbit.py`. Run the orbit-generation notebook first, from the repository root, so the files are saved in this directory.

## Required files

### `truth_orbit.npz`

Reference trajectory propagated with SGP4. The EKF script reads:

- `states`: array of shape `(N, 6)`, with each row `[x, y, z, vx, vy, vz]`.
- `times_minutes`: array of shape `(N,)`, containing minutes from the start epoch.

The notebook also saves positions, velocities and Julian dates for inspection.

### `noisy_orbit.npz`

A saved noisy measurement sequence. The EKF script reads:

- `states`: array of shape `(N, 6)`, in the same state order and units as the reference.
- `times_minutes`: array of shape `(N,)`, matching the reference timestamps.

The saved noisy sequence is used for the final single-run filtered output. The Monte Carlo evaluation generates fresh Gaussian measurements for each evaluation seed and each configured noise level.

## Data and reproducibility notes

- The expected experiment uses 433 states at a 10-minute cadence over 72 hours.
- Position is in kilometres and velocity is in kilometres per second.
- The notebook fetches current TLE data from CelesTrak. Because TLEs can change over time, regenerating the data later may produce a different reference trajectory and different numerical results.
- The generated `.npz` and `.csv` files are excluded from Git by `.gitignore`. Regenerate them locally by running `sgp4_orbit.py` in order.
