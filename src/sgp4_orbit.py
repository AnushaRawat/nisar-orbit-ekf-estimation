import requests
import numpy as np
import matplotlib.pyplot as plt
from sgp4.api import Satrec
from pathlib import Path
import pandas as pd
from pathlib import Path

data_dir = Path("data")
data_dir.mkdir(exist_ok=True)

EARTH_RADIUS_KM = 6378.137

def fetch_tle_from_celestrak(satellite_name="NISAR", group="sar"):
    url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP=sar&FORMAT=tle"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()

    lines = [l.strip() for l in resp.text.strip().splitlines() if l.strip()]

    # TLE list is 3 lines per satellite: name, line1, line2
    for i in range(0, len(lines) - 2, 3):
        name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
        if satellite_name.upper() in name.upper():
            return name, l1, l2

    raise ValueError(f"Satellite '{satellite_name}' not found in '{group}' group TLE list.")


name, line1, line2 = fetch_tle_from_celestrak("NISAR", group="sar")
print(name)
print(line1)
print(line2)

sat = Satrec.twoline2rv(line1, line2)

epoch_jd = sat.jdsatepoch
epoch_fr = sat.jdsatepochF
epoch_jd_full = epoch_jd + epoch_fr

# Convert Julian Date -> UTC datetime 
from datetime import datetime, timedelta
JD_UNIX_EPOCH = 2440587.5  # JD at 1970-01-01 00:00 UTC
unix_seconds = (epoch_jd_full - JD_UNIX_EPOCH) * 86400.0
epoch_utc = datetime(1970, 1, 1) + timedelta(seconds=unix_seconds)

print(f"TLE epoch (UTC): {epoch_utc}")
print(f"NORAD catalog number: {sat.satnum}")
print(f"Mean motion: {sat.no_kozai * (1440 / (2*np.pi)):.4f} rev/day")
print(f"Orbital period: {1440 / (sat.no_kozai * (1440 / (2*np.pi))):.2f} minutes")
print(f"Inclination: {np.degrees(sat.inclo):.2f} deg")

def propagate_orbit(line1, line2, days=3, step_minutes=10):
    sat = Satrec.twoline2rv(line1, line2)
    epoch_jd, epoch_fr = sat.jdsatepoch, sat.jdsatepochF

    n_steps = int(days * 24 * 60 / step_minutes) + 1
    times_min = np.arange(n_steps) * step_minutes

    positions = np.zeros((n_steps, 3))
    velocities = np.zeros((n_steps, 3))
    jds = np.zeros(n_steps)

    for i, t_min in enumerate(times_min):
        jd = epoch_jd
        fr = epoch_fr + (t_min / (24.0 * 60.0))
        error_code, r, v = sat.sgp4(jd, fr)
        if error_code != 0:
            raise RuntimeError(f"SGP4 propagation error {error_code} at t={t_min} min")
        positions[i] = r
        velocities[i] = v
        jds[i] = jd + fr

    return times_min, positions, velocities, jds


DAYS = 3
STEP_MINUTES = 10   

times_min, positions, velocities, jds = propagate_orbit(line1, line2, days=DAYS, step_minutes=STEP_MINUTES)
truth_states = np.hstack((positions, velocities))

# Save truth orbit
np.savez(
    "truth_orbit.npz",
    times=times_min,
    positions=positions,
    velocities=velocities,
    states=truth_states,
    julian_dates=jds
)
print(f"Saved {len(truth_states)} truth states.")
print(f"Propagated {len(times_min)} points over {DAYS} days at {STEP_MINUTES}-min cadence")

truth_states = np.hstack((positions, velocities))
np.savez(data_dir / "truth_orbit.npz",states=truth_states,positions=positions,velocities=velocities,times_minutes=times_min,julian_dates=jds)
truth_df = pd.DataFrame({"time_min": times_min,"jd": jds,"x_km": positions[:, 0],"y_km": positions[:, 1],"z_km": positions[:, 2],"vx_km_s": velocities[:, 0],"vy_km_s": velocities[:, 1],"vz_km_s": velocities[:, 2]})
truth_df.to_csv(data_dir / "truth_orbit.csv", index=False)
print(f"Saved {len(truth_df)} truth states to {data_dir/'truth_orbit.csv'}")

altitude_km = np.linalg.norm(positions, axis=1) - EARTH_RADIUS_KM
speed_km_s = np.linalg.norm(velocities, axis=1)
print(f"Altitude: min={altitude_km.min():.1f} km, max={altitude_km.max():.1f} km, mean={altitude_km.mean():.1f} km")
print(f"Speed:    min={speed_km_s.min():.3f} km/s, max={speed_km_s.max():.3f} km/s")

# Plotting the orbit in 3D
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
# Earth wireframe
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
ex = EARTH_RADIUS_KM * np.cos(u) * np.sin(v)
ey = EARTH_RADIUS_KM * np.sin(u) * np.sin(v)
ez = EARTH_RADIUS_KM * np.cos(v)
ax.plot_wireframe(ex, ey, ez, color="lightblue", alpha=0.3, linewidth=0.5)
# Orbit path
ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], color="crimson", linewidth=0.8)
ax.set_title(f"{name} orbit ({DAYS} days, TEME frame)")
ax.set_xlabel("X (km)"); ax.set_ylabel("Y (km)"); ax.set_zlabel("Z (km)")
ax.set_box_aspect([1, 1, 1])
plt.tight_layout()
plt.show()

#Adding noise to the truth states
#=========================================
data_dir = Path("data")
truth = np.load(data_dir / "truth_orbit.npz")

states = truth["states"]              # (N,6)
times = truth["times_minutes"]
jds = truth["julian_dates"]

np.random.seed(42)      

POSITION_STD = 0.1      # km (100 m)
VELOCITY_STD = 0.001    # km/s (1 m/s)

noise = np.zeros_like(states)

# Position noise
noise[:, 0:3] = np.random.normal(
    loc=0.0,
    scale=POSITION_STD,
    size=(states.shape[0], 3)
)

# Velocity noise
noise[:, 3:6] = np.random.normal(
    loc=0.0,
    scale=VELOCITY_STD,
    size=(states.shape[0], 3)
)

noisy_states = states + noise

np.savez(
    data_dir / "noisy_orbit.npz",
    states=noisy_states,
    times_minutes=times,
    julian_dates=jds
)

noisy_df = pd.DataFrame({
    "time_min": times,
    "jd": jds,
    "x_km": noisy_states[:,0],
    "y_km": noisy_states[:,1],
    "z_km": noisy_states[:,2],
    "vx_km_s": noisy_states[:,3],
    "vy_km_s": noisy_states[:,4],
    "vz_km_s": noisy_states[:,5],
})

noisy_df.to_csv(data_dir / "noisy_orbit.csv", index=False)

print(f"Saved {len(noisy_df)} noisy measurements.")
print(f"Position noise σ = {POSITION_STD:.3f} km")
print(f"Velocity noise σ = {VELOCITY_STD:.3f} km/s")

print("\nFirst truth state:")
print(states[0])

print("\nFirst noisy state:")
print(noisy_states[0])
