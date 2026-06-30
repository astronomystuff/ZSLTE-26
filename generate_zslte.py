import rebound
import numpy as np

# -----------------------------------------
# CONFIG
# -----------------------------------------
START_JD = 2451545.0          # J2000
YEARS = 10000000              # ±10 Myr
DT = 1.0                      # 1 day
STEPS = int(YEARS * 365.25 / DT)

# -----------------------------------------
# CREATE SIMULATION
# -----------------------------------------
sim = rebound.Simulation()
sim.units = ('AU', 'day', 'Msun')
sim.integrator = "whfast"
sim.dt = DT

# -----------------------------------------
# ADD SUN
# -----------------------------------------
sim.add(m=1.0)

# -----------------------------------------
# ADD PLANETS (J2000 barycentric)
# -----------------------------------------
# You will replace these with Horizons initial states
planets = [
    ("Mercury", 1.651e-7),
    ("Venus",   2.447e-6),
    ("Earth",   3.003e-6),
    ("Mars",    3.227e-7),
    ("Jupiter", 9.545e-4),
    ("Saturn",  2.857e-4),
    ("Uranus",  4.365e-5),
    ("Neptune", 5.150e-5)
]

for name, mass in planets:
    sim.add(m=mass)  # placeholder; you will add full state vectors

# -----------------------------------------
# ADD BIG ASTEROIDS
# -----------------------------------------
asteroids = [
    ("Ceres",      4.7e-10),
    ("Vesta",      1.3e-10),
    ("Pallas",     1.1e-10),
    ("Hygiea",     8.6e-11),
    ("Interamnia", 3.5e-11)
]

for name, mass in asteroids:
    sim.add(m=mass)  # placeholder; you will add full state vectors

# -----------------------------------------
# MOVE TO BARYCENTER
# -----------------------------------------
sim.move_to_com()

# -----------------------------------------
# ENABLE 1PN RELATIVITY (EIH)
# -----------------------------------------
sim.add_force("gr_pn")

# -----------------------------------------
# OUTPUT STORAGE
# -----------------------------------------
out = {
    "jd": [],
    "x": [],
    "y": [],
    "z": []
}

# -----------------------------------------
# MAIN LOOP
# -----------------------------------------
current_jd = START_JD

for i in range(STEPS):
    sim.step(DT)
    current_jd += DT

    earth = sim.particles[3]  # example: Earth
    out["jd"].append(current_jd)
    out["x"].append(earth.x)
    out["y"].append(earth.y)
    out["z"].append(earth.z)

# -----------------------------------------
# SAVE OUTPUT
# -----------------------------------------
np.savez("zslte_raw_xyz.npz", **out)
print("Done.")
