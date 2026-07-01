import numpy as np
import rebound
import json

G = 6.67430e-11
c = 299792458.0
c2 = c*c

# ------------------------------------------------------------
# 1PN GR correction (Sun-only Schwarzschild)
# ------------------------------------------------------------
def gr_accel_sun_only(sim):
    ps = sim.particles
    sun = ps[0]
    GM = G * sun.m

    sx, sy, sz = sun.x, sun.y, sun.z
    svx, svy, svz = sun.vx, sun.vy, sun.vz

    for i in range(1, len(ps)):
        p = ps[i]

        rx = p.x - sx
        ry = p.y - sy
        rz = p.z - sz
        r2 = rx*rx + ry*ry + rz*rz
        if r2 == 0.0:
            continue

        r = np.sqrt(r2)
        vx = p.vx - svx
        vy = p.vy - svy
        vz = p.vz - svz
        v2 = vx*vx + vy*vy + vz*vz
        rv = rx*vx + ry*vy + rz*vz

        inv_r3 = 1.0 / (r2 * r)
        fac0 = GM / c2 * inv_r3
        term = 4.0 * GM / r - v2

        p.ax += fac0 * (term * rx + 4.0 * rv * vx)
        p.ay += fac0 * (term * ry + 4.0 * rv * vy)
        p.az += fac0 * (term * rz + 4.0 * rv * vz)


# ------------------------------------------------------------
# Chebyshev helpers
# ------------------------------------------------------------
def cheb_nodes(t0, t1, n):
    k = np.arange(n, dtype=np.float64)
    theta = np.pi * (2*k + 1) / (2*n)
    tau = np.cos(theta)
    return 0.5*(t0 + t1) + 0.5*(t1 - t0)*tau

def cheb_fit(values, degree):
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    k = np.arange(n, dtype=np.float64)
    theta = np.pi * (2*k + 1) / (2*n)
    coeffs = np.zeros(degree+1, dtype=np.float64)
    for m in range(degree+1):
        coeffs[m] = (2.0/n) * np.sum(values * np.cos(m * theta))
    coeffs[0] *= 0.5
    return coeffs

# ------------------------------------------------------------
# J2000 initial conditions (barycentric, DE430/DE440)
# Units: meters, meters/second
# ------------------------------------------------------------
def get_initial_bodies_j2000():
    return [

        # Sun (barycenter offset is tiny; DE sets Sun near origin)
        {'name':'Sun','m':1.98847e30,
         'x':-1.068638e6,'y':-4.195e5,'z':1.880e4,
         'vx':0.007,'vy':-0.011,'vz':-0.0002},

        # Mercury
        {'name':'Mercury','m':3.3011e23,
         'x':-3.713e10,'y':-5.569e10,'z':-2.727e9,
         'vx':2.986e4,'vy':-1.908e4,'vz':-4.889e3},

        # Venus
        {'name':'Venus','m':4.8675e24,
         'x':-6.950e10,'y':-2.273e10,'z':4.086e9,
         'vx':7.783e3,'vy':-3.453e4,'vz':-1.573e3},

        # Earth
        {'name':'Earth','m':5.9722e24,
         'x':-2.632e10,'y':1.445e11,'z':-1.258e7,
         'vx':-2.978e4,'vy':-5.453e3,'vz':0.0007},

        # Mars
        {'name':'Mars','m':6.4171e23,
         'x':1.876e11,'y':-1.012e11,'z':-6.725e9,
         'vx':1.020e4,'vy':2.403e4,'vz':1.017e3},

        # Jupiter
        {'name':'Jupiter','m':1.89813e27,
         'x':3.815e11,'y':3.951e11,'z':-1.015e10,
         'vx':-7.350e3,'vy':6.297e3,'vz':1.667e2},

        # Saturn
        {'name':'Saturn','m':5.6834e26,
         'x':6.503e11,'y':-3.747e11,'z':-1.046e10,
         'vx':4.665e3,'vy':8.987e3,'vz':-2.273e2},

        # Uranus
        {'name':'Uranus','m':8.6810e25,
         'x':1.434e12,'y':1.754e11,'z':-1.508e10,
         'vx':-3.725e3,'vy':6.046e3,'vz':6.727e1},

        # Neptune
        {'name':'Neptune','m':1.02413e26,
         'x':1.683e12,'y':-2.989e11,'z':-1.781e10,
         'vx':3.919e3,'vy':5.393e3,'vz':-1.511e2},

        # Pluto
        {'name':'Pluto','m':1.303e22,
         'x':-1.136e12,'y':-4.953e11,'z':2.093e10,
         'vx':1.150e4,'vy':-1.080e4,'vz':-2.748e2},

        # Ceres
        {'name':'Ceres','m':9.393e20,
         'x':-2.707e11,'y':1.956e11,'z':1.322e10,
         'vx':-1.073e4,'vy':-1.531e4,'vz':-1.046e3},

        # Vesta
        {'name':'Vesta','m':2.59076e20,
         'x':-1.353e11,'y':-3.545e11,'z':-1.725e10,
         'vx':1.587e4,'vy':-6.291e3,'vz':-1.066e3},

        # Pallas
        {'name':'Pallas','m':2.11e20,
         'x':3.195e11,'y':-1.689e11,'z':-2.118e10,
         'vx':-1.046e4,'vy':1.581e4,'vz':1.092e3},

        # Hygiea
        {'name':'Hygiea','m':8.67e19,
         'x':-2.998e11,'y':-1.238e11,'z':-1.497e10,
         'vx':1.287e4,'vy':-1.853e4,'vz':-1.011e3},

        # Interamnia
        {'name':'Interamnia','m':3.9e19,
         'x':-3.548e11,'y':-1.044e11,'z':-1.221e10,
         'vx':1.178e4,'vy':-1.944e4,'vz':-9.88e2},
    ]

# ------------------------------------------------------------
# Ephemeris generator
# ------------------------------------------------------------
def generate_ephemeris(bodies, dt, steps, block_size, cheb_degree):
    sim = rebound.Simulation()
    sim.units = ('m','s','kg')
    sim.integrator = "ias15"
    sim.dt = dt

    for b in bodies:
        sim.add(m=b['m'], x=b['x'], y=b['y'], z=b['z'],
                vx=b['vx'], vy=b['vy'], vz=b['vz'])

    sim.move_to_com()
    sim.additional_forces = gr_accel_sun_only

    N = len(bodies)
    result = {'bodies':[{'name':b['name'],'segments':[]} for b in bodies]}
    current_step = 0

    while current_step < steps:
        block_steps = min(block_size, steps - current_step)
        t_start = sim.t
        t_end = t_start + block_steps * dt

        n_nodes = cheb_degree + 1
        node_times = cheb_nodes(t_start, t_end, n_nodes)

        samples_x = [np.zeros(n_nodes) for _ in range(N)]
        samples_y = [np.zeros(n_nodes) for _ in range(N)]
        samples_z = [np.zeros(n_nodes) for _ in range(N)]

        for k, t_node in enumerate(node_times):
            sim.integrate(t_node)
            for i in range(N):
                p = sim.particles[i]
                samples_x[i][k] = p.x
                samples_y[i][k] = p.y
                samples_z[i][k] = p.z

        current_step += block_steps

        for i in range(N):
            cx = cheb_fit(samples_x[i], cheb_degree)
            cy = cheb_fit(samples_y[i], cheb_degree)
            cz = cheb_fit(samples_z[i], cheb_degree)
            result['bodies'][i]['segments'].append({
                't0':float(t_start),
                't1':float(t_end),
                'coeffsX':cx.tolist(),
                'coeffsY':cy.tolist(),
                'coeffsZ':cz.tolist()
            })

    return result

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    bodies = get_initial_bodies_j2000()

    dt = 86400.0          # 1 day
    steps = 365 * 50      # 50 years
    block_size = 73       # Chebyshev segment size
    cheb_degree = 12      # polynomial degree

    eph = generate_ephemeris(bodies, dt, steps, block_size, cheb_degree)

    with open("zs_ephemeris.json","w") as f:
        json.dump(eph,f,indent=2)

if __name__ == "__main__":
    main()
