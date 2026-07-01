// cartesian_wh_ephemeris.js
function generateCartesianWHEphemeris(config) {
    // config = {
    //   bodies: [ { name, m, x, y, z, vx, vy, vz }, ... ], // Sun first
    //   dt: timestep_seconds,          // e.g. 86400 / 4
    //   steps: total_steps,            // total integration steps
    //   central: 0,                    // index of Sun
    //   blockSize: steps_per_block,    // steps per Chebyshev segment
    //   chebDegree: polynomial_degree  // e.g. 12–20
    // }

    const G  = 6.67430e-11;
    const c2 = 299792458 ** 2;

    const N = config.bodies.length;
    const central = config.central;
    const dt = config.dt;

    // ---------- STATE ----------
    let state = config.bodies.map(b => ({
        name: b.name,
        m:    b.m,
        x:    b.x,
        y:    b.y,
        z:    b.z,
        vx:   b.vx,
        vy:   b.vy,
        vz:   b.vz
    }));

    // ---------- STUMPFF FUNCTIONS ----------
    function stumpffC(z) {
        if (Math.abs(z) < 1e-8) return 0.5 - z/24 + z*z/720;
        if (z > 0) {
            const s = Math.sqrt(z);
            return (1 - Math.cos(s)) / z;
        } else {
            const s = Math.sqrt(-z);
            return (1 - Math.cosh(s)) / z;
        }
    }

    function stumpffS(z) {
        if (Math.abs(z) < 1e-8) return 1/6 - z/120 + z*z/5040;
        if (z > 0) {
            const s = Math.sqrt(z);
            return (s - Math.sin(s)) / (s*z);
        } else {
            const s = Math.sqrt(-z);
            return (Math.sinh(s) - s) / (s*z);
        }
    }

    // ---------- KEPLER DRIFT (universal variables, heliocentric) ----------
    function keplerDriftUV(state, i, h) {
        if (i === central) return;

        const sun = state[central];
        const bi  = state[i];

        const M0 = sun.m;
        const GM = G * M0;

        // heliocentric position/velocity
        let rx = bi.x - sun.x;
        let ry = bi.y - sun.y;
        let rz = bi.z - sun.z;
        let vx = bi.vx - sun.vx;
        let vy = bi.vy - sun.vy;
        let vz = bi.vz - sun.vz;

        const r0    = Math.sqrt(rx*rx + ry*ry + rz*rz);
        const v0sq  = vx*vx + vy*vy + vz*vz;
        const vr0   = (rx*vx + ry*vy + rz*vz) / r0;
        const alpha = 2*GM/r0 - v0sq;

        const sqrtGM = Math.sqrt(GM);
        let X = sqrtGM * h * alpha;

        for (let iter = 0; iter < 30; iter++) {
            const z = alpha * X * X;
            const C = stumpffC(z);
            const S = stumpffS(z);

            const F  = r0*vr0/sqrtGM * X*X*C + (1 - alpha*r0)*X*X*X*S + r0*X - sqrtGM*h;
            const dF = r0*vr0/sqrtGM * X*(1 - z*S) + (1 - alpha*r0)*X*X*C + r0;

            const dX = F/dF;
            X -= dX;
            if (Math.abs(dX) < 1e-13) break;
        }

        const z = alpha * X * X;
        const C = stumpffC(z);
        const S = stumpffS(z);

        const f = 1 - (X*X*C)/r0;
        const g = h - (X*X*X*S)/sqrtGM;

        const rx_new = f*rx + g*vx;
        const ry_new = f*ry + g*vy;
        const rz_new = f*rz + g*vz;

        const r_new = Math.sqrt(rx_new*rx_new + ry_new*ry_new + rz_new*rz_new);

        const fdot = (sqrtGM/(r0*r_new)) * (z*S - 1);
        const gdot = 1 - (X*X*C)/r_new;

        const vx_new = fdot*rx + gdot*vx;
        const vy_new = fdot*ry + gdot*vy;
        const vz_new = fdot*rz + gdot*vz;

        // back to barycentric
        bi.x  = sun.x + rx_new;
        bi.y  = sun.y + ry_new;
        bi.z  = sun.z + rz_new;
        bi.vx = sun.vx + vx_new;
        bi.vy = sun.vy + vy_new;
        bi.vz = sun.vz + vz_new;
    }

    // ---------- PERTURBATION KICK (planet–planet + GR, Cartesian) ----------
    const ax = new Float64Array(N);
    const ay = new Float64Array(N);
    const az = new Float64Array(N);

    function computePerturbations(state) {
        for (let i = 0; i < N; i++) {
            ax[i] = 0; ay[i] = 0; az[i] = 0;
        }

        // planet–planet Newtonian (Sun–planet handled by Kepler drift)
        for (let i = 0; i < N; i++) {
            if (i === central) continue;
            const bi = state[i];
            for (let j = 0; j < N; j++) {
                if (j === central || j === i) continue;
                const bj = state[j];

                const dx = bj.x - bi.x;
                const dy = bj.y - bi.y;
                const dz = bj.z - bi.z;
                const r2 = dx*dx + dy*dy + dz*dz;
                if (r2 === 0) continue;

                const invR  = 1 / Math.sqrt(r2);
                const invR3 = invR * invR * invR;

                const fac = G * bj.m * invR3;
                ax[i] += fac * dx;
                ay[i] += fac * dy;
                az[i] += fac * dz;
            }
        }

        // GR 1PN relative to Sun, heliocentric
        const sun = state[central];
        const M0  = sun.m;
        const GM  = G * M0;

        for (let i = 0; i < N; i++) {
            if (i === central) continue;

            const bi = state[i];

            const rx = bi.x - sun.x;
            const ry = bi.y - sun.y;
            const rz = bi.z - sun.z;
            const r2 = rx*rx + ry*ry + rz*rz;
            if (r2 === 0) continue;

            const invR  = 1 / Math.sqrt(r2);
            const invR3 = invR * invR * invR;
            const r     = 1 / invR;

            const vx = bi.vx - sun.vx;
            const vy = bi.vy - sun.vy;
            const vz = bi.vz - sun.vz;
            const v2 = vx*vx + vy*vy + vz*vz;
            const rv = rx*vx + ry*vy + rz*vz;

            const fac0 = GM / c2 * invR3;
            const term = 4 * GM / r - v2;

            ax[i] += fac0 * ( term * rx + 4 * rv * vx );
            ay[i] += fac0 * ( term * ry + 4 * rv * vy );
            az[i] += fac0 * ( term * rz + 4 * rv * vz );
        }
    }

    function kick(state, h) {
        computePerturbations(state);
        for (let i = 0; i < N; i++) {
            state[i].vx += 0.5 * h * ax[i];
            state[i].vy += 0.5 * h * ay[i];
            state[i].vz += 0.5 * h * az[i];
        }
    }

    // ---------- SIMPLE SYMPLECTIC CORRECTOR HOOK ----------
    // For now this is a no-op; you can later implement a WH-style corrector here.
    function applyCorrector(state, h) {
        // placeholder for symplectic corrector
    }

    // ---------- CHEBYSHEV NODES + TRANSFORM ----------
    function chebNodes(t0, t1, n) {
        const times = new Float64Array(n);
        for (let k = 0; k < n; k++) {
            const theta = Math.PI * (2*k + 1) / (2*n);
            const tau   = Math.cos(theta);
            times[k] = 0.5*(t0 + t1) + 0.5*(t1 - t0)*tau;
        }
        return times; // Chebyshev order
    }

    function chebFitFromSamples(values, degree) {
        const n = values.length;
        const coeffs = new Float64Array(degree + 1);
        for (let m = 0; m <= degree; m++) {
            let sum = 0;
            for (let k = 0; k < n; k++) {
                const theta = Math.PI * (2*k + 1) / (2*n);
                sum += values[k] * Math.cos(m * theta);
            }
            coeffs[m] = (2 / n) * sum;
        }
        coeffs[0] *= 0.5;
        return coeffs;
    }

    // ---------- MAIN LOOP ----------
    const result = {
        bodies: state.map(b => ({ name: b.name, segments: [] }))
    };

    let globalStep = 0;

    while (globalStep < config.steps) {
        const blockSteps = Math.min(config.blockSize, config.steps - globalStep);
        const tStart = globalStep * dt;
        const tEnd   = (globalStep + blockSteps) * dt;

        // block-start snapshot
        const blockStartState = state.map(b => ({ ...b }));

        // apply corrector at block start (currently no-op)
        applyCorrector(state, dt);

        // integrate global state through block: KDK
        for (let s = 0; s < blockSteps; s++) {
            kick(state, dt);
            for (let i = 0; i < N; i++) {
                keplerDriftUV(state, i, dt);
            }
            kick(state, dt);

            globalStep++;
        }

        // apply inverse corrector at block end (currently no-op)
        applyCorrector(state, -dt);

        // Chebyshev sampling from block-start
        const nNodes = config.chebDegree + 1;
        const nodeTimes = chebNodes(tStart, tEnd, nNodes);

        const samplesX = Array(N).fill(null).map(() => new Float64Array(nNodes));
        const samplesY = Array(N).fill(null).map(() => new Float64Array(nNodes));
        const samplesZ = Array(N).fill(null).map(() => new Float64Array(nNodes));

        let localState = blockStartState.map(b => ({ ...b }));
        let localTime  = tStart;

        for (let k = 0; k < nNodes; k++) {
            const targetT = nodeTimes[k];

            while (localTime < targetT) {
                const h = Math.min(dt, targetT - localTime);

                kick(localState, h);
                for (let i = 0; i < N; i++) {
                    keplerDriftUV(localState, i, h);
                }
                kick(localState, h);

                localTime += h;
            }

            for (let i = 0; i < N; i++) {
                samplesX[i][k] = localState[i].x;
                samplesY[i][k] = localState[i].y;
                samplesZ[i][k] = localState[i].z;
            }
        }

        // fit Chebyshev per body per coordinate
        for (let i = 0; i < N; i++) {
            const coeffsX = chebFitFromSamples(samplesX[i], config.chebDegree);
            const coeffsY = chebFitFromSamples(samplesY[i], config.chebDegree);
            const coeffsZ = chebFitFromSamples(samplesZ[i], config.chebDegree);

            result.bodies[i].segments.push({
                t0: tStart,
                t1: tEnd,
                coeffsX: Array.from(coeffsX),
                coeffsY: Array.from(coeffsY),
                coeffsZ: Array.from(coeffsZ)
            });
        }
    }

    return result;
}
