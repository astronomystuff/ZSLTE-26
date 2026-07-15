class solarSystemEngine {
    constructor() {
        this.gConstant = 2.959122082855911e-4; // AU^3 / (day^2 * M_sun)
        this.cLight = 173.144543265; // AU / day
        this.cSquare = this.cLight * this.cLight;
        this.bodies = [];
    }

    addBody(name, mass, x, y, z, vx, vy, vz, isPerturberOnly = false) {
        this.bodies.push({
            name, mass, isPerturberOnly,
            pos: new Float64Array([x, y, z]),
            vel: new Float64Array([vx, vy, vz])
        });
    }

    addBodyFromElements(name, mass, semiMajorAxis, eccentricity, inclinationDeg, omegaBigDeg, omegaSmallDeg, meanAnomalyDeg, meanMotionDeg, nativeEpochJD, isPerturberOnly = false) {
        const simStartEpochJD = 2451545.0; // J2000
        
        const deltaDays = simStartEpochJD - nativeEpochJD;
        let correctedMeanAnomalyDeg = meanAnomalyDeg + (meanMotionDeg * deltaDays);
        
        correctedMeanAnomalyDeg = correctedMeanAnomalyDeg % 360;
        if (correctedMeanAnomalyDeg < 0) {
            correctedMeanAnomalyDeg += 360;
        }

        const inclinationRad = inclinationDeg * Math.PI / 180;
        const omegaBigRad = omegaBigDeg * Math.PI / 180;
        const omegaSmallRad = omegaSmallDeg * Math.PI / 180;
        const meanAnomalyRad = correctedMeanAnomalyDeg * Math.PI / 180;
        const meanMotionRad = meanMotionDeg * Math.PI / 180;

        let eccentricAnomaly = meanAnomalyRad;
        for (let iteration = 0; iteration < 100; iteration++) {
            const deltaE = (eccentricAnomaly - eccentricity * Math.sin(eccentricAnomaly) - meanAnomalyRad) / (1 - eccentricity * Math.cos(eccentricAnomaly));
            eccentricAnomaly -= deltaE;
            if (Math.abs(deltaE) < 1e-12) {
                break;
            }
        }

        const xPlane = semiMajorAxis * (Math.cos(eccentricAnomaly) - eccentricity);
        const yPlane = semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * Math.sin(eccentricAnomaly);
        
        const rDistance = semiMajorAxis * (1 - eccentricity * Math.cos(eccentricAnomaly));
        const vFactor = (semiMajorAxis * semiMajorAxis * meanMotionRad) / rDistance;
        
        const vxPlane = -vFactor * Math.sin(eccentricAnomaly);
        const vyPlane = vFactor * Math.sqrt(1 - eccentricity * eccentricity) * Math.cos(eccentricAnomaly);

        const cosOmegaBig = Math.cos(omegaBigRad);
        const sinOmegaBig = Math.sin(omegaBigRad);
        const cosOmegaSmall = Math.cos(omegaSmallRad);
        const sinOmegaSmall = Math.sin(omegaSmallRad);
        const cosInclination = Math.cos(inclinationRad);
        const sinInclination = Math.sin(inclinationRad);

        const pX = cosOmegaSmall * cosOmegaBig - sinOmegaSmall * sinOmegaBig * cosInclination;
        const pY = cosOmegaSmall * sinOmegaBig + sinOmegaSmall * cosOmegaBig * cosInclination;
        const pZ = sinOmegaSmall * sinInclination;

        const qX = -sinOmegaSmall * cosOmegaBig - cosOmegaSmall * sinOmegaBig * cosInclination;
        const qY = -sinOmegaSmall * sinOmegaBig + cosOmegaSmall * cosOmegaBig * cosInclination;
        const qZ = cosOmegaSmall * sinInclination;

        const xPos = xPlane * pX + yPlane * qX;
        const yPos = xPlane * pY + yPlane * qY;
        const zPos = xPlane * pZ + yPlane * qZ;

        const vxVel = vxPlane * pX + vyPlane * qX;
        const vyVel = vxPlane * pY + vyPlane * qY;
        const vzVel = vxPlane * pZ + vyPlane * qZ;

        this.addBody(name, mass, xPos, yPos, zPos, vxVel, vyVel, vzVel, isPerturberOnly);
    }

    balanceBarycenter() {
        let totalMass = 0;
        const centerPos = new Float64Array([0, 0, 0]);
        const centerVel = new Float64Array([0, 0, 0]);

        this.bodies.forEach(body => {
            totalMass += body.mass;
            centerPos[0] += body.pos[0] * body.mass;
            centerPos[1] += body.pos[1] * body.mass;
            centerPos[2] += body.pos[2] * body.mass;
            centerVel[0] += body.vel[0] * body.mass;
            centerVel[1] += body.vel[1] * body.mass;
            centerVel[2] += body.vel[2] * body.mass;
        });

        this.bodies.forEach(body => {
            body.pos[0] -= centerPos[0] / totalMass;
            body.pos[1] -= centerPos[1] / totalMass;
            body.pos[2] -= centerPos[2] / totalMass;
            body.vel[0] -= centerVel[0] / totalMass;
            body.vel[1] -= centerVel[1] / totalMass;
            body.vel[2] -= centerVel[2] / totalMass;
        });
    }

    // Universal Variable Kepler Drift: Resolves the e=0 singularity and supports e>=1 analytically
    driftKeplerUniversal(body, sun, dt) {
        const mu = this.gConstant * (sun.mass + body.mass);
        
        // Relative state vectors
        const rx = body.pos[0] - sun.pos[0];
        const ry = body.pos[1] - sun.pos[1];
        const rz = body.pos[2] - sun.pos[2];
        const r0 = Math.sqrt(rx*rx + ry*ry + rz*rz);

        const vx = body.vel[0] - sun.vel[0];
        const vy = body.vel[1] - sun.vel[1];
        const vz = body.vel[2] - sun.vel[2];
        
        const v2 = vx*vx + vy*vy + vz*vz;
        const rDotV = rx*vx + ry*vy + rz*vz;

        // Reciprocal of semi-major axis (alpha)
        const alpha = 2.0 / r0 - v2 / mu;

        // Solve Kepler's equation using Universal Variables (Danby's method)
        let chi = Math.sqrt(mu) * Math.abs(alpha) * dt; // Initial guess
        let chiNew = chi;
        const tolerance = 1e-14;
        let c0 = 0, c1 = 0, c2 = 0, c3 = 0;

        for (let iter = 0; iter < 100; iter++) {
            const psi = chi * chi * alpha;
            
            // Compute Stumpff functions c2 and c3
            if (psi > 1e-5) {
                const sqrtPsi = Math.sqrt(psi);
                c2 = (1.0 - Math.cos(sqrtPsi)) / psi;
                c3 = (sqrtPsi - Math.sin(sqrtPsi)) / (psi * sqrtPsi);
            } else if (psi < -1e-5) {
                const sqrtNegPsi = Math.sqrt(-psi);
                c2 = (1.0 - Math.cosh(sqrtNegPsi)) / psi;
                c3 = (Math.sinh(sqrtNegPsi) - sqrtNegPsi) / (psi * sqrtNegPsi);
            } else {
                c2 = 1.0 / 2.0 - psi / 24.0 + (psi * psi) / 720.0;
                c3 = 1.0 / 6.0 - psi / 120.0 + (psi * psi) / 5040.0;
            }

            const f = r0 * chi * (1.0 - psi * c3) + rDotV * chi * chi * c2 + Math.sqrt(mu) * chi * chi * chi * c3 - Math.sqrt(mu) * dt;
            const fPrime = r0 * (1.0 - psi * c2) + rDotV * chi * (1.0 - psi * c3) + Math.sqrt(mu) * chi * chi * c2;
            
            const delta = f / fPrime;
            chiNew = chi - delta;
            
            if (Math.abs(delta) < tolerance) {
                chi = chiNew;
                break;
            }
            chi = chiNew;
        }

        const psi = chi * chi * alpha;
        const f = 1.0 - (chi * chi / r0) * c2;
        const g = dt - (chi * chi * chi / Math.sqrt(mu)) * c3;

        const xNew = f * rx + g * vx;
        const yNew = f * ry + g * vy;
        const zNew = f * rz + g * vz;
        const rNew = Math.sqrt(xNew*xNew + yNew*yNew + zNew*zNew);

        const fDot = (Math.sqrt(mu) / (rNew * r0)) * chi * (psi * c3 - 1.0);
        const gDot = 1.0 - (chi * chi / rNew) * c2;

        // Absolute state update relative to Sun's motion
        body.pos[0] = sun.pos[0] + xNew;
        body.pos[1] = sun.pos[1] + yNew;
        body.pos[2] = sun.pos[2] + zNew;

        body.vel[0] = sun.vel[0] + fDot * rx + gDot * vx;
        body.vel[1] = sun.vel[1] + fDot * ry + gDot * vy;
        body.vel[2] = sun.vel[2] + fDot * rz + gDot * vz;
    }

    // Complete Einstein-Infeld-Hoffmann (EIH) N-Body Equations of Motion
    getEIHAccelerations() {
        const n = this.bodies.length;
        const acc = Array.from({ length: n }, () => new Float64Array([0, 0, 0]));
        const G = this.gConstant;
        const c2 = this.cSquare;

        // Precompute Cartesian separation vectors and magnitudes
        const r_vec = Array.from({ length: n }, () => Array(n));
        const r_mag = Array.from({ length: n }, () => new Float64Array(n));
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (i === j) continue;
                const dx = this.bodies[j].pos[0] - this.bodies[i].pos[0];
                const dy = this.bodies[j].pos[1] - this.bodies[i].pos[1];
                const dz = this.bodies[j].pos[2] - this.bodies[i].pos[2];
                r_vec[i][j] = new Float64Array([dx, dy, dz]);
                r_mag[i][j] = Math.sqrt(dx*dx + dy*dy + dz*dz);
            }
        }

        // Apply EIH equations of motion symmetrically across all active bodies
        for (let i = 0; i < n; i++) {
            const mi = this.bodies[i].mass;
            const vi = this.bodies[i].vel;
            const vi2 = vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2];

            for (let j = 0; j < n; j++) {
                if (i === j) continue;

                const mj = this.bodies[j].mass;
                const vj = this.bodies[j].vel;
                const vj2 = vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2];
                
                const rij_vec = r_vec[i][j];
                const rij = r_mag[i][j];
                const rij3 = rij * rij * rij;

                // Classical Mutual Newtonian baseline (no Sun exclusion)
                const newtonianFactor = G * mj / rij3;
                let ax_N = rij_vec[0] * newtonianFactor;
                let ay_N = rij_vec[1] * newtonianFactor;
                let az_N = rij_vec[2] * newtonianFactor;

                // EIH Relativistic Terms (includes three-body potentials)
                const vij_dot_rij = (vi[0] - vj[0])*rij_vec[0] + (vi[1] - vj[1])*rij_vec[1] + (vi[2] - vj[2])*rij_vec[2];
                
                let sum_potential_i = 0;
                let sum_three_body_k = 0;

                for (let k = 0; k < n; k++) {
                    if (k !== i) sum_potential_i += G * this.bodies[k].mass / r_mag[i][k];
                    if (k !== j) {
                        const rjk_vec = r_vec[j][k];
                        const rjk = r_mag[j][k];
                        const rij_dot_rjk = rij_vec[0]*rjk_vec[0] + rij_vec[1]*rjk_vec[1] + rij_vec[2]*rjk_vec[2];
                        sum_three_body_k += (G * this.bodies[k].mass / rjk) * (1.0 - rij_dot_rjk / (2.0 * rjk * rjk));
                    }
                }

                const scalarPN = 1.0 - (4.0 / c2) * sum_potential_i - (1.0 / c2) * sum_three_body_k +
                                 (1.0 / c2) * (vi2 + 2.0*vj2 - 4.0*(vi[0]*vj[0] + vi[1]*vj[1] + vi[2]*vj[2])) -
                                 (1.5 / c2) * (vij_dot_rij * vij_dot_rij / (rij * rij));

                const vectorFactor = (1.0 / c2) * (4.0 * (rij_vec[0]*vi[0] + rij_vec[1]*vi[1] + rij_vec[2]*vi[2]) - 
                                                  3.0 * (rij_vec[0]*vj[0] + rij_vec[1]*vj[1] + rij_vec[2]*vj[2]));

                acc[i][0] += ax_N * scalarPN + (G * mj / rij3) * vectorFactor * (vi[0] - vj[0]);
                acc[i][1] += ay_N * scalarPN + (G * mj / rij3) * vectorFactor * (vi[1] - vj[1]);
                acc[i][2] += az_N * scalarPN + (G * mj / rij3) * vectorFactor * (vi[2] - vj[2]);
            }
        }
        return acc;
    }

    // Symmetric Drift-Kick-Drift execution sequence preserving complete linear momentum
    step(dt) {
        const n = this.bodies.length;
        const sun = this.bodies[0];

        // 1. Kick (Half-step velocity update to both Sun and planets)
        const accHalf1 = this.getEIHAccelerations();
        for (let i = 0; i < n; i++) {
            this.bodies[i].vel[0] += accHalf1[i][0] * (dt * 0.5);
            this.bodies[i].vel[1] += accHalf1[i][1] * (dt * 0.5);
            this.bodies[i].vel[2] += accHalf1[i][2] * (dt * 0.5);
        }

        // 2. Full Drift (Barycentric linear Sun drift + Universal relative Keplerian drift)
        const sunDx = sun.vel[0] * dt;
        const sunDy = sun.vel[1] * dt;
        const sunDz = sun.vel[2] * dt;

        for (let i = 1; i < n; i++) {
            this.driftKeplerUniversal(this.bodies[i], sun, dt);
        }

        sun.pos[0] += sunDx;
        sun.pos[1] += sunDy;
        sun.pos[2] += sunDz;

        // 3. Kick (Half-step velocity update to both Sun and planets)
        const accHalf2 = this.getEIHAccelerations();
        for (let i = 0; i < n; i++) {
            this.bodies[i].vel[0] += accHalf2[i][0] * (dt * 0.5);
            this.bodies[i].vel[1] += accHalf2[i][1] * (dt * 0.5);
            this.bodies[i].vel[2] += accHalf2[i][2] * (dt * 0.5);
        }
    }

    downloadEphemeris(historyData, filename) {
        const dataString = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(historyData, null, 2));
        const downloadAnchor = document.createElement('a');
        downloadAnchor.setAttribute("href", dataString);
        downloadAnchor.setAttribute("download", filename);
        document.body.appendChild(downloadAnchor);
        downloadAnchor.click();
        downloadAnchor.remove();
        console.log(`Saved successfully as: ${filename}`);
    }

    // Integrates based on true continuous days (avoids floating-point year accumulator drift)
    runSimulation(totalYears, deltaTimeDays, snapshotIntervalYears, outputFilename = "ZSLTE-26-ephemeris.json") {
        this.balanceBarycenter(); 

        const totalSimDays = totalYears * 365.25002; 
        let elapsedDays = 0;
        let nextSnapshotDays = 0;
        const snapshotStepDays = snapshotIntervalYears * 365.25002;

        const historyData = {};
        this.bodies.forEach(body => { 
            if (!body.isPerturberOnly) {
                historyData[body.name] = []; 
            }
        });

        console.log(`Simulating ${totalYears.toLocaleString()} years via Cartesian WHFast-style EIH Engine...`);
        
        while (elapsedDays <= totalSimDays) {
            this.step(deltaTimeDays);
            elapsedDays += deltaTimeDays;

            if (elapsedDays >= nextSnapshotDays) {
                const sun = this.bodies[0];
                const currentSimYear = Math.round(elapsedDays / 365.25002);
                
                this.bodies.forEach(body => {
                    if (body.isPerturberOnly) return;
                    historyData[body.name].push({
                        t: currentSimYear,
                        x: Number((body.pos[0] - sun.pos[0]).toFixed(8)), 
                        y: Number((body.pos[1] - sun.pos[1]).toFixed(8)),
                        z: Number((body.pos[2] - sun.pos[2]).toFixed(8))
                    });
                });
                
                nextSnapshotDays += snapshotStepDays;
                console.log(`Progress: ${((elapsedDays / totalSimDays) * 100).toFixed(2)}%`);
            }
        }
        this.downloadEphemeris(historyData, outputFilename);
    }
}

// --- INITIALIZE TIMELINE INJECTION ---
const sim = new solarSystemEngine();
const J2000 = 2451545.0;
const JD2026 = 2461000.5;

// Sun Anchor (S - Exported)
sim.addBody("Sun", 1.0, 0, 0, 0, 0, 0, 0, false);

// 1. Planets (J2000 Epoch - Exported)
sim.addBodyFromElements("Mercury", 1.66013e-7,  0.387, 0.2056, 7.00, 48.33,  29.12,  174.8, 4.092,  J2000, false); 
sim.addBodyFromElements("Venus",   2.447838e-6, 0.723, 0.0068, 3.39, 76.68,  54.88,  50.1,  1.602,  J2000, false); 
sim.addBodyFromElements("Earth",   3.003489e-6, 1.000, 0.0167, 0.00, 0.00,   102.94, 357.5, 0.9856, J2000, false); 
sim.addBodyFromElements("Mars",    3.227151e-7, 1.524, 0.0934, 1.85, 49.56,  286.50, 19.4,  0.524,  J2000, false); 
sim.addBodyFromElements("Jupiter", 9.547919e-4, 5.203, 0.0489, 1.30, 100.46, 273.87, 20.0,  0.0831, J2000, false); 
sim.addBodyFromElements("Saturn",  2.85886e-4,  9.537, 0.0565, 2.49, 113.72, 339.39, 317.0, 0.0335, J2000, false); 
sim.addBodyFromElements("Uranus",  4.366244e-5, 19.191, 0.0472, 0.77, 74.00,  96.88,  142.2, 0.0117, J2000, false); 
sim.addBodyFromElements("Neptune", 5.151389e-5, 30.069, 0.0086, 1.77, 131.79, 272.85, 256.2, 0.0060, J2000, false); 

// 2. Specific Dwarfs/Asteroids (JD2026 Epoch - Backward shifted & Exported)
sim.addBodyFromElements("Ceres",  4.726e-10,  2.77, 0.0796, 10.6, 80.2,  73.3,  232,  0.214,   JD2026, false); 
sim.addBodyFromElements("Vesta",  1.303e-10,  2.36, 0.0902, 7.14, 104.0, 152.0, 26.8, 0.272,   JD2026, false); 
sim.addBodyFromElements("Pluto",  6.58e-9,    39.6, 0.2520, 17.1, 110.0, 114.0, 38.7, 0.00396, JD2026, false); 
sim.addBodyFromElements("Juno",   1.44e-11,   2.67, 0.2560, 13.0, 170.0, 248.0, 218.0, 0.226,  JD2026, false); 
sim.addBodyFromElements("Pallas", 1.03e-10,   2.77, 0.2310, 34.9, 173.0, 311.0, 212.0, 0.214,  JD2026, false); 

// 3. Gravitational Perturbers (JD2026 Epoch - Backward shifted & Omitted from Export)
sim.addBodyFromElements("Hygiea",     4.34e-11, 3.15, 0.1080, 3.83, 283.0, 313.0, 217.0, 0.176, JD2026, true); 
sim.addBodyFromElements("Interamnia", 1.76e-11, 3.06, 0.1550, 17.3, 280.0, 94.1,  184.0, 0.184, JD2026, true); 

// Run
sim.runSimulation(10000000, 0.20, 10000, "ZSLTE-26-ephemeris.json");
