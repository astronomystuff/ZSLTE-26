class solarSystemEngine {
    constructor() {
        this.gConstant = 2.959122082855911e-4;
        this.cLight = 173.144543265;
        this.cSquare = this.cLight * this.cLight;
        this.bodies = [];
    }

    // Standard Cartesian Insertion
    addBody(name, mass, x, y, z, vx, vy, vz, isPerturberOnly = false) {
        this.bodies.push({
            name, mass, isPerturberOnly,
            pos: new Float64Array([x, y, z]),
            vel: new Float64Array([vx, vy, vz]),
            acc: new Float64Array([0, 0, 0])
        });
    }

    // Automatically converts Keplerian elements from tables to Cartesian vectors
    addBodyFromElements(name, mass, a, e, iDeg, omegaBigDeg, omegaSmallDeg, mDeg, nDegPerDay, isPerturberOnly = false) {
        // Convert angles from degrees to radians
        const i = iDeg * Math.PI / 180;
        const oBig = omegaBigDeg * Math.PI / 180;
        const oSmall = omegaSmallDeg * Math.PI / 180;
        const M = mDeg * Math.PI / 180;
        const n = nDegPerDay * Math.PI / 180; // Mean motion in rad/day

        // 1. Solve Kepler's Equation: E - e*sin(E) = M using Newton-Raphson iteration
        let E = M;
        for (let iter = 0; iter < 100; iter++) {
            const deltaE = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
            E -= deltaE;
            if (Math.abs(deltaE) < 1e-12) break;
        }

        // 2. Positions and velocities in the orbital plane
        const xPlane = a * (Math.cos(E) - e);
        const yPlane = a * Math.sqrt(1 - e * e) * Math.sin(E);

        const r = a * (1 - e * Math.cos(E));
        const vFactor = (a * a * n) / r;
        const vxPlane = -vFactor * Math.sin(E);
        const vyPlane = vFactor * Math.sqrt(1 - e * e) * Math.cos(E);

        // 3. Transformation matrices to convert to 3D Cartesian coordinates
        const cosO = Math.cos(oBig);
        const sinO = Math.sin(oBig);
        const cosW = Math.cos(oSmall);
        const sinW = Math.sin(oSmall);
        const cosI = Math.cos(i);
        const sinI = Math.sin(i);

        const px = cosW * cosO - sinW * sinO * cosI;
        const py = cosW * sinO + sinW * cosO * cosI;
        const pz = sinW * sinI;

        const qx = -sinW * cosO - cosW * sinO * cosI;
        const qy = -sinW * sinO + cosW * cosO * cosI;
        const qz = cosW * sinI;

        // Final Cartesian Vectors (Relative to the Sun)
        const x = xPlane * px + yPlane * qx;
        const y = xPlane * py + yPlane * qy;
        const z = xPlane * pz + yPlane * qz;

        const vx = vxPlane * px + vyPlane * qx;
        const vy = vxPlane * py + vyPlane * qy;
        const vz = vxPlane * pz + vyPlane * qz;

        this.addBody(name, mass, x, y, z, vx, vy, vz, isPerturberOnly);
    }

    // Shifts all coordinates to ensure the entire system interacts in a net-zero momentum frame
    balanceBarycenter() {
        let totalMass = 0;
        const centerPos = new Float64Array([0, 0, 0]);
        const centerVel = new Float64Array([0, 0, 0]);

        this.bodies.forEach(b => {
            totalMass += b.mass;
            centerPos[0] += b.pos[0] * b.mass;
            centerPos[1] += b.pos[1] * b.mass;
            centerPos[2] += b.pos[2] * b.mass;
            centerVel[0] += b.vel[0] * b.mass;
            centerVel[1] += b.vel[1] * b.mass;
            centerVel[2] += b.vel[2] * b.mass;
        });

        // Shift everything so the system barycenter sits cleanly at (0,0,0)
        this.bodies.forEach(b => {
            b.pos[0] -= centerPos[0] / totalMass;
            b.pos[1] -= centerPos[1] / totalMass;
            b.pos[2] -= centerPos[2] / totalMass;
            b.vel[0] -= centerVel[0] / totalMass;
            b.vel[1] -= centerVel[1] / totalMass;
            b.vel[2] -= centerVel[2] / totalMass;
        });
    }

    computeAccelerations() {
        const totalBodies = this.bodies.length;
        for (let i = 0; i < totalBodies; i++) this.bodies[i].acc.fill(0);

        for (let i = 0; i < totalBodies; i++) {
            const b1 = this.bodies[i];
            for (let j = i + 1; j < totalBodies; j++) {
                const b2 = this.bodies[j];
                const dx = b2.pos[0] - b1.pos[0];
                const dy = b2.pos[1] - b1.pos[1];
                const dz = b2.pos[2] - b1.pos[2];
                const distSq = dx * dx + dy * dy + dz * dz;
                const dist = Math.sqrt(distSq);
                if (dist === 0) continue;
                const forceMagnitude = this.gConstant / (distSq * dist);

                b1.acc[0] += dx * (forceMagnitude * b2.mass);
                b1.acc[1] += dy * (forceMagnitude * b2.mass);
                b1.acc[2] += dz * (forceMagnitude * b2.mass);

                b2.acc[0] -= dx * (forceMagnitude * b1.mass);
                b2.acc[1] -= dy * (forceMagnitude * b1.mass);
                b2.acc[2] -= dz * (forceMagnitude * b1.mass);
            }
        }

        const sun = this.bodies[0];
        const gMassSun = this.gConstant * sun.mass;
        for (let i = 1; i < totalBodies; i++) {
            const planet = this.bodies[i];
            const rx = planet.pos[0] - sun.pos[0];
            const ry = planet.pos[1] - sun.pos[1];
            const rz = planet.pos[2] - sun.pos[2];
            const rSq = rx * rx + ry * ry + rz * rz;
            const rScalar = Math.sqrt(rSq);
            const vx = planet.vel[0] - sun.vel[0];
            const vy = planet.vel[1] - sun.vel[1];
            const vz = planet.vel[2] - sun.vel[2];
            const vSq = vx * vx + vy * vy + vz * vz;
            const rDotV = rx * vx + ry * vy + rz * vz;
            const commonFactor = gMassSun / (this.cSquare * rSq * rScalar);
            const scalarPositionTerm = (4.0 * gMassSun / rScalar) - vSq;
            const scalarVelocityTerm = 4.0 * rDotV;

            const grAx = commonFactor * (scalarPositionTerm * rx + scalarVelocityTerm * vx);
            const grAy = commonFactor * (scalarPositionTerm * ry + scalarVelocityTerm * vy);
            const grAz = commonFactor * (scalarPositionTerm * rz + scalarVelocityTerm * vz);

            planet.acc[0] += grAx; planet.acc[1] += grAy; planet.acc[2] += grAz;
            sun.acc[0] -= grAx * (planet.mass / sun.mass);
            sun.acc[1] -= grAy * (planet.mass / sun.mass);
            sun.acc[2] -= grAz * (planet.mass / sun.mass);
        }
    }

    step(dt) {
        const totalBodies = this.bodies.length;
        for (let i = 0; i < totalBodies; i++) {
            const b = this.bodies[i];
            b.pos[0] += b.vel[0] * dt + 0.5 * b.acc[0] * dt * dt;
            b.pos[1] += b.vel[1] * dt + 0.5 * b.acc[1] * dt * dt;
            b.pos[2] += b.vel[2] * dt + 0.5 * b.acc[2] * dt * dt;
        }
        const oldAcc = this.bodies.map(b => new Float64Array(b.acc));
        this.computeAccelerations();
        for (let i = 0; i < totalBodies; i++) {
            const b = this.bodies[i];
            b.vel[0] += 0.5 * (oldAcc[i][0] + b.acc[0]) * dt;
            b.vel[1] += 0.5 * (oldAcc[i][1] + b.acc[1]) * dt;
            b.vel[2] += 0.5 * (oldAcc[i][2] + b.acc[2]) * dt;
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
    }

    runSimulation(totalYears, dtDays, snapshotIntervalYears, outputFilename = "ZSLTE-26-ephemeris.json") {
        this.balanceBarycenter(); // Converts heliocentric inputs to proper mathematical barycentric space
        this.computeAccelerations(); 

        const stepsPerYear = 365.25 / dtDays;
        const historyData = {};
        this.bodies.forEach(b => { if (!b.isPerturberOnly) historyData[b.name] = []; });

        for (let year = 0; year <= totalYears; year++) {
            for (let step = 0; step < stepsPerYear; step++) this.step(dtDays);

            if (year % snapshotIntervalYears === 0) {
                const sun = this.bodies[0];
                this.bodies.forEach(b => {
                    if (b.isPerturberOnly) return;
                    historyData[b.name].push({
                        t: year,
                        x: Number((b.pos[0] - sun.pos[0]).toFixed(8)), 
                        y: Number((b.pos[1] - sun.pos[1]).toFixed(8)),
                        z: Number((b.pos[2] - sun.pos[2]).toFixed(8))
                    });
                });
            }
        }
        this.downloadEphemeris(historyData, outputFilename);
    }
}
