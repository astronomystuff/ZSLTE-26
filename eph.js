class solarSystemEngine {
    constructor() {
        this.gConstant = 2.959122082855911e-4;
        this.cLight = 173.144543265;
        this.cSquare = this.cLight * this.cLight;
        this.bodies = [];
    }

    addBody(name, mass, x, y, z, vx, vy, vz, isPerturberOnly = false) {
        this.bodies.push({
            name,
            mass,
            isPerturberOnly,
            pos: new Float64Array([x, y, z]),
            vel: new Float64Array([vx, vy, vz]),
            acc: new Float64Array([0, 0, 0])
        });
    }

    computeAccelerations() {
        const totalBodies = this.bodies.length;
        for (let i = 0; i < totalBodies; i++) {
            this.bodies[i].acc.fill(0);
        }

        // Newtonian Mutual Perturbations
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

        // 1PN General Relativity
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

            planet.acc[0] += grAx;
            planet.acc[1] += grAy;
            planet.acc[2] += grAz;

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

    runTest(totalYears, dtDays, snapshotIntervalYears) {
        this.computeAccelerations();
        const stepsPerYear = 365.25 / dtDays;
        const historyData = {};

        this.bodies.forEach(b => {
            if (!b.isPerturberOnly) historyData[b.name] = [];
        });

        console.log(`Running ZSLTE-26 test simulation...`);

        for (let year = 0; year <= totalYears; year++) {
            for (let step = 0; step < stepsPerYear; step++) {
                this.step(dtDays);
            }

            if (year % snapshotIntervalYears === 0) {
                const sun = this.bodies[0];
                this.bodies.forEach(b => {
                    if (b.isPerturberOnly) return;
                    const helioX = b.pos[0] - sun.pos[0];
                    const helioY = b.pos[1] - sun.pos[1];
                    const helioZ = b.pos[2] - sun.pos[2];

                    historyData[b.name].push({
                        t: year,
                        x: Number(helioX.toFixed(8)),
                        y: Number(helioY.toFixed(8)),
                        z: Number(helioZ.toFixed(8))
                    });
                });
            }
        }
        
        console.log("Simulation complete. Output Data Structure:");
        console.log(historyData);
    }
}

// Instantiate the test engine
const sim = new solarSystemEngine();

// Seed baseline values (Sun, Earth, Mars, and Hygeia as a pure perturber)
sim.addBody("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
sim.addBody("Earth", 3.003489e-6, -0.1685, 0.9684, -0.00004, -0.0172, -0.0030, 0.0);
sim.addBody("Mars", 3.227151e-7, 1.324, -0.543, -0.021, 0.002, 0.014, 0.0001);
sim.addBody("Hygeia", 4.34e-11, 2.5, 1.2, -0.3, -0.008, 0.007, 0.001, true);

// Run for 100 years, 0.25 day steps, logging every 25 years
sim.runTest(100, 0.25, 25);
