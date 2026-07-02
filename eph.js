export function generateEphemeris(initStateAtJ2000) {
  const DAY   = 86400;
  const YEAR  = 365.25 * DAY;
  const G     = 6.67430e-11;
  const c2    = 299792458 ** 2;

  const T0    = -5_000_000 * YEAR;
  const T1    =  5_000_000 * YEAR;
  const STEP  = 2 * DAY;

  const bodies = [
    { name:"Sun",     m:1.98847e30 },
    { name:"Mercury", m:3.3011e23 },
    { name:"Venus",   m:4.8675e24 },
    { name:"Earth",   m:5.97237e24 },
    { name:"Mars",    m:6.4171e23 },
    { name:"Jupiter", m:1.8982e27 },
    { name:"Saturn",  m:5.6834e26 },
    { name:"Uranus",  m:8.6810e25 },
    { name:"Neptune", m:1.02413e26 },
    { name:"Pluto",   m:1.303e22 },
    { name:"Ceres",   m:9.393e20 },
    { name:"Vesta",   m:2.590e20 },
    { name:"Pallas",  m:2.110e20 },
    { name:"Juno",    m:2.670e19 },
    { name:"Hygieia", m:8.67e19 },
    { name:"Interamnia", m:3.5e19 },
    { name:"Eunomia", m:3.12e19 }
  ];
  const N = bodies.length;

  const cadence = {
    Sun:4*DAY, Mercury:8*DAY, Venus:12*DAY, Earth:12*DAY,
    Mars:16*DAY, Jupiter:16*DAY, Saturn:28*DAY,
    Uranus:32*DAY, Neptune:36*DAY, Pluto:36*DAY,
    Ceres:16*DAY, Vesta:16*DAY, Pallas:16*DAY,
    Juno:16*DAY, Hygieia:16*DAY, Interamnia:16*DAY, Eunomia:16*DAY
  };

  const degree = {
    Mercury:48, Venus:48, Earth:48, Mars:48,
    Jupiter:32, Saturn:32, Uranus:32, Neptune:32, Pluto:32,
    Sun:32, Ceres:32, Vesta:32, Pallas:32,
    Juno:32, Hygieia:32, Interamnia:32, Eunomia:32
  };

  const blockLen = {
    Mercury:10*YEAR, Venus:10*YEAR, Earth:10*YEAR, Mars:10*YEAR,
    Jupiter:50*YEAR, Saturn:50*YEAR, Uranus:100*YEAR,
    Neptune:100*YEAR, Pluto:100*YEAR,
    Sun:50*YEAR,
    Ceres:20*YEAR, Vesta:20*YEAR, Pallas:20*YEAR,
    Juno:20*YEAR, Hygieia:20*YEAR, Interamnia:20*YEAR, Eunomia:20*YEAR
  };

  let state = initStateAtJ2000.map(s => ({ ...s })); // {x,y,z,vx,vy,vz}
  const iSun = bodies.findIndex(b => b.name === "Sun");

  // ---------- acceleration (Newtonian + Sun 1PN) ----------
  function computeAcc(st) {
    const acc = Array.from({ length:N }, () => ({ x:0,y:0,z:0 }));
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) if (i !== j) {
        const dx = st[j].x - st[i].x;
        const dy = st[j].y - st[i].y;
        const dz = st[j].z - st[i].z;
        const r2 = dx*dx + dy*dy + dz*dz;
        const r  = Math.sqrt(r2);
        const mu = G * bodies[j].m;
        const fN = mu / (r2 * r);
        let ax = fN * dx;
        let ay = fN * dy;
        let az = fN * dz;

        if (bodies[j].name === "Sun") {
          const v2 = st[i].vx**2 + st[i].vy**2 + st[i].vz**2;
          const rv = dx*st[i].vx + dy*st[i].vy + dz*st[i].vz;
          const fac = mu / (c2 * r2 * r);
          ax += fac * ((4*mu/r - v2)*dx + 4*rv*st[i].vx);
          ay += fac * ((4*mu/r - v2)*dy + 4*rv*st[i].vy);
          az += fac * ((4*mu/r - v2)*dz + 4*rv*st[i].vz);
        }

        acc[i].x += ax;
        acc[i].y += ay;
        acc[i].z += az;
      }
    }
    return acc;
  }

  // ---------- RK4 startup ----------
  function rk4Step(st, dt) {
    const k1v = computeAcc(st);
    const k1x = st.map((s,i) => ({ x:s.vx, y:s.vy, z:s.vz }));

    const st2 = st.map((s,i) => ({
      x: s.x + 0.5*dt*k1x[i].x,
      y: s.y + 0.5*dt*k1x[i].y,
      z: s.z + 0.5*dt*k1x[i].z,
      vx: s.vx + 0.5*dt*k1v[i].x,
      vy: s.vy + 0.5*dt*k1v[i].y,
      vz: s.vz + 0.5*dt*k1v[i].z
    }));
    const k2v = computeAcc(st2);
    const k2x = st2.map((s,i) => ({ x:s.vx, y:s.vy, z:s.vz }));

    const st3 = st.map((s,i) => ({
      x: s.x + 0.5*dt*k2x[i].x,
      y: s.y + 0.5*dt*k2x[i].y,
      z: s.z + 0.5*dt*k2x[i].z,
      vx: s.vx + 0.5*dt*k2v[i].x,
      vy: s.vy + 0.5*dt*k2v[i].y,
      vz: s.vz + 0.5*dt*k2v[i].z
    }));
    const k3v = computeAcc(st3);
    const k3x = st3.map((s,i) => ({ x:s.vx, y:s.vy, z:s.vz }));

    const st4 = st.map((s,i) => ({
      x: s.x + dt*k3x[i].x,
      y: s.y + dt*k3x[i].y,
      z: s.z + dt*k3x[i].z,
      vx: s.vx + dt*k3v[i].x,
      vy: s.vy + dt*k3v[i].y,
      vz: s.vz + dt*k3v[i].z
    }));
    const k4v = computeAcc(st4);
    const k4x = st4.map((s,i) => ({ x:s.vx, y:s.vy, z:s.vz }));

    const out = st.map((s,i) => ({
      x: s.x + dt*(k1x[i].x + 2*k2x[i].x + 2*k3x[i].x + k4x[i].x)/6,
      y: s.y + dt*(k1x[i].y + 2*k2x[i].y + 2*k3x[i].y + k4x[i].y)/6,
      z: s.z + dt*(k1x[i].z + 2*k2x[i].z + 2*k3x[i].z + k4x[i].z)/6,
      vx: s.vx + dt*(k1v[i].x + 2*k2v[i].x + 2*k3v[i].x + k4v[i].x)/6,
      vy: s.vy + dt*(k1v[i].y + 2*k2v[i].y + 2*k3v[i].y + k4v[i].y)/6,
      vz: s.vz + dt*(k1v[i].z + 2*k2v[i].z + 2*k3v[i].z + k4v[i].z)/6
    }));
    return out;
  }

  // ---------- Gauss–Jackson history ----------
  const K = 7;

  const hist = {
    t:   new Array(K),
    pos: new Array(K),
    vel: new Array(K),
    acc: new Array(K)
  };

  let st = initStateAtJ2000.map(s => ({ ...s }));
  for (let i = 0; i < K; i++) st = rk4Step(st, -STEP);

  let tHist = T0 - K*STEP;
  for (let k = 0; k < K; k++) {
    hist.t[k]   = tHist;
    hist.pos[k] = st.map(s => ({ x:s.x, y:s.y, z:s.z }));
    hist.vel[k] = st.map(s => ({ vx:s.vx, vy:s.vy, vz:s.vz }));
    const a     = computeAcc(st);
    hist.acc[k] = a.map(v => ({ x:v.x, y:v.y, z:v.z }));

    st    = rk4Step(st, STEP);
    tHist += STEP;
  }
  state = st;

  // ---------- GJ8 coefficients ----------
  const aP = [
    -9/840,
     128/840,
    -1008/840,
     4032/840,
    -1008/840,
     128/840,
    -9/840
  ];

  const aC = [
     9/840,
    -128/840,
     1008/840,
    -4032/840,
     1008/840,
    -128/840,
     9/840
  ];

  function gjStep() {
    const posK = hist.pos[K-1];
    const velK = hist.vel[K-1];
    const acc  = hist.acc;

    // predictor
    const pred = posK.map((p,i) => {
      let ax = 0, ay = 0, az = 0;
      for (let j = 0; j < K; j++) {
        ax += aP[j] * acc[K-1-j][i].x;
        ay += aP[j] * acc[K-1-j][i].y;
        az += aP[j] * acc[K-1-j][i].z;
      }
      return {
        x: p.x + velK[i].vx*STEP + ax*STEP*STEP,
        y: p.y + velK[i].vy*STEP + ay*STEP*STEP,
        z: p.z + velK[i].vz*STEP + az*STEP*STEP,
        vx: velK[i].vx + ax*STEP,
        vy: velK[i].vy + ay*STEP,
        vz: velK[i].vz + az*STEP
      };
    });

    const accPred = computeAcc(pred);

    // corrector
    const corr = posK.map((p,i) => {
      let ax = 0, ay = 0, az = 0;
      for (let j = 0; j < K; j++) {
        ax += aC[j] * acc[K-1-j][i].x;
        ay += aC[j] * acc[K-1-j][i].y;
        az += aC[j] * acc[K-1-j][i].z;
      }
      ax += aC[0] * accPred[i].x;
      ay += aC[0] * accPred[i].y;
      az += aC[0] * accPred[i].z;

      return {
        x: p.x + velK[i].vx*STEP + ax*STEP*STEP,
        y: p.y + velK[i].vy*STEP + ay*STEP*STEP,
        z: p.z + velK[i].vz*STEP + az*STEP*STEP,
        vx: velK[i].vx + ax*STEP,
        vy: velK[i].vy + ay*STEP,
        vz: velK[i].vz + az*STEP
      };
    });

    // shift history
    for (let k = 0; k < K-1; k++) {
      hist.t[k]   = hist.t[k+1];
      hist.pos[k] = hist.pos[k+1];
      hist.vel[k] = hist.vel[k+1];
      hist.acc[k] = hist.acc[k+1];
    }
    const tNext = hist.t[K-2] + STEP;
    hist.t[K-1]   = tNext;
    hist.pos[K-1] = corr.map(s => ({ x:s.x, y:s.y, z:s.z }));
    hist.vel[K-1] = corr.map(s => ({ vx:s.vx, vy:s.vy, vz:s.vz }));
    hist.acc[K-1] = accPred.map(a => ({ x:a.x, y:a.y, z:a.z }));

    state = corr;
  }

  // ---------- interpolation + Chebyshev ----------
  function interpolateHermite(samples, nodes, key) {
    const vals = [];
    for (const tn of nodes) {
      let i = 0;
      while (i+1 < samples.length && samples[i+1].t < tn) i++;
      const p0 = samples[Math.max(0,i)];
      const p1 = samples[Math.min(samples.length-1,i+1)];
      const t0 = p0.t, t1 = p1.t;
      if (t1 === t0) { vals.push(p0[key]); continue; }
      const u = (tn - t0)/(t1 - t0);
      const v0 = p0[key], v1 = p1[key];
      const h = v0*(2*u*u*u - 3*u*u + 1) + v1*(-2*u*u*u + 3*u*u);
      vals.push(h);
    }
    return vals;
  }

  function dctCheby(vals) {
    const n = vals.length;
    const c = new Array(n).fill(0);
    for (let k = 0; k < n; k++) {
      let sum = 0;
      for (let j = 0; j < n; j++) {
        const theta = Math.PI * (j + 0.5) / n;
        sum += vals[j] * Math.cos(k * theta);
      }
      c[k] = (2 / n) * sum;
    }
    c[0] *= 0.5;
    return c;
  }

  // ---------- per-body streaming ----------
  const perBody = bodies.map(b => ({
    name: b.name,
    nextSave: T0,
    blockStart: T0,
    blockEnd: T0 + blockLen[b.name],
    samples: [],
    blocks: []
  }));

  for (let t = T0; t <= T1; t += STEP) {
    gjStep();

    // sampling
    for (const B of perBody) {
      if (t >= B.nextSave) {
        const iBody = bodies.findIndex(bb => bb.name === B.name);
        const hx = state[iBody].x - state[iSun].x;
        const hy = state[iBody].y - state[iSun].y;
        const hz = state[iBody].z - state[iSun].z;
        B.samples.push({ t, x:hx, y:hy, z:hz });
        B.nextSave += cadence[B.name];
      }
    }

    // block completion
    for (const B of perBody) {
      if (t >= B.blockEnd && B.samples.length) {
        const deg = degree[B.name];
        const n   = deg + 1;
        const t0  = B.blockStart;
        const t1  = B.blockEnd;
        const tMid = 0.5 * (t0 + t1);
        const dt   = 0.5 * (t1 - t0);

        const nodes = [];
        for (let k = 0; k < n; k++) {
          const theta = Math.PI * (k + 0.5) / n;
          const s = Math.cos(theta);
          nodes.push(tMid + dt * s);
        }

        const xVals = interpolateHermite(B.samples, nodes, "x");
        const yVals = interpolateHermite(B.samples, nodes, "y");
        const zVals = interpolateHermite(B.samples, nodes, "z");

        const cx = dctCheby(xVals);
        const cy = dctCheby(yVals);
        const cz = dctCheby(zVals);

        B.blocks.push({ t0, t1, degree:deg, cx, cy, cz });

        B.blockStart = B.blockEnd;
        B.blockEnd   = B.blockStart + blockLen[B.name];
        B.samples    = [];
      }
    }
  }

  return perBody.map(B => ({
    body: B.name,
    frame: "heliocentric",
    blocks: B.blocks
  }));
}
