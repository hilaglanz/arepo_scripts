// Offline sensitivity audit: reads saved text profiles, never snapshots or AREPO source.
// Run: node mean_profile_optical_depth.mjs [profile_directory] [output_directory]
import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';

function interpolate(rows, x, logY = false, logX = false) {
  if (x < rows[0][0] || x > rows.at(-1)[0]) throw new Error('Interpolation outside input domain');
  let lo = 0, hi = rows.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (rows[mid][0] <= x) lo = mid; else hi = mid;
  }
  const f = logX ? Math.log(x / rows[lo][0]) / Math.log(rows[hi][0] / rows[lo][0])
    : (x - rows[lo][0]) / (rows[hi][0] - rows[lo][0]);
  return logY ? Math.exp(Math.log(rows[lo][1]) + f * Math.log(rows[hi][1] / rows[lo][1]))
    : rows[lo][1] + f * (rows[hi][1] - rows[lo][1]);
}

export function calculateAudit(density, gasTemperature, grainTemperature, mie) {
  const Rstar = 1.748342799e13, gate = 1e-15, Tc = 1100, width = 50, dustRatio = 0.001;
  const outer = density.at(-1)[0];
  const radii = [1, 1.25, 1.4, 1.5, 1.5855, 1.75, 2, 2.08182, 3, 4, 5, 10];
  const starts = radii.map(x => x * Rstar);
  const coefficients = {
    scatteringTransport: mie.coefficient * dustRatio,
    rawScattering: mie.rawScatteringCoefficient * dustRatio,
    absorptionStellarSpectrum: mie.absorptionCoefficient * dustRatio,
    legacyBowenScattering: 14.821
  };
  const activation = temperature => 1 / (1 + Math.exp((temperature - Tc) / width));

  function integrate(logRho, subdivisions) {
    let grid = [Rstar, outer, ...starts, ...density.map(v => v[0]),
      ...gasTemperature.map(v => v[0]), ...grainTemperature.map(v => v[0] * Rstar)]
      .filter(x => x >= Rstar && x <= outer);
    // Split exactly at each density-gate crossing for the chosen interpolant.
    for (let i = 0; i < density.length - 1; i++) {
      const [x0, y0] = density[i], [x1, y1] = density[i + 1];
      if ((y0 - gate) * (y1 - gate) < 0) {
        const f = logRho ? Math.log(gate / y0) / Math.log(y1 / y0) : (gate - y0) / (y1 - y0);
        const x = x0 + f * (x1 - x0);
        if (x >= Rstar && x <= outer) grid.push(x);
      }
    }
    grid = [...new Set(grid)].sort((a, b) => a - b);
    const sums = { grain: Array(grid.length).fill(0), gas: Array(grid.length).fill(0),
      allGas: Array(grid.length).fill(0) };
    function value(x, kind) {
      const rho = interpolate(density, x, logRho);
      if (kind === 'allGas') return rho;
      const temperature = kind === 'grain' ? interpolate(grainTemperature, x / Rstar, true, true)
        : interpolate(gasTemperature, x);
      return rho * activation(temperature);
    }
    for (let i = grid.length - 2; i >= 0; i--) {
      const a = grid[i], b = grid[i + 1];
      const active = interpolate(density, (a + b) / 2, logRho) > gate;
      for (const kind of Object.keys(sums)) {
        let integral = 0;
        if (active || kind === 'allGas') {
          const dx = (b - a) / subdivisions;
          for (let j = 0; j < subdivisions; j++) {
            const l = a + j * dx, h = l + dx, m = (l + h) / 2;
            integral += dx * (value(l, kind) + 4 * value(m, kind) + value(h, kind)) / 6;
          }
        }
        sums[kind][i] = sums[kind][i + 1] + integral;
      }
    }
    const total = sums.grain[0];
    return starts.map((r, n) => {
      const i = grid.indexOf(r), column = sums.grain[i], legacyColumn = sums.gas[i];
      return {
        rOverRstar: radii[n], radiusCm: r, meanDensity: interpolate(density, r, logRho),
        grainTemperatureK: interpolate(grainTemperature, r / Rstar, true, true),
        activationFraction: activation(interpolate(grainTemperature, r / Rstar, true, true)),
        gasColumnGcm2: sums.allGas[i], grainWeightedColumnGcm2: column,
        legacyGasTemperatureWeightedColumnGcm2: legacyColumn,
        tauScatteringTransport: column * coefficients.scatteringTransport,
        tauAbsorptionStellarSpectrum: column * coefficients.absorptionStellarSpectrum,
        tauRawScattering: column * coefficients.rawScattering,
        tauLegacyBowenScattering: legacyColumn * coefficients.legacyBowenScattering,
        incomingTauScatteringTransport: (total - column) * coefficients.scatteringTransport,
        incomingTauAbsorptionStellarSpectrum: (total - column) * coefficients.absorptionStellarSpectrum
      };
    });
  }
  const runs = [1, 4, 16].map(n => ({ subdivisionsPerUnionInterval: n, rows: integrate(false, n) }));
  const linear = runs.at(-1).rows, logRho = integrate(true, 16);
  return {
    scope: 'Toy radial optical depths on a frozen time-averaged profile; not measured simulation optical depths or a self-consistent new atmosphere.',
    formula: 'Sigma_f(r) = integral_r^r_outer rho_bar(r) * logistic[(Td(r)-1100)/50] * 1[rho_bar(r)>1e-15] dr',
    constants: { radiusCm: Rstar, outerRadiusCm: outer, outerRadiusRstar: outer / Rstar,
      condensationTemperatureK: Tc, logisticWidthK: width, strictDensityGateGcm3: gate,
      dustGasMassRatio: dustRatio, coefficientsCm2PerGgas: coefficients },
    interpolation: 'Main result: linear rho(r), linear gas T(r), log(Td) versus log(r). Sensitivity: log(rho) versus linear r. Composite Simpson integration on the union of input knots, requested radii and exact gate crossings.',
    weights: 'The logistic is used at all r>=Rstar, including partial condensation before Td=1100 K. It is 0.5 near 1.5855 Rstar; no hard temperature or 1.5855-Rstar radial cutoff is added. Cells represented by rho<=1e-15 receive zero dust weight.',
    linearDensityResults: linear,
    logDensitySensitivity: logRho,
    quadratureConvergence: runs.map(run => ({ subdivisionsPerUnionInterval: run.subdivisionsPerUnionInterval,
      totalWeightedColumnGcm2: run.rows[0].grainWeightedColumnGcm2,
      totalTauScatteringTransport: run.rows[0].tauScatteringTransport,
      totalLegacyTauScattering: run.rows[0].tauLegacyBowenScattering,
      outwardFrom2RWeightedColumnGcm2: run.rows.find(v => v.rOverRstar === 2).grainWeightedColumnGcm2 })),
    sensitivityPercent: 100 * (logRho[0].grainWeightedColumnGcm2 / linear[0].grainWeightedColumnGcm2 - 1),
    cautions: [
      'Opacity evaluated on averaged gas temperature/density is not the average cellwise opacity; shocks, clumps and correlations are absent.',
      'The atmosphere is frozen while the new grain-temperature assignment is imposed; dynamics and grain growth could change the density and dust distribution.',
      'Stellar-spectrum flux-mean coefficients are held fixed; spectral transfer, diffuse radiation and backwarming are not solved.',
      'Large scattering optical depth changes angular intensity; conserved spherical bolometric flux is not exponentially attenuated by scattering alone.',
      'The thin-temperature lookup fails this frozen-profile thinness check. Use it as a diagnostic closure, not an already validated production dust-temperature model.'
    ]
  };
}

export function markdown(result) {
  const rows = result.linearDensityResults.filter(v => [1, 1.5855, 2, 3].includes(v.rOverRstar));
  const total = result.linearDensityResults[0], logTotal = result.logDensitySensitivity[0];
  return `# Frozen mean-profile optical-depth sensitivity audit\n\n` +
    `This is an offline toy calculation on the saved snapshots-500–571 averaged atmosphere, not an actual snapshot optical depth or a relaxed new donor. No snapshots or AREPO source are modified.\n\n` +
    `For dust/gas=0.001, integrate rho_bar times the logistic grain activation and the strict rho_bar>1e-15 gate. The logistic width is 50 K around 1100 K. Partial activation **inside** the nominal 1.5855 Rstar condensation radius is retained; there is no hard cutoff there. The upper integration radius is ${result.constants.outerRadiusRstar.toFixed(3)} Rstar.\n\n` +
    `| Outward from r/Rstar | Weighted gas column (g/cm²) | Scattering transport tau | Stellar-spectrum absorption tau | Raw scattering tau |\n|---:|---:|---:|---:|---:|\n` +
    rows.map(v => `| ${v.rOverRstar} | ${v.grainWeightedColumnGcm2.toPrecision(5)} | ${v.tauScatteringTransport.toPrecision(5)} | ${v.tauAbsorptionStellarSpectrum.toPrecision(5)} | ${v.tauRawScattering.toPrecision(5)} |`).join('\n') +
    `\n\nThe old gas-temperature logistic on the same mean profile has weighted column ${total.legacyGasTemperatureWeightedColumnGcm2.toPrecision(6)} g/cm² and tau=${total.tauLegacyBowenScattering.toPrecision(6)} with kappa_max=14.821 cm²/g. This also is opacity evaluated on averages, not mean cellwise opacity.\n\n` +
    `Linear-density interpolation gives total weighted column ${total.grainWeightedColumnGcm2.toPrecision(7)} g/cm²; log-density interpolation gives ${logTotal.grainWeightedColumnGcm2.toPrecision(7)} g/cm² (${result.sensitivityPercent.toFixed(2)}%).\n\n` +
    `| Simpson subintervals per union-knot interval | Total weighted column (g/cm²) | Total transport tau | Legacy scattering tau |\n|---:|---:|---:|---:|\n` +
    result.quadratureConvergence.map(v => `| ${v.subdivisionsPerUnionInterval} | ${v.totalWeightedColumnGcm2.toFixed(8)} | ${v.totalTauScatteringTransport.toFixed(8)} | ${v.totalLegacyTauScattering.toFixed(9)} |`).join('\n') +
    `\n\nThe input profile spacing is about 0.0711 Rstar. Integration splits at input profile and lookup knots and density-gate crossings before Simpson refinement; convergence here establishes quadrature accuracy only, not physical resolution.\n\n` +
    result.cautions.map(v => `- ${v}`).join('\n') +
    `\n\n## Reproduction and provenance\n\nRun \`node mean_profile_optical_depth.mjs [profile_directory] [output_directory]\`. Defaults use the original Dropbox profile directory and this script's directory. The adjacent JSON records constants, input paths and SHA-256 hashes, detailed inward/outward columns, and convergence.\n\n` +
    result.provenance.inputs.map(v => `- ${v.role}: \`${v.path}\`; SHA-256 \`${v.sha256}\`.`).join('\n') + '\n';
}

// CLI: all mutations are restricted to the two named audit outputs.
if (typeof process !== 'undefined') {
const scriptDirectory = path.dirname(fileURLToPath(import.meta.url));
const profileDirectory = process.argv[2] || String.raw`C:\Users\hila\Dropbox\RMHD\IDORT\AGB\1.5\1e-8\dust`;
const outputDirectory = process.argv[3] || scriptDirectory;
const inputFiles = [
  ['density', path.join(profileDirectory, 'rho_profile__500_to_571_time_avg.txt')],
  ['gasTemperature', path.join(profileDirectory, 'temp_profile__500_to_571_time_avg.txt')],
  ['grainTemperature', path.join(scriptDirectory, 'mie_Td_Mg2SiO4_a0p3um_T3061.dat')],
  ['mieCoefficients', path.join(scriptDirectory, 'results_subdivisions_4.json')]
];
const sources = await Promise.all(inputFiles.map(async ([role, filename]) => ({ role, path: filename,
  text: await fs.readFile(filename, 'utf8') })));
const csv = text => text.trim().split(/\r?\n/).slice(1).map(line => line.split(',').map(Number));
const table = text => text.split(/\r?\n/).filter(line => line.trim() && !line.startsWith('#'))
  .map(line => line.trim().split(/\s+/).map(Number));
const findGrain = obj => {
  if (obj && typeof obj === 'object' && obj.radiusUm === 0.3) return obj;
  if (obj && typeof obj === 'object') for (const value of Object.values(obj)) {
    const result = findGrain(value); if (result) return result;
  }
};
const mie = findGrain(JSON.parse(sources[3].text));
if (!mie) throw new Error('Missing 0.3-micron Mie result');
const result = calculateAudit(csv(sources[0].text), csv(sources[1].text), table(sources[2].text), mie);
result.provenance = { script: fileURLToPath(import.meta.url), inputs: sources.map(({ role, path, text }) =>
  ({ role, path, sha256: createHash('sha256').update(text).digest('hex') })) };
await fs.mkdir(outputDirectory, { recursive: true });
await fs.writeFile(path.join(outputDirectory, 'mean_profile_optical_depth.json'), JSON.stringify(result, null, 2) + '\n');
await fs.writeFile(path.join(outputDirectory, 'mean_profile_optical_depth.md'), markdown(result));
console.log(JSON.stringify({ totalTransportTau: result.linearDensityResults[0].tauScatteringTransport,
  legacyScatteringTau: result.linearDensityResults[0].tauLegacyBowenScattering,
  logDensityColumnDifferencePercent: result.sensitivityPercent }));
}
