#!/usr/bin/env python3
"""Compile focused tests against the actual opacity function bodies (no AREPO run).

Run with Python 3 and a C99 compiler: python3 tests/check_bowen_opacity.py
Use --openmp to test the parallel reduction with a compiler supporting -fopenmp.
Use --emit-only DIR to inspect the generated, self-contained C harnesses.
Exit 77 means no supported compiler was available; it is not a test pass.
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


def extract_function(source, name):
    """Find a definition, then balance braces while ignoring comments/strings."""
    pattern = rf"^[ \t]*(?:static\s+)?(?:double|void)\s+{name}\s*\([^;{{]*\)\s*{{"
    match = re.search(pattern, source, re.MULTILINE)
    if not match:
        raise ValueError(f"Function definition not found: {name}")
    opening = source.index("{", match.start())
    tokens = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}]', re.DOTALL)
    depth = 0
    for token in tokens.finditer(source, opening):
        if token.group() == "{":
            depth += 1
        elif token.group() == "}":
            depth -= 1
            if depth == 0:
                return source[match.start():token.end()]
    raise ValueError(f"Unbalanced function body: {name}")


PREAMBLE = r'''
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define RAD_CONST 7.5657e-15
struct sph_particle_data {
  double Density, EOSTemperature, E_rad;
  double KappaRoss, KappaPlanck, KappaScatter, tau;
};
struct {
  int BowenTemperatureMode;
  double BowenOpacityRhomin, BowenOpacityKappamax;
  double DustCondensationTemp, BowenOpacitydT, UnitDensity_in_cgs;
} All;
static struct sph_particle_data cells[3], *SphP = cells;
static int NumGas = 3;
static int active[4] = {-1, 0, 1, 2};
static struct { int NActiveParticles; int *ActiveParticleList; }
  TimeBinsHydro = {4, active};
struct opal_opacity_result { double kappa_ross; };
struct planck_opacity_result { double kappa_planck; };
static int opal_calls, terminate_calls, invalid_count, terminated_in_parallel;
static int calls_at_termination;
static void get_opal_opacity(double rho, double temperature,
                             struct opal_opacity_result *result) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
  opal_calls++;
  result->kappa_ross = rho + temperature * 0.001;
}
static void get_planck_opacity(double rho, double temperature,
                               struct planck_opacity_result *result) {
  result->kappa_planck = rho + temperature * 0.002;
}
static double get_cell_radius(int i) { return 2.0 + i; }
static void terminate(const char *format, ...) {
  va_list arguments;
  (void)format;
  va_start(arguments, format);
  invalid_count = va_arg(arguments, int);
  va_end(arguments);
  terminate_calls++;
  calls_at_termination = opal_calls;
#ifdef _OPENMP
  terminated_in_parallel = omp_in_parallel();
#endif
}
static int checks;
#define CHECK(condition) do { checks++; if(!(condition)) { \
  fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #condition); return 1; \
} } while(0)
static int close_enough(double a, double b) {
  return fabs(a-b) <= 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));
}
static void reset(void) {
  int i;
  All.BowenTemperatureMode = 0;
  All.BowenOpacityRhomin = 1e-15;
  All.BowenOpacityKappamax = 4.5363373545;
  All.DustCondensationTemp = 1100.0;
  All.BowenOpacitydT = 50.0;
  All.UnitDensity_in_cgs = 1.0;
  opal_calls = terminate_calls = invalid_count = terminated_in_parallel = 0;
  calls_at_termination = 0;
  TimeBinsHydro.NActiveParticles = 4;
  for(i=0; i<3; i++) {
    cells[i].Density = 1e-12;
    cells[i].EOSTemperature = 2000.0;
    cells[i].E_rad = RAD_CONST * pow(1100.0, 4.0);
    cells[i].KappaScatter = 999.0;
    cells[i].KappaRoss = 100.0 + i;
    cells[i].KappaPlanck = 200.0 + i;
    cells[i].tau = -99.0;
  }
}
'''

TESTS = r'''
int main(void) {
  int i;
  reset();
#ifdef BOWEN_OPACITY
  {
    const double temperatures[] = {50.0, 1000.0, 1100.0, 2000.0, 1e8};
    for(i=0; i<5; i++) {
      cells[0].EOSTemperature = temperatures[i];
      cells[0].E_rad = NAN;
      CHECK(get_bowen_opacity(&cells[0], 0) ==
        All.BowenOpacityKappamax /
        (1.0 + exp((temperatures[i] - All.DustCondensationTemp) / All.BowenOpacitydT)));
    }
    cells[0].Density = All.BowenOpacityRhomin;
    CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
    cells[0].Density *= 0.5;
    CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
  }
#ifdef IDORT
  reset();
  All.BowenTemperatureMode = 1;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  cells[0].EOSTemperature = 1e7;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  cells[0].E_rad = 0.0;
  cells[0].EOSTemperature = 1100.0;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  CHECK(isnan(get_bowen_opacity(&cells[0], 0)));
  {
    const double invalid[] = {-1.0, NAN, INFINITY};
    for(i=0; i<3; i++) {
      cells[0].E_rad = invalid[i];
      CHECK(isnan(get_bowen_opacity(&cells[0], 0)));
      CHECK(isnan(get_bowen_opacity(&cells[0], 1)));
    }
  }
  cells[0].E_rad = RAD_CONST * pow(1e8, 4.0);
  CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
  cells[0].E_rad = RAD_CONST * pow(1.0, 4.0);
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
    All.BowenOpacityKappamax / (1.0 + exp((1.0-1100.0)/50.0))));
  reset();
  All.BowenTemperatureMode = 1;
  update_opacity(0);
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 3); /* Negative active-list entry was skipped. */
  for(i=0; i<3; i++) {
    CHECK(cells[i].KappaRoss == cells[i].Density + cells[i].EOSTemperature * 0.001);
    CHECK(cells[i].KappaPlanck == cells[i].Density + cells[i].EOSTemperature * 0.002);
    CHECK(close_enough(cells[i].KappaScatter, All.BowenOpacityKappamax / 2));
    CHECK(cells[i].tau == -99.0); /* Startup-only tau refresh is in init.c. */
  }
  reset();
  All.BowenTemperatureMode = 1;
  cells[0].E_rad = -1.0;
  cells[1].E_rad = NAN;
  update_opacity(0);
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(calls_at_termination == 3); /* Errors aggregated after all cells. */
  CHECK(!terminated_in_parallel);
  reset();
  All.BowenTemperatureMode = 1;
  All.BowenOpacityKappamax = 12.0;
  TimeBinsHydro.NActiveParticles = 1; /* Only entry is -1; all cells inactive. */
  refresh_bowen_opacity();
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 0);
  for(i=0; i<3; i++) {
    CHECK(close_enough(cells[i].KappaScatter, 6.0));
    CHECK(cells[i].KappaRoss == 100.0 + i);
    CHECK(cells[i].KappaPlanck == 200.0 + i);
    CHECK(cells[i].tau == (cells[i].KappaRoss + cells[i].KappaScatter) *
          cells[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i));
  }
  cells[0].E_rad = -1.0;
  cells[1].E_rad = NAN;
  cells[2].KappaScatter = 999.0;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(close_enough(cells[2].KappaScatter, 6.0));
  CHECK(!terminated_in_parallel);
#endif
  reset();
  All.BowenOpacityKappamax = 17.0;
  TimeBinsHydro.NActiveParticles = 1;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 0);
  for(i=0; i<3; i++) {
    CHECK(cells[i].KappaScatter == 17.0 /
          (1.0 + exp((cells[i].EOSTemperature-1100.0)/50.0)));
    CHECK(cells[i].KappaRoss == 100.0 + i);
    CHECK(cells[i].KappaPlanck == 200.0 + i);
    CHECK(cells[i].tau == (cells[i].KappaRoss + cells[i].KappaScatter) *
          cells[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i));
  }
#endif
  reset();
  update_opacity(0);
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 3);
  for(i=0; i<3; i++) {
    CHECK(cells[i].tau == -99.0); /* Legacy/non-Bowen update is unchanged. */
    CHECK(cells[i].KappaRoss == cells[i].Density + cells[i].EOSTemperature * 0.001);
    CHECK(cells[i].KappaPlanck == cells[i].Density + cells[i].EOSTemperature * 0.002);
#ifndef BOWEN_OPACITY
    CHECK(cells[i].KappaScatter == 0.0);
#else
    CHECK(cells[i].KappaScatter == All.BowenOpacityKappamax /
          (1.0 + exp((cells[i].EOSTemperature-All.DustCondensationTemp)/All.BowenOpacitydT)));
#endif
  }
  printf("PASS: %d checks\n", checks);
  return 0;
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).resolve().parents[1] / "src" / "mesa_opacity.c")
    parser.add_argument("--emit-only", type=Path, metavar="DIR")
    parser.add_argument("--cc", default=os.environ.get("CC"))
    parser.add_argument("--openmp", action="store_true")
    args = parser.parse_args()
    source = args.source.read_text(encoding="utf-8")
    functions = ("#ifdef BOWEN_OPACITY\n" + extract_function(source, "get_bowen_opacity")
                 + "\n#endif\n" + extract_function(source, "set_opacity") + "\n"
                 + extract_function(source, "update_opacity") + "\n#ifdef BOWEN_OPACITY\n"
                 + extract_function(source, "refresh_bowen_opacity") + "\n#endif\n")
    variants = {"radiation": ("BOWEN_OPACITY", "IDORT"),
                "legacy_no_idort": ("BOWEN_OPACITY",), "no_bowen": ()}
    harnesses = {name: "#define OPACITY_PLANCK\n" + "".join(f"#define {flag}\n" for flag in flags)
                 + PREAMBLE + functions + TESTS for name, flags in variants.items()}
    if args.emit_only:
        args.emit_only.mkdir(parents=True, exist_ok=True)
        for name, body in harnesses.items():
            (args.emit_only / f"{name}.c").write_text(body, encoding="utf-8")
        print(f"Emitted {len(harnesses)} harnesses from {args.source}; no tests run.")
        return 0
    compiler = args.cc or next((shutil.which(cc) for cc in ("cc", "gcc", "clang") if shutil.which(cc)), None)
    if not compiler:
        print("UNAVAILABLE: no C99 compiler found; tests were not run (exit 77).", file=sys.stderr)
        return 77
    command = [compiler] if Path(compiler).is_file() else shlex.split(compiler)
    if not command or not shutil.which(command[0]):
        print(f"UNAVAILABLE: compiler cannot be found: {compiler} (exit 77).", file=sys.stderr)
        return 77
    with tempfile.TemporaryDirectory(prefix="bowen-opacity-check-") as directory:
        root = Path(directory)
        for name, body in harnesses.items():
            c_file, executable = root / f"{name}.c", root / (name + (".exe" if os.name == "nt" else ""))
            c_file.write_text(body, encoding="utf-8")
            compile_command = command + ["-std=c99", "-O2"]
            if args.openmp:
                compile_command.append("-fopenmp")
            compile_command += [str(c_file), "-lm", "-o", str(executable)]
            subprocess.run(compile_command, check=True)
            print(f"{name}:", flush=True)
            subprocess.run([str(executable)], check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
