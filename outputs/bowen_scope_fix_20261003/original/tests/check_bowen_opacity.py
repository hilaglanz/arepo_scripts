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
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define RAD_CONST 7.5657e-15
#define IDORT_NUM_ANGLES 4
struct sph_particle_data {
  double Density, EOSTemperature, E_rad;
  double KappaRoss, KappaPlanck, KappaScatter, tau;
  double Volume, Intensity[IDORT_NUM_ANGLES], IV[IDORT_NUM_ANGLES];
};
static struct particle_data {
  unsigned long long ID;
  int Type;
  double Mass, Pos[3];
} P[3];
static int ThisTask = 29, RestartFlag = 1;
struct {
  int BowenTemperatureMode;
  double BowenOpacityRhomin, BowenOpacityKappamax;
  double DustCondensationTemp, BowenOpacitydT, UnitDensity_in_cgs, Time;
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
static char diagnostic_text[8192], terminate_text[1024];
static int diagnostic_calls, diagnostic_in_parallel;
static int capture_fprintf(FILE *stream, const char *format, ...) {
  va_list arguments;
  const size_t used = strlen(diagnostic_text);
  int written;
  (void)stream;
#ifdef _OPENMP
  diagnostic_in_parallel |= omp_in_parallel();
#endif
  diagnostic_calls++;
  va_start(arguments, format);
  written = vsnprintf(diagnostic_text + used, sizeof(diagnostic_text) - used,
                      format, arguments);
  va_end(arguments);
  return written;
}
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
  va_start(arguments, format);
  vsnprintf(terminate_text, sizeof(terminate_text), format, arguments);
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
  All.Time = 200050555.8;
  opal_calls = terminate_calls = invalid_count = terminated_in_parallel = 0;
  calls_at_termination = 0;
  diagnostic_calls = diagnostic_in_parallel = 0;
  diagnostic_text[0] = terminate_text[0] = '\0';
  TimeBinsHydro.NActiveParticles = 4;
  for(i=0; i<3; i++) {
    cells[i].Density = 1e-12;
    cells[i].EOSTemperature = 2000.0;
    cells[i].E_rad = RAD_CONST * pow(1100.0, 4.0);
    cells[i].KappaScatter = 999.0;
    cells[i].KappaRoss = 100.0 + i;
    cells[i].KappaPlanck = 200.0 + i;
    cells[i].tau = -99.0;
    cells[i].Volume = 16.0 + i;
    P[i].ID = 10000000000ULL + i;
    P[i].Type = 0;
    P[i].Mass = 123.0 + i;
    P[i].Pos[0] = 1.0 + i;
    P[i].Pos[1] = 2.0 + i;
    P[i].Pos[2] = 3.0 + i;
    for(int angle=0; angle<IDORT_NUM_ANGLES; angle++) {
      cells[i].Intensity[angle] = 1.0;
      cells[i].IV[angle] = cells[i].Volume;
    }
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
  cells[0].EOSTemperature = NAN;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  cells[0].E_rad = 0.0;
  cells[0].EOSTemperature = 1100.0;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
    All.BowenOpacityKappamax / (1.0 + exp(-1100.0/50.0))));
  {
    const double invalid_gas[] = {0.0, -1.0, NAN, INFINITY};
    for(i=0; i<4; i++) {
      cells[0].EOSTemperature = invalid_gas[i];
      CHECK(isnan(get_bowen_opacity(&cells[0], 1)));
      CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
        All.BowenOpacityKappamax / (1.0 + exp(-1100.0/50.0))));
    }
  }
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
  cells[0].E_rad = DBL_MAX; /* E_rad/RAD_CONST would overflow. */
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
  cells[0].Intensity[0] = 0.0;
  cells[0].Intensity[1] = -2.0;
  cells[0].Intensity[2] = NAN;
  cells[0].Intensity[3] = 3.0;
  cells[0].IV[0] = INFINITY;
  update_opacity(0);
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(calls_at_termination == 3); /* Errors aggregated after all cells. */
  CHECK(!terminated_in_parallel);
  CHECK(!diagnostic_in_parallel);
  CHECK(diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "context=update_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "RestartFlag=1 bootstrap=0 i=0 ID=10000000000 type=0 mass=123") != NULL);
  CHECK(strstr(diagnostic_text, "pos=(1,2,3)") != NULL);
  CHECK(strstr(diagnostic_text, "rho=") != NULL && strstr(diagnostic_text, "Tgas=2000") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=-1 Volume=16") != NULL);
  CHECK(strstr(diagnostic_text, "stored_I_min_finite=-2 stored_I_max_finite=3 angles=4") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=1 negative_angles=1 nonfinite_angles=1 nonfinite_IV=1") != NULL);
  CHECK(strstr(terminate_text, "finite and nonnegative") != NULL);
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
  for(i=0; i<IDORT_NUM_ANGLES; i++) cells[0].Intensity[i] = NAN;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(close_enough(cells[2].KappaScatter, 6.0));
  CHECK(!terminated_in_parallel);
  CHECK(!diagnostic_in_parallel);
  CHECK(diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "context=refresh_bowen_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "i=0 ID=10000000000") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=-1 Volume=16") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=0 negative_angles=0 nonfinite_angles=4 nonfinite_IV=0") != NULL);

  /* Excluded/bootstrap zeros must not consume the one-time warning. */
  reset();
  All.BowenTemperatureMode = 1;
  cells[0].E_rad = 0.0;
  cells[0].Density = All.BowenOpacityRhomin;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);
  CHECK(cells[0].KappaScatter == 0.0);
  cells[0].Density = 1e-12;
  cells[0].EOSTemperature = 1100.0;
  update_opacity(1);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);
  CHECK(close_enough(cells[0].KappaScatter, All.BowenOpacityKappamax / 2));
  All.BowenTemperatureMode = 0;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);

  reset();
  All.BowenTemperatureMode = 1;
  cells[1].E_rad = cells[2].E_rad = 0.0;
  for(i=0; i<IDORT_NUM_ANGLES; i++) cells[1].Intensity[i] = 0.0;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 2);
  CHECK(!diagnostic_in_parallel);
  CHECK(strstr(diagnostic_text, "finite zero radiation accepted; first report on this rank") != NULL);
  CHECK(strstr(diagnostic_text, "context=update_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "i=1 ID=10000000001") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=0 Volume=17") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=4 negative_angles=0 nonfinite_angles=0") != NULL);
  update_opacity(0);
  CHECK(diagnostic_calls == 2); /* No duplicate warning on subsequent updates. */
  refresh_bowen_opacity();
  CHECK(diagnostic_calls == 2); /* Shared suppression across both callers. */
  CHECK(terminate_calls == 0);

  /* Invalid radiation remains diagnostic and fatal after a zero warning. */
  reset();
  All.BowenTemperatureMode = 1;
  cells[1].E_rad = INFINITY;
  update_opacity(0);
  CHECK(terminate_calls == 1 && invalid_count == 1 && diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "invalid dust temperature") != NULL);
  CHECK(strstr(diagnostic_text, "i=1 ID=10000000001") != NULL);
  CHECK(!diagnostic_in_parallel && !terminated_in_parallel);
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
    functions = ("#define fprintf capture_fprintf\n#ifdef BOWEN_OPACITY\n"
                 + extract_function(source, "get_bowen_opacity") + "\n"
                 + extract_function(source, "report_bowen_cell")
                 + "\n#endif\n" + extract_function(source, "set_opacity") + "\n"
                 + extract_function(source, "update_opacity") + "\n#ifdef BOWEN_OPACITY\n"
                 + extract_function(source, "refresh_bowen_opacity") + "\n#endif\n#undef fprintf\n")
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
