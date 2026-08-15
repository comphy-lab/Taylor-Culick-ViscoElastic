/**
# Axisymmetric Newtonian Taylor--Culick retraction

Axisymmetric hole-opening retraction of a Newtonian liquid film, derived from
the elastic case in [TaylorCulick.c](TaylorCulick.c) with the non-Newtonian
machinery removed rather than parked at zero: no scalar log-conformation
solver, no `two-phaseVE.h` `Gp`/`lambda` coupling, no conformation-tensor
components in `adapt_wavelet`. It is the axisymmetric counterpart of
[TaylorCulickPlanarNewtonian.c](TaylorCulickPlanarNewtonian.c), sharing its
Newtonian two-phase path, its three defect fixes, and its in-code tip
diagnostic adapted to the axisymmetric orientation.

## Phase convention

- `f = 1`: dense retracting film, `rho1`, `mu1`.
- `f = 0`: surrounding gas, `rho2`, `mu2`.

## Geometry and orientation

`axi.h` makes `y` the radial coordinate; the bottom boundary `y = 0` is the
symmetry axis. The film's own midplane symmetry is at `x = 0` (the flat
axial midplane of the sheet), which is Basilisk's default symmetry
condition and is left alone. A circular hole of initial radius `hole0`
opens about the axis; the film occupies `y > R(t)` and retracts outward as
`R(t)` grows. The top and right boundaries are open outflow.

## Runtime parameters

Parameters are loaded from a `key=value` file through `src-local/params.h`;
the runner passes the copied `case.params` file as `argv[1]`.
*/

#include "axi.h"
#include "navier-stokes/centered.h"

/**
`two-phase-generic.h` (included by `two-phase.h`) tests `FILTERED` with
`#if`, so it needs a value, not just a definition.
*/
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#define PARSE_PARAMS_IMPLEMENTATION
#include "params.h"
#undef PARSE_PARAMS_IMPLEMENTATION

/**
## Numerical controls

Hardcoded exactly as in the parent elastic case: these tolerances describe
the discretisation, not the constitutive model, so they stay fixed across a
refinement study.
*/
#define FERR 1e-3
#define VELERR 1e-6
#define KERR 1e-6

/**
## Outer boundary conditions

Top and right are open outflow. The bottom (`y = 0`) is the axis and keeps
`axi.h`'s built-in axis condition; the film's own midplane at `x = 0` keeps
Basilisk's default symmetry condition.
*/
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

int CaseNo, MAXlevel, MINlevel;
double tmax, Ldomain;
char dumpFile[128], logFile[128], tipFile[128], snapshotFile[160];

/**
## Periodic-event increments need a non-zero static initialiser

`qcc` registers events inside `_init_solver()`, which runs before any
statement of `main()`. `init_event()` classifies each event expression there
by calling it twice and watching whether `i` or `t` change; an expression
that leaves both unchanged is taken to be a *condition* rather than an
*increment*, and that classification is never redone.

So `event e (t = 0.; t += tsnap)` with `tsnap` assigned only inside `main()`
is registered while `tsnap == 0.`, is misclassified as a condition, and
fires exactly once, at `t = 0`. Nothing warns: the run completes normally
and `intermediate/` (and the tip file) hold a single row. A non-zero
file-scope initialiser fixes the classification; the runtime value is still
what gets used, because `init_event()` re-runs at `iter == 0`, after
`main()`, so `tsnap` stays settable from `case.params`.
*/
double tsnap = 1.;

/**
### main()
*/
int main (int argc, char const * argv[])
{
  params_init_from_argv(argc, argv);

  CaseNo = param_int("CaseNo", 1000);
  MAXlevel = param_int("MAXlevel", 12);
  MINlevel = param_int("MINlevel", max(6, MAXlevel - 4));
  Ldomain = param_double("Ldomain", 100.);
  tmax = param_double("tmax", 40.);
  tsnap = param_double("tsnap", 1.0);
  dtmax = param_double("dtmax", 1e-5);

  rho1 = param_double("rho1", 1.);
  mu1 = param_double("mu1", 5e-2);
  rho2 = param_double("rho2", 1e-3);
  mu2 = param_double("mu2", 1e-5);

  if (CaseNo < 1000 || MAXlevel < 1 || MAXlevel > 20 ||
      MINlevel < 1 || MINlevel > MAXlevel || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || dtmax <= 0. || dtmax > tmax ||
      rho1 <= 0. || rho2 <= 0. || mu1 < 0. || mu2 < 0.) {
    fprintf(ferr, "ERROR: invalid runtime parameters.\n");
    return 1;
  }

  L0 = Ldomain;
  X0 = 0.;
  Y0 = 0.;
  init_grid(1 << MINlevel);

  if (system("mkdir -p intermediate") != 0) {
    fprintf(ferr, "ERROR: unable to create intermediate output directory.\n");
    return 1;
  }

  sprintf(dumpFile, "dump");
  sprintf(logFile, "c%d-log", CaseNo);
  sprintf(tipFile, "c%d-tip.dat", CaseNo);

  f.sigma = 1.;
  TOLERANCE = 1e-4;
  CFL = 0.5;

  /**
  ## The timestep cap must be written to `DT`, not to `dtmax`

  `navier-stokes/centered.h` contains

  ~~~literatec
  event set_dtmax (i++,last) dtmax = DT;
  ~~~

  so `dtmax` is re-read from the global `DT` at the top of *every*
  timestep. Assigning `dtmax` here would bind the first step only; from
  `i = 1` onwards the cap silently reverts to `DT` (`HUGE` by default). The
  failure is resolution-dependent -- a coarse run survives and looks
  converged while the same script blows up at higher `MAXlevel` -- which is
  exactly backwards from how a refinement study should behave. Setting `DT`
  makes the cap persistent, which is what `dtmax=` in the parameter file is
  meant to express.

  ## The capillary constraint is `tension.h`'s job, and it does it

  `tension.h` contributes its own `stability` event that narrows `dtmax` to
  the capillary limit

  $$\Delta t_\sigma = \sqrt{\rho_m\Delta_{min}^3/(\pi\sigma)}$$

  using the actual finest cell carrying an interface, not a worst-case
  guess; this was checked against this solver stack with an instrumented
  `tension.h` and does apply the cap from `i = 0`. So the constraint is
  reported below rather than re-imposed here -- clamping `DT` to the
  `MAXlevel` value would pin every run to the finest cell the adaptation
  *could* produce even while the interface still sits on coarser cells.
  `dt_sigma` is printed at `MAXlevel` as the worst case the run can reach.
  */
  const double rho_m = (rho1 + rho2)/2.;
  const double Delta_min = Ldomain/(1 << MAXlevel);
  const double dt_sigma = sqrt(rho_m*cube(Delta_min)/(pi*f.sigma));
  DT = dtmax;

  if (pid() == 0) {
    fprintf(ferr,
            "AXISYMMETRIC NEWTONIAN Taylor-Culick\n"
            "CaseNo=%d MAXlevel=%d MINlevel=%d Ldomain=%g "
            "tmax=%g dtmax_requested=%g dt_sigma=%g DT=%g\n",
            CaseNo, MAXlevel, MINlevel, Ldomain, tmax, dtmax,
            dt_sigma, DT);
    fprintf(ferr,
            "liquid: rho=%g mu=%g ; gas: rho=%g mu=%g\n",
            rho1, mu1, rho2, mu2);
  }

  run();
}

/**
## Initial interface

The initial condition is the circular cap and retracting sheet used by the
elastic case: a hole of radius `hole0` about the axis, closed by a
semicircular rim of radius `h0/2`, opening onto a flat film of thickness
`h0/2` (from the axis to `x = h0/2`) that extends to `y = Ldomain`.
*/
event init (t = 0)
{
  const double hole0 = 1.;
  const double h0 = 1.;

  if (!restore(file = dumpFile)) {
    /**
    The `- Delta` terms are load-bearing. `refine()` evaluates its
    condition at cell *centres*, so a plain `x < h0/2 + pad` test never
    fires on the initial grid whenever `Ldomain/2^MINlevel` exceeds about
    `h0`: the first column of cells has its centre outside the band it is
    meant to resolve, nothing is refined, and `fraction()` then writes one
    smeared value across the whole sheet. Testing the cell's inner edge
    instead of its centre makes the pre-refinement independent of
    `MINlevel`.
    */
    const double pad = 0.1*h0;
    refine(x - Delta < h0/2. + pad && level < MAXlevel - 3);
    refine(x - Delta < h0/2. + pad &&
           y - Delta < hole0 + 2.*h0 && level < MAXlevel);
    fraction(f, y < hole0 + h0/2.
             ? sq(h0/2.) - (sq(x) + sq(y - h0/2. - hole0))
             : h0/2. - x);

    /**
    ## Fail loudly on a wrong initial condition

    In `axi.h` the metric is `cm = y`, so `sum f dv()` is the liquid volume
    divided by `2*pi`. For the flat sheet `0 < x < h0/2`,
    `hole0 + h0/2 < y < L0`, that is
    `(h0/2)(L0^2 - (hole0 + h0/2)^2)/2`.

    The rim is a QUARTER disc, not a half one. The circle of radius `h0/2`
    about `(0, yrim)` is cut twice: by the midplane symmetry at `x = 0`, and
    by the `y < yrim` branch of the `fraction()` expression above, since for
    `y >= yrim` the flat-sheet branch already covers that band at full
    thickness. Its area is therefore `pi h0^2/16`, and the centroid of a
    quarter disc of radius `a` lying below `yrim` sits at `yrim - 4a/(3 pi)`,
    not at `yrim`, so by Pappus the contribution is
    `(pi h0^2/16)(yrim - 4(h0/2)/(3 pi))`.

    Getting this wrong is how the check read `relative error 1.35e-4` for
    every run at every resolution: the half-disc form overstates the rim by
    `0.336`, which is 0.013% of the total and so sat quietly inside the 5%
    tolerance while looking like a small mesh error. With the quarter-disc
    form the residual is the discretisation error and nothing else, which is
    what makes this a guard rather than a decoration. The tolerance stays at
    5% deliberately: it exists to catch a sheet smeared across one cell on a
    coarse initial grid, not to police quadrature.
    */
    const double yrim = hole0 + h0/2.;
    const double arim = h0/2.;
    const double expected = (h0/2.)*(sq(L0) - sq(yrim))/2.
      + (pi*sq(h0)/16.)*(yrim - 4.*arim/(3.*pi));
    double vol = 0.;
    foreach (reduction(+:vol))
      vol += f[]*dv();

    if (pid() == 0)
      fprintf(ferr, "initial liquid volume/(2 pi) = %g (expected %g, "
              "relative error %.3g)\n",
              vol, expected, fabs(vol - expected)/expected);
    if (fabs(vol - expected) > 0.05*expected) {
      fprintf(ferr, "ERROR: initial volume fraction is wrong; the sheet is "
              "probably unresolved on the initial grid.\n");
      /**
      `return 1` from an event only stops the time loop; the process still
      exits 0 and a batch runner or systemd unit records success. A guard
      that reports a wrong initial condition as a clean early finish is
      worse than no guard, so fail the process itself.
      */
      fflush(ferr);
      exit(1);
    }
  }
}

/**
## Adaptive mesh refinement

Interface, both velocity components and curvature -- the conformation
components of the parent elastic case are gone.
*/
scalar KAPPA[];

event adapt_mesh (i++)
{
  curvature(f, KAPPA);
  adapt_wavelet((scalar *) {f, u.x, u.y, KAPPA},
                (double[]) {FERR, VELERR, VELERR, KERR},
                MAXlevel);
}

/**
## Tip (hole-radius) tracking

The axisymmetric counterpart of `TaylorCulickPlanarNewtonian.c`'s
`tip_output` event, with the roles of `x` and `y` swapped: here `y` is the
radial direction the hole grows in, and `x = 0` is the film's midplane, so
the near-midplane band is `x < h0/10` rather than `y < h0/10`.

Two independent estimates of the hole radius are written, both sub-cell
accurate:

- `R_tip` is the smallest `y` over the reconstructed VOF facets in the
  near-midplane band `x < h0/10` (the film occupies `y > R(t)`, so the
  interface crosses the midplane at `y = R`).
- `R_tip_vof` is the gas length along the near-axis-midplane column,
  `integral (1-f) dy` at `x -> 0`, an independent cross-check that agrees
  with `R_tip` while the midplane is crossed exactly once.

A disagreement between the two flags rim pinch-off or an entrained bubble
on the midplane, which is why both are recorded, alongside `R_tip_global`,
the same minimum taken over the whole interface.
*/
event tip_output (t = 0.; t += tsnap)
{
  const double h0 = 1.;
  double Rtip = HUGE, Rtipglobal = HUGE, Rvof = 0.;
  const double band = h0/10.;

  foreach (reduction(min:Rtip) reduction(min:Rtipglobal)
           reduction(+:Rvof)) {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2)
        for (int k = 0; k < 2; k++) {
          double xf = x + segment[k].x*Delta;
          double yf = y + segment[k].y*Delta;
          if (yf < Rtipglobal)
            Rtipglobal = yf;
          if (xf < band && yf < Rtip)
            Rtip = yf;
        }
    }
    if (x < 0.75*Delta)
      Rvof += (1. - clamp(f[], 0., 1.))*Delta;
  }

  if (pid() == 0) {
    FILE * fp = fopen(tipFile, t == 0. ? "w" : "a");
    if (fp) {
      if (t == 0.)
        fprintf(fp, "# axisymmetric Newtonian Taylor-Culick, CaseNo %d\n"
                    "# mu1 %g rho1 %g mu2 %g rho2 %g "
                    "MAXlevel %d Ldomain %g\n"
                    "# t R_tip R_tip_vof R_tip_global\n",
                CaseNo, mu1, rho1, mu2, rho2, MAXlevel, Ldomain);
      /**
      `Rtip` and `Rtipglobal` are reduction minima seeded at `HUGE`. If a
      band ever contains no interface -- after pinch-off, or if the rim
      leaves the midplane band entirely -- the seed survives the reduction
      and would be written as a radius of order `1e308`. Emit `nan` instead,
      so a gap in the record reads as missing rather than as a real number
      that would silently poison any downstream fit or plot.
      */
      fprintf(fp, "%.8e %.8e %.8e %.8e\n", t,
              Rtip == HUGE ? nan("") : Rtip, Rvof,
              Rtipglobal == HUGE ? nan("") : Rtipglobal);
      fclose(fp);
    }
  }
}

/**
## Restart and snapshots
*/
event writing_files (t = 0.; t += tsnap)
{
  dump(file = dumpFile);
  sprintf(snapshotFile, "intermediate/snapshot-%5.4f", t);
  dump(file = snapshotFile);
}

/**
## Progress log

The logged kinetic energy is the axisymmetric integral, unchanged from the
parent elastic case. A run is stopped early only for the same two
deterministic failure conditions: energy blow-up or complete decay after
the initial transient.
*/
event log_writing (i++)
{
  double ke = 0.;
  int stop = 0;
  int dump_state = 0;
  foreach (reduction(+:ke))
    ke += (2.*pi*y)*(0.5*rho(f[])*
                     (sq(u.x[]) + sq(u.y[])))*sq(Delta);

  if (pid() == 0) {
    FILE * fp = fopen(logFile, i == 0 ? "w" : "a");
    if (!fp) {
      fprintf(ferr, "ERROR: cannot open %s\n", logFile);
      stop = 1;
    }
    else {
      if (i == 0) {
        fprintf(fp, "CaseNo %d, MAXlevel %d, mu1 %g, Ldomain %g\n",
                CaseNo, MAXlevel, mu1, Ldomain);
        fprintf(fp, "i dt t ke\n");
      }
      fprintf(fp, "%d %.8e %.8e %.8e\n", i, dt, t, ke);
      fclose(fp);
      fprintf(ferr, "%d %.8e %.8e %.8e\n", i, dt, t, ke);

      // The template's 1e2 threshold is wrong for this case: the axisymmetric
      // rim mass grows like R^2, so the physical kinetic energy grows without
      // bound (~2e4 by t = 40 at mu1 = 0.05) and crossed 1e2 at t = 6.4 on a
      // perfectly healthy run (dt steady, R(t) smooth). A genuine blow-up in
      // this family reaches 1e6+ within a few steps (recorded defect history:
      // ke hit 5.5e6/3.2e8 by i = 2). Guard at 1e6.
      if (ke > 1e6 && i > 10) {
        fprintf(ferr, "ERROR: kinetic energy blew up.\n");
        stop = 1;
        dump_state = 1;
      }
      if (ke < 1e-8 && i > 10) {
        fprintf(ferr, "Kinetic energy decayed below the stopping threshold.\n");
        stop = 1;
        dump_state = 1;
      }
    }
  }

  mpi_all_reduce(stop, MPI_INT, MPI_MAX);
  mpi_all_reduce(dump_state, MPI_INT, MPI_MAX);
  if (dump_state)
    dump(file = dumpFile);
  assert(ke > -1e-10);
  return stop;
}

/**
## Completion
*/
event stop_simulation (t = tmax)
{
  if (pid() == 0)
    fprintf(ferr, "Case %d complete at t=%g.\n", CaseNo, t);
  return 1;
}
