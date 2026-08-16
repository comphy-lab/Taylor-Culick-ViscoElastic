/**
# Planar Newtonian Taylor--Culick retraction

Two-dimensional **planar** (not axisymmetric) retraction of a semi-infinite
Newtonian liquid sheet, used as an independent VOF cross-check of a
co-moving sharp-interface FEM (pyoomph) campaign against
Savva & Bush, *JFM* **626** (2009).

This case is derived from `simulationCases/TaylorCulick.c` with two
deliberate changes:

1. **Newtonian.** The scalar log-conformation solver
   (`log-conform-viscoelastic-scalar-2D.h`), the `two-phaseVE.h`
   `Gp`/`lambda` coupling, the conformation-tensor entries in
   `adapt_wavelet` and the polymer stress `acceleration` are all removed.
   Setting `Ec = 0` would already give a bit-exact Newtonian answer, but it
   still advects and relaxes the conformation tensor; here the machinery is
   gone, so nothing is computed for it.
2. **Planar.** `axi.h` is removed, so `cm = fm = 1` and the metric terms
   disappear from the momentum, viscous and VOF operators.

## Geometry and orientation

The sheet lies along $x$ and is symmetric about the midplane $y = 0$:

- $y = 0$ (bottom) is the **midplane**; Basilisk's default symmetry
  boundary condition is exactly the mirror condition required there.
- The liquid occupies $y < h_0/2$ for $x$ beyond the retracting edge, and
  the free edge is closed by a semicircular cap of radius $h_0/2$ centred
  on the midplane.
- The edge retracts towards $+x$, i.e. the gas region $x < x_{tip}$ grows.
- $x = L_0$ (right) is the quiescent far end of the semi-infinite sheet.

This is the source case rotated by ninety degrees; the axisymmetric hole of
initial radius `hole0` becomes a planar edge at initial position `xtip0`.

## Non-dimensionalisation

$\rho_l = 1$, $\sigma = 1$, $h_0 = 1$ (**full** thickness, half-thickness
$1/2$), so $\mu_l = Oh = \mu/\sqrt{\rho\sigma h_0}$ and the Taylor--Culick
speed is
$$V_{TC} = \sqrt{2\sigma/(\rho h_0)} = \sqrt{2}.$$
Savva & Bush use $Oh_{SB} = \mu/\sqrt{2 h_0 \rho \sigma}$, hence
$\mu = \sqrt{2}\,Oh_{SB}$, and their viscous clock is
$\tau_{vis} = \mu h_0/(2\sigma) = \mu/2$ with $t^* = t/\tau_{vis}$.

## Runtime parameters

Parameters are loaded from a `key=value` file through `src-local/params.h`,
exactly as in the parent case; the runner passes `case.params` as `argv[1]`.
*/

#include "navier-stokes/centered.h"

/**
The parent case writes `#define FILTERED` because `two-phaseVE.h` tests it
with `#ifdef`.  Basilisk's own `two-phase-generic.h` tests it with `#if`,
so it needs a value.  Both select the same vertex-averaged `sf` smearing.
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

`FERR`, `KERR` are the parent case values.  `VELERR` is a runtime parameter
because the parent value (`1e-6`) refines essentially the whole moving
region to `MAXlevel` and is not affordable over a domain that must be a
hundred sheet thicknesses long.
*/
#define FERR 1e-3
#define KERR 1e-6

double VELERR;

/**
## Outer boundary conditions

The bottom boundary is the midplane and keeps Basilisk's default symmetry
condition.  The left (gas), right (far end of the sheet) and top (gas)
boundaries are open: the quiescent far field of a semi-infinite sheet at
rest has zero pressure and zero normal velocity gradient.
*/
u.n[left] = neumann(0.);
p[left] = dirichlet(0.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);

int CaseNo, MAXlevel, MINlevel;
double tmax, Ldomain, xtip0, h0, tauvis;
char dumpFile[128], logFile[128], tipFile[128], snapshotFile[160];

/**
## Periodic-event increments must have a non-zero static initialiser

`qcc` registers events inside `_init_solver()`, which runs at the very top
of `main()` *before* any user statement.  `init_event()` classifies each
event expression there by calling it twice and watching whether `i`/`t`
change: an expression that leaves them unchanged is taken to be a
*condition*, not an *increment*.  The classification is done once
(`ev->nexpr` is zeroed) and is never redone.

So `event e (t = 0.; t += tsnap)` where `tsnap` is only assigned inside
`main()` is registered while `tsnap == 0.`, is misclassified as a
condition, and the event fires exactly once, at `t = 0`.  This is silent:
there is no warning and the run otherwise completes normally.

Giving these two variables a non-zero file-scope initialiser makes the
classification correct.  `ev->t` and the increment are then evaluated from
the *runtime* values (`init_event()` is re-run at `iter == 0`, after
`main()`), so `tsnap` and `tout` remain fully settable from `case.params`.

Absolute-time events such as `event stop_simulation (t = tmax)` are not
affected: `t = tmax` is still classified as an initialiser and `ev->t` is
recomputed from it at `iter == 0`.
*/
double tsnap = 1.;
double tout = 0.1;

/**
### main()
*/
int main (int argc, char const * argv[])
{
  params_init_from_argv(argc, argv);

  CaseNo = param_int("CaseNo", 1000);
  MAXlevel = param_int("MAXlevel", 13);
  MINlevel = param_int("MINlevel", 6);
  Ldomain = param_double("Ldomain", 128.);
  tmax = param_double("tmax", 40.);
  tsnap = param_double("tsnap", 0.5);
  tout = param_double("tout", 0.02);
  dtmax = param_double("dtmax", 1e-3);
  VELERR = param_double("VELERR", 1e-3);

  /**
  ## The timestep cap must be written to `DT`, not just `dtmax`

  `navier-stokes/centered.h` contains

      event set_dtmax (i++,last) dtmax = DT;

  so `dtmax` is re-read from the global `DT` at the top of *every* timestep.
  Assigning `dtmax` here therefore binds the first step only; from `i = 1`
  onwards the cap silently reverts to `DT`, whose default is effectively
  unbounded, and the timestep is limited by the CFL and surface-tension
  criteria alone.

  That matters because the first `tip_output` event forces the solver to land
  exactly on `t = tout`.  With no cap the solver takes a single step of
  `dt ~ 2e-2` straight after two steps of `dt ~ 1e-6` -- a jump of four
  orders of magnitude.  At `MAXlevel = 11` this is a survivable transient
  (`ke` spikes to ~1e-3 and decays), which is why it went unnoticed; at
  `MAXlevel = 13` with `Oh_SB = 1` the same step is far beyond the capillary
  and viscous stability limits, the implicit viscous solve fails to converge
  and `ke` reaches 5.5e6 by `i = 2`, i.e. the run is destroyed.

  Setting `DT` makes the cap persistent, which is what `dtmax=` in the
  parameter file was always meant to express.
  */
  DT = dtmax;

  h0 = param_double("h0", 1.);
  xtip0 = param_double("xtip0", 1.);

  rho1 = param_double("rho1", 1.);
  mu1 = param_double("mu1", 5e-2);
  rho2 = param_double("rho2", 1e-3);
  mu2 = param_double("mu2", 1e-5);

  /**
  `tauvis = mu1/2` is the Savva & Bush viscous time; it is only used to
  report `t* = t/tauvis` in the tip file.
  */
  tauvis = mu1/2.;

  if (CaseNo < 1000 || MAXlevel < 1 || MAXlevel > 20 ||
      MINlevel < 1 || MINlevel > MAXlevel || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || tout <= 0. || dtmax <= 0. ||
      dtmax > tmax || h0 <= 0. || xtip0 <= 0. ||
      xtip0 + h0 >= Ldomain || VELERR <= 0. ||
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

  if (pid() == 0) {
    fprintf(ferr,
            "PLANAR NEWTONIAN Taylor-Culick\n"
            "CaseNo=%d MAXlevel=%d MINlevel=%d Ldomain=%g "
            "tmax=%g tsnap=%g tout=%g dtmax=%g VELERR=%g\n",
            CaseNo, MAXlevel, MINlevel, Ldomain,
            tmax, tsnap, tout, dtmax, VELERR);
    fprintf(ferr,
            "liquid: rho=%g mu=%g ; gas: rho=%g mu=%g\n",
            rho1, mu1, rho2, mu2);
    fprintf(ferr,
            "h0=%g xtip0=%g V_TC=%g tau_vis=%g Delta_min=%g\n",
            h0, xtip0, sqrt(2.*f.sigma/(rho1*h0)), tauvis,
            Ldomain/(1 << MAXlevel));
  }

  run();
}

/**
## Initial interface

A flat sheet of full thickness `h0` (half-thickness `h0/2` above the
midplane `y = 0`) closed by a semicircular cap of radius `h0/2` whose
centre sits at `x = xtip0 + h0/2`.  The interface therefore meets the
midplane at `x = xtip0` and is everywhere normal to it.  The fluid starts
from rest.
*/
event init (t = 0)
{
  const double xc = xtip0 + h0/2.;

  if (!restore(file = dumpFile)) {
    /**
    Resolve the whole sheet band well enough that the flat interface is
    represented properly, and the cap neighbourhood at full resolution,
    *before* calling `fraction()`.  `adapt_wavelet` takes over from the
    first timestep.

    The `- Delta` terms are load-bearing.  `refine()` evaluates its
    condition at cell *centres*, so a plain `y < h0/2 + pad` test never
    fires on the initial coarse grid whenever `L0/2^MINlevel` is larger
    than about `h0`: the bottom row of cells has its centre above the
    band it is supposed to resolve, nothing is refined, and `fraction()`
    then writes a single smeared value (`f = 0.25` for `h0/2 = 0.5` in a
    cell of size 2) across the whole sheet.  The run still starts and
    still looks healthy, but the initial interface is wrong -- with
    `L0 = 128` and `MINlevel = 6` this put the tip at `x = 1.573`
    instead of `x = xtip0 = 1`.  Testing the cell's lower edge instead of
    its centre makes the pre-refinement independent of `MINlevel`.

    The padding is a fraction of `h0`, rather than a fixed length, so the
    resolved band stays proportionate when the sheet thickness changes.
    */
    const double pad = 0.1*h0;
    refine(y - Delta < h0/2. + pad && level < MAXlevel - 3);
    refine(x - Delta < xc + 2.*h0 && y - Delta < h0/2. + pad &&
           level < MAXlevel);
    fraction(f, x < xc
             ? sq(h0/2.) - (sq(x - xc) + sq(y))
             : h0/2. - y);

    /**
    Check the initial half-sheet area against the planar analytic geometry:
    the flat sheet contributes `(h0/2)(L0 - xc)`.  In the simulated
    half-domain, `x < xc` and `y >= 0` retain one quarter of the circular rim,
    which contributes `pi h0^2/16`.  This turns the smeared-interface
    failure into an immediate error instead of a plausible-looking run.
    */
    const double expected = (h0/2.)*(L0 - xc) + pi*sq(h0)/16.;
    double area = 0.;
    foreach (reduction(+:area))
      area += f[]*dv();

    if (pid() == 0)
      fprintf(ferr, "initial half-sheet area = %g (expected %g, "
              "relative error %.3g)\n",
              area, expected, fabs(area - expected)/expected);
    if (fabs(area - expected) > 0.05*expected) {
      fprintf(ferr, "ERROR: initial volume fraction is wrong; the sheet is "
              "probably unresolved on the initial grid.\n");
      /**
      Returning non-zero from an event stops the time loop but leaves the
      process status successful.  Fail the process so batch runners cannot
      accept a rejected initial condition as a completed run.
      */
      fflush(ferr);
      exit(1);
    }
  }
}

/**
## Adaptive mesh refinement

Interface, both velocity components and curvature.  The conformation
components of the parent case are gone.
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
## Tip tracking

Two independent estimates of the retracting-edge position are written, both
sub-cell accurate:

- `x_tip` is the smallest $x$ over the reconstructed VOF facets in the
  midplane band $y < h_0/10$.  For a convex cap this is the interface
  position on the midplane.
- `x_tip_vof` is the gas length along the bottom row of cells,
  $\int (1-f)\,\mathrm{d}x$ at $y \to 0$.  It agrees with `x_tip` while the
  midplane is crossed exactly once.

A disagreement between the two flags rim pinch-off or an entrained bubble
on the midplane, and is the reason both are recorded.
*/
event tip_output (t = 0.; t += tout)
{
  double xtip = HUGE, xtipglobal = HUGE, xvof = 0.;
  const double band = h0/10.;

  foreach (reduction(min:xtip) reduction(min:xtipglobal)
           reduction(+:xvof)) {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2)
        for (int k = 0; k < 2; k++) {
          double xf = x + segment[k].x*Delta;
          double yf = y + segment[k].y*Delta;
          if (xf < xtipglobal)
            xtipglobal = xf;
          if (yf < band && xf < xtip)
            xtip = xf;
        }
    }
    if (y < 0.75*Delta)
      xvof += (1. - clamp(f[], 0., 1.))*Delta;
  }

  if (pid() == 0) {
    FILE * fp = fopen(tipFile, t == 0. ? "w" : "a");
    if (!fp)
      fprintf(ferr, "tip_output: failed to open %s at t = %g, "
              "sample dropped\n", tipFile, t);
    else {
      if (t == 0.)
        fprintf(fp, "# planar Newtonian Taylor-Culick, CaseNo %d\n"
                    "# mu1 %g rho1 %g mu2 %g rho2 %g "
                    "MAXlevel %d Ldomain %g\n"
                    "# tau_vis %g V_TC %g\n"
                    "# t tstar x_tip x_tip_vof x_tip_global\n",
                CaseNo, mu1, rho1, mu2, rho2, MAXlevel, Ldomain,
                tauvis, sqrt(2.*f.sigma/(rho1*h0)));
      /**
      `xtip` and `xtipglobal` are reduction minima seeded at `HUGE`. If a
      band ever holds no interface -- after pinch-off, or if the rim leaves
      the midplane band -- the seed survives the reduction and would be
      written as a tip position of order `1e308`. Emit `nan` instead, so a
      gap in the record reads as missing rather than as a real number that
      would silently poison a downstream fit. `nan`, not Basilisk's
      `nodata`, because `nodata` is itself a large finite sentinel and would
      land in the data column looking like a coordinate; the axisymmetric
      counterpart does the same.
      */
      fprintf(fp, "%.8e %.8e %.8e %.8e %.8e\n",
              t, tauvis > 0. ? t/tauvis : 0.,
              xtip == HUGE ? nan("") : xtip, xvof,
              xtipglobal == HUGE ? nan("") : xtipglobal);
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

The kinetic energy is the planar (per unit span) integral; the parent
case's $2\pi y$ axisymmetric weight is removed.  The decay stop is guarded
by `t > 1.` because the planar sheet starts exactly at rest, so the parent
case's `i > 10` guard would trip on the initial transient.
*/
event log_writing (i++)
{
  double ke = 0.;
  int stop = 0;
  int dump_state = 0;
  foreach (reduction(+:ke))
    ke += (0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);

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
        fprintf(fp, "i dt t ke cells\n");
      }
      fprintf(fp, "%d %.8e %.8e %.8e %ld\n", i, dt, t, ke, grid->tn);
      fclose(fp);
      if (i % 100 == 0)
        fprintf(ferr, "%d %.8e %.8e %.8e %ld\n", i, dt, t, ke, grid->tn);

      if (!(ke == ke)) {
        fprintf(ferr, "ERROR: kinetic energy is NaN.\n");
        stop = 1;
        dump_state = 1;
      }
      if (ke > 1e3 && i > 10) {
        fprintf(ferr, "ERROR: kinetic energy blew up.\n");
        stop = 1;
        dump_state = 1;
      }
      if (ke < 1e-8 && t > 1.) {
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
