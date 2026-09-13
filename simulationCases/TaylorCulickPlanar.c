/**
# Planar Taylor--Culick retraction

Two-dimensional **planar** retraction of a semi-infinite liquid sheet, the
sibling of the axisymmetric hole-opening case in
[TaylorCulick.c](TaylorCulick.c).  It uses the same CoMPhy scalar
log-conformation solver and the same `two-phaseVE.h` phase mapping, so it
supports the same elastic and viscoelastic physics through `G1`, `lambda1`,
`G2` and `lambda2`, and reduces to a Newtonian sheet by parameter
(`Ec = 0`, equivalently `G1 = G2 = 0`) rather than by deleting code.

The *only* structural differences from the axisymmetric case are:

- no `#include "axi.h"`, so `cm = fm = 1` and the metric terms drop out of
  the momentum, viscous and VOF operators;
- a planar initial condition and matching boundary conditions;
- the conformation tensor has no $\theta\theta$ component, so `AThTh` is
  absent from `adapt_wavelet` (the solver already guards it with `#if AXI`);
- the kinetic-energy integral loses its $2\pi y$ weight;
- an in-code tip diagnostic (see below).

## Geometry and orientation

The sheet lies along $x$ and is symmetric about the midplane $y = 0$:

- $y = 0$ (bottom) is the **midplane**.  Basilisk's default symmetry
  boundary condition is exactly the mirror condition required there, so it
  is deliberately left alone.
- The liquid occupies $y < h_0/2$ beyond the retracting edge, and the free
  edge is closed by a semicircular cap of radius $h_0/2$ centred on the
  midplane, so the interface meets the midplane normally.
- The edge retracts towards $+x$: the gas region $x < x_{tip}$ grows.
- $x = L_0$ (right) is the quiescent far end of the semi-infinite sheet.

The axisymmetric case's initial hole radius `hole0` becomes the planar
initial edge position `xtip0`.  Keep `Ldomain` large enough that the edge
never approaches the right boundary: the edge travels at most
$\sqrt{2}\,t$, so `Ldomain` should exceed `xtip0 + sqrt(2)*tmax` with a
comfortable margin for the capillary waves that run ahead of it.

## Non-dimensionalisation

Lengths are scaled with the **full** sheet thickness $h_0$, so the
half-thickness is $h_0/2$; densities with $\rho_l$; and stresses with
$\sigma/h_0$.  Setting $\rho_l = \sigma = h_0 = 1$ gives

$$\mu_l = \frac{\mu}{\sqrt{\rho\sigma h_0}}, \qquad
  V_{TC} = \sqrt{\frac{2\sigma}{\rho h_0}} = \sqrt{2}.$$

Savva & Bush (*JFM* **626**, 2009) instead use the half-thickness in their
Ohnesorge number,
$$Oh = \frac{\mu}{\sqrt{2 h_0 \rho \sigma}} = \frac{\mu_l}{\sqrt{2}},$$
which is this project's only $Oh$. The coded `mu1` is $\mu_l$, NOT $Oh$; the
two differ by $\sqrt{2}$. A case at a quoted $Oh$ is run with
`mu1 = sqrt(2)*Oh`.
Their viscous time is $\tau_{vis} = \mu h_0 / (2\sigma) = \mu_l/2$ and their
reduced time is $t^* = t/\tau_{vis}$; both are reported in the tip file so
results can be plotted on either clock.

## Runtime parameters

Parameters are loaded from a `key=value` file through `src-local/params.h`.
The root runner passes the copied `case.params` file to the executable.
*/

#include "navier-stokes/centered.h"

#define FILTERED
#include "log-conform-viscoelastic-scalar-2D.h"
#include "two-phaseVE.h"

#include "navier-stokes/conserving.h"
#include "tension.h"
#define PARSE_PARAMS_IMPLEMENTATION
#include "params.h"
#undef PARSE_PARAMS_IMPLEMENTATION

/**
## Numerical controls

The tolerances are deliberately explicit so that a refinement study can
change them in one place without changing the constitutive model.
`VELERR` is a runtime parameter here because the planar domain is long and
the useful value depends strongly on `Ldomain`.
*/
#define FERR 1e-3
#define AERR 1e-4
#define KERR 1e-6

double VELERR = 1e-3;

/**
## Outer boundary conditions

The bottom boundary is the midplane and keeps the default symmetry
condition.  The left (gas), right (far end of the sheet) and top (gas)
boundaries are open: the quiescent far field of a sheet at rest carries no
curvature, hence zero pressure, and no normal velocity gradient.
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
## Periodic-event increments need a non-zero static initialiser

`qcc` registers events inside `_init_solver()`, which runs at the very top
of `main()`, *before* any user statement.  `init_event()` classifies each
event expression there by calling it twice and watching whether `i` or `t`
change; an expression that leaves both unchanged is taken to be a
*condition* rather than an *increment*.  That classification happens once
and is never redone.

So `event e (t = 0.; t += tsnap)` where `tsnap` is assigned only inside
`main()` is registered while `tsnap == 0.`, is misclassified as a
condition, and fires exactly once, at `t = 0`.  Nothing warns; the run
completes normally with a single snapshot.

A non-zero file-scope initialiser makes the classification correct.  The
runtime value is still what gets used, because `init_event()` re-runs at
`iter == 0` -- after `main()` -- so `tsnap` and `tout` remain fully
settable from `case.params`.

Absolute-time events such as `event stop_simulation (t = tmax)` are
unaffected: `t = tmax` is still classified as an initialiser and `ev->t` is
recomputed from it at `iter == 0`.
*/
double tsnap = 1.;
double tout = 0.1;

/**
### main()

Loads runtime parameters, assigns both phase rheologies, and starts the
Basilisk event loop.
*/
int main (int argc, char const * argv[])
{
  params_init_from_argv(argc, argv);

  CaseNo = param_int("CaseNo", 1000);
  MAXlevel = param_int("MAXlevel", 12);
  MINlevel = param_int("MINlevel", 6);
  Ldomain = param_double("Ldomain", 128.);
  tmax = param_double("tmax", 25.);
  tsnap = param_double("tsnap", 0.5);
  tout = param_double("tout", 0.02);
  dtmax = param_double("dtmax", 5e-3);
  VELERR = param_double("VELERR", 1e-3);

  h0 = param_double("h0", 1.);
  xtip0 = param_double("xtip0", 1.);

  rho1 = param_double("rho1", 1.);
  mu1 = param_double("mu1", 5e-2);
  rho2 = param_double("rho2", 1e-3);
  mu2 = param_double("mu2", 1e-5);

  G1 = param_double("G1", param_double("Ec", 0.));
  lambda1 = param_double("lambda1", param_double("De", 0.));
  G2 = param_double("G2", 0.);
  lambda2 = param_double("lambda2", 0.);
  TOLelastic = param_double("TOLelastic", 1e-2);

  /**
  Savva & Bush's viscous clock, reported alongside `t` in the tip file.
  */
  tauvis = mu1/2.;

  if (CaseNo < 1000 || MAXlevel < 1 || MAXlevel > 20 ||
      MINlevel < 1 || MINlevel > MAXlevel || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || tout <= 0. || dtmax <= 0. ||
      dtmax > tmax || h0 <= 0. || xtip0 <= 0. ||
      xtip0 + h0 >= Ldomain || VELERR <= 0. ||
      rho1 <= 0. || rho2 <= 0. || mu1 < 0. || mu2 < 0. ||
      G1 < 0. || G2 < 0. || lambda1 < 0. || lambda2 < 0. ||
      TOLelastic < 0. || TOLelastic >= 0.5) {
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
  timestep.  Assigning `dtmax` here would bind the first step only; from
  `i = 1` onwards the cap silently reverts to `DT`.  The failure is
  resolution-dependent -- a coarse run survives and looks converged while
  the same script blows up at higher `MAXlevel` -- which is exactly
  backwards from how a refinement study should behave.  Setting `DT` makes
  the cap persistent, which is what `dtmax=` in the parameter file is
  meant to express.

  ## The capillary constraint is `tension.h`'s job, and it does it

  `tension.h` contributes its own `stability` event that narrows `dtmax`
  to the capillary limit

  $$\Delta t_\sigma = \sqrt{\rho_m\Delta_{min}^3/(\pi\sigma)}$$

  using the *actual* finest cell carrying an interface, not a worst-case
  guess.  That narrowing was checked against this exact solver stack with
  an instrumented `tension.h`: at `i = 0` it sees finite
  `dmin`, `rho_m` and `sigma` and does apply the cap.  So the constraint
  is reported below rather than re-imposed here -- clamping `DT` to the
  `MAXlevel` value would pin every run to the finest cell the adaptation
  *could* produce even while the interface still sits on coarser cells.

  `dt_sigma` is printed at `MAXlevel` as the worst case the run can reach,
  which is what to compare an achieved `dt` against when a high-resolution
  run misbehaves.
  */
  const double rho_m = (rho1 + rho2)/2.;
  const double Delta_min = Ldomain/(1 << MAXlevel);
  const double dt_sigma = sqrt(rho_m*cube(Delta_min)/(pi*f.sigma));
  DT = dtmax;

  if (pid() == 0) {
    fprintf(ferr,
            "PLANAR Taylor-Culick\n"
            "CaseNo=%d MAXlevel=%d MINlevel=%d Ldomain=%g "
            "tmax=%g tsnap=%g tout=%g dtmax_requested=%g VELERR=%g\n",
            CaseNo, MAXlevel, MINlevel, Ldomain,
            tmax, tsnap, tout, dtmax, VELERR);
    fprintf(ferr,
            "phase1: rho=%g mu=%g G=%g lambda=%g; "
            "phase2: rho=%g mu=%g G=%g lambda=%g\n",
            rho1, mu1, G1, lambda1, rho2, mu2, G2, lambda2);
    fprintf(ferr,
            "TOLelastic=%g h0=%g xtip0=%g V_TC=%g tau_vis=%g Delta_min=%g "
            "dt_sigma=%g DT=%g\n",
            TOLelastic, h0, xtip0, sqrt(2.*f.sigma/(rho1*h0)), tauvis,
            Delta_min, dt_sigma, DT);
  }

  run();
}

/**
## Initial interface

A flat sheet of full thickness `h0` (half-thickness `h0/2` above the
midplane `y = 0`) closed by a semicircular cap of radius `h0/2` centred at
`x = xtip0 + h0/2`.  The interface therefore meets the midplane at
`x = xtip0`.  The fluid starts from rest.
*/
event init (t = 0)
{
  const double xc = xtip0 + h0/2.;

  if (!restore(file = dumpFile)) {
    /**
    Resolve the sheet band before calling `fraction()`, and the cap
    neighbourhood at full resolution; `adapt_wavelet` takes over from the
    first timestep.

    The `- Delta` terms are load-bearing.  `refine()` evaluates its
    condition at cell *centres*, so a plain `y < h0/2 + 0.1` test never
    fires on the initial grid whenever `Ldomain/2^MINlevel` exceeds about
    `h0`: the bottom row of cells has its centre above the band it is
    meant to resolve, nothing is refined, and `fraction()` then writes one
    smeared value across the whole sheet.  The run still starts and still
    looks healthy, but the initial interface is wrong.  With
    `Ldomain = 128` and `MINlevel = 6` this put the edge at `x = 1.573`
    instead of `x = xtip0 = 1`.  Testing the cell's lower edge instead of
    its centre makes the pre-refinement independent of `MINlevel`.

    The padding beyond `h0/2` is a fraction of `h0` itself (not a fixed
    `0.1`), so the band stays proportionate to the sheet thickness whether
    `h0` is 1 or 0.1.
    */
    const double pad = 0.1*h0;
    refine(y - Delta < h0/2. + pad && level < MAXlevel - 3);
    refine(x - Delta < xc + 2.*h0 && y - Delta < h0/2. + pad &&
           level < MAXlevel);
    fraction(f, x < xc
             ? sq(h0/2.) - (sq(x - xc) + sq(y))
             : h0/2. - y);

    /**
    ## Fail loudly on a wrong initial condition

    The smeared-initial-condition failure described above costs one wasted
    campaign precisely because it is silent, so the initial state is
    checked against its analytic value before the first timestep rather
    than trusted.

    The metric is unity here, so `sum f dv()` is the simulated half-sheet
    area: the flat part `0 < y < h0/2`, `xc < x < L0` contributes
    `(h0/2)(L0 - xc)` and the semicircular cap of radius `h0/2` adds
    `pi h0^2/8`.

    The tolerance is loose enough that ordinary VOF discretisation error
    never trips it and tight enough that the smearing failure -- which put
    the edge at `x = 1.573` instead of `x = 1` -- cannot get past.
    */
    const double expected = (h0/2.)*(L0 - xc) + pi*sq(h0)/8.;
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
      return 1;
    }
  }
}

/**
## Adaptive mesh refinement

Interface, both velocity components, curvature and the planar conformation
components.  There is no $\theta\theta$ component in planar geometry.
*/
scalar KAPPA[];

event adapt_mesh (i++)
{
  curvature(f, KAPPA);
  adapt_wavelet((scalar *) {f, u.x, u.y, KAPPA, A11, A12, A22},
                (double[]) {FERR, VELERR, VELERR, KERR, AERR, AERR, AERR},
                MAXlevel);
}

/**
## Tip diagnostic

The retraction speed is the quantity of interest, so it is measured in the
run rather than reconstructed offline.  Two independent, sub-cell-accurate
estimates of the edge position are written:

- `x_tip`: the smallest $x$ over the reconstructed VOF facets inside the
  midplane band $y < h_0/10$.  For a convex rim this is the interface
  position on the midplane.
- `x_tip_vof`: the gas length along the bottom row of cells,
  $\int (1-f)\,\mathrm{d}x$ at $y \to 0$.

They agree while the midplane is crossed exactly once; a growing
discrepancy is a useful flag for rim pinch-off or an entrained bubble on
the midplane.  `x_tip_global` is the same minimum taken over the whole
interface and catches anything that escapes the band.

The retraction speed is `d(x_tip)/dt`; differentiating in post-processing
(see `postProcess/tip_to_csv.py`) keeps this event cheap and lets the
smoothing window be chosen after the fact.
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
    if (fp) {
      if (t == 0.)
        fprintf(fp, "# planar Taylor-Culick, CaseNo %d\n"
                    "# mu1 %g rho1 %g G1 %g lambda1 %g "
                    "mu2 %g rho2 %g\n"
                    "# MAXlevel %d Ldomain %g tau_vis %g V_TC %g\n"
                    "# t tstar x_tip x_tip_vof x_tip_global\n",
                CaseNo, mu1, rho1, G1, lambda1, mu2, rho2,
                MAXlevel, Ldomain, tauvis,
                sqrt(2.*f.sigma/(rho1*h0)));
      fprintf(fp, "%.8e %.8e %.8e %.8e %.8e\n",
              t, tauvis > 0. ? t/tauvis : 0.,
              xtip, xvof, xtipglobal);
      fclose(fp);
    }
    else
      fprintf(stderr, "tip_output: failed to open %s at t = %g, "
                       "sample dropped\n", tipFile, t);
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

The kinetic energy is the planar (per unit span) integral; the
axisymmetric $2\pi y$ weight is absent.  The decay stop is guarded by
`t > 1.` rather than `i > 10` because the planar sheet starts exactly at
rest, so an iteration-count guard trips on the initial transient.
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
        fprintf(fp, "CaseNo %d, MAXlevel %d, mu1 %g, G1 %g, lambda1 %g, "
                    "Ldomain %g\n",
                CaseNo, MAXlevel, mu1, G1, lambda1, Ldomain);
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
