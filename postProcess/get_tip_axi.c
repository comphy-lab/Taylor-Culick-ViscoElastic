/**
# Snapshot-based axisymmetric tip cell velocities

Restores one snapshot from `TaylorCulickAxiNewtonian.c` and recomputes the
same `R_tip` as that case's `tip_output` event (smallest facet `y` in the
near-midplane band `x < h0/10`). Then records the **cell** radial velocity
`u.y` in three places, none of which is a time derivative of `R_tip`:

- `u_tip`: `u.y[]` in the mixed cell that owns that facet.
- `u_gas1`, `u_gas2`: `u.y` interpolated one and two local cell widths into
  the hole (`y` decreasing, the gas), at the same axial station as the tip
  facet.

Must be compiled with `axi.h` so the restored dump's axisymmetric metrics
match the run.

Usage: `get_tip_axi snapshot-file [h0]`

Prints one line: `t R_tip u_tip u_gas1 u_gas2 Delta`
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];

int main (int argc, char const * argv[])
{
  if (argc < 2) {
    fprintf (stderr, "usage: %s snapshot-file [h0]\n", argv[0]);
    return 1;
  }
  double h0 = argc > 2 ? atof (argv[2]) : 1.;
  if (!restore (file = argv[1])) {
    fprintf (stderr, "%s: cannot restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
  f.prolongation = fraction_refine;

  double Rtip = HUGE, xf_tip = 0., yf_tip = 0.;
  double utip = nodata, Delta_tip = nodata;
  const double band = h0/10.;

  foreach() {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
        for (int k = 0; k < 2; k++) {
          double xf = x + segment[k].x*Delta;
          double yf = y + segment[k].y*Delta;
          if (xf < band && yf < Rtip) {
            Rtip = yf;
            xf_tip = xf;
            yf_tip = yf;
            utip = u.y[];
            Delta_tip = Delta;
          }
        }
    }
  }

  double ugas1 = nodata, ugas2 = nodata;
  if (Rtip < HUGE/2. && Delta_tip > 0.) {
    double y1 = yf_tip - Delta_tip;
    double y2 = yf_tip - 2.*Delta_tip;
    if (y1 > 0.)
      ugas1 = interpolate (u.y, xf_tip, y1);
    if (y2 > 0.)
      ugas2 = interpolate (u.y, xf_tip, y2);
  }

  fprintf (stdout, "%.8e %.8e %.8e %.8e %.8e %.8e\n",
           t,
           Rtip == HUGE ? nan ("") : Rtip,
           utip,
           ugas1,
           ugas2,
           Delta_tip);
  return 0;
}
