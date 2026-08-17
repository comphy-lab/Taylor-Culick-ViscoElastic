# Taylor-Culick-ViscoElastic

Taylor--Culick retraction cases in Basilisk C, organised using the CoMPhy
project layout. There are two geometries -- axisymmetric hole opening and a
planar semi-infinite sheet -- sharing one constitutive path: the scalar
2D/axisymmetric log-conformation solver from
[MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow).

## Repository layout

```
├── simulationCases/ - Active Basilisk entry points and generated case folders
│   ├── TaylorCulick.c - Axisymmetric elastic and viscoelastic hole-opening case
│   └── TaylorCulickPlanar.c - Planar semi-infinite sheet, same constitutive path
├── src-local/ - Stable rheology headers and typed parameter parser
│   ├── log-conform-viscoelastic-scalar-2D.h - Upstream scalar log-conformation solver
│   ├── two-phaseVE.h - Per-phase VE material-property coupling
│   ├── parse_params.h - key=value file loader
│   └── params.h - typed runtime accessors
├── postProcess/ - Basilisk-native and offline post-processing
│   ├── tip_to_csv.py - In-code tip log to t,tstar,x_tip,v_tip,v_over_VTC CSV
│   ├── get_facets.c - Interface segments from a snapshot (output_facets)
│   ├── get_tip.c - Independent tip position/velocity, recomputed from a
│   │   snapshot rather than read from the in-code log
│   ├── get_fields.c - f/u/KAPPA/PHI (dissipation) (+ A11..T22 under
│   │   -DVISCOELASTIC=1) interpolated onto a grid from a snapshot
│   ├── run_postprocess.sh - Compiles the three above and runs them over a
│   │   case directory's intermediate/snapshot-* files
│   └── plot_fields.py - Renders per-snapshot PNG frames and an mp4
├── scripts/params.sh - Shared shell parameter helpers
├── default.params - Purely elastic axisymmetric default (lambda1 = 1e30)
├── default-viscoelastic.params - Finite-relaxation axisymmetric default
├── default-planar.params - Planar Newtonian default (Oh_SB = 0.1)
├── sweep.params - Two-combination sweep definition
├── runSimulation.sh - Single-case compile/run driver
├── runParameterSweep.sh - Deterministic Cartesian sweep driver
├── AGENTS.md - Repository guidance
└── README.md - Project documentation
```

Generated case directories such as simulationCases/c1000/ are ignored and
contain the copied source, case.params, executable, logs, dumps, and snapshots
for that run.

## Requirements

- Basilisk C with qcc available in PATH
- bash, awk, and standard POSIX utilities

The local source snapshot was copied from MultiRheoFlow commit
4695a434f0750c8476c9e094bb8561384092299a on 2026-08-06. The scalar solver
matches upstream blob cea628698212c4172780b35ec8949fbfb1bc4570 exactly.
The companion two-phaseVE.h was whitespace-normalized from upstream blob
537123b1390f798ad7262321f7505b8c83efa1db.

## Single-case runs

Compile and run the elastic default:

```bash
bash runSimulation.sh --input default.params
```

Run the finite-relaxation case:

```bash
bash runSimulation.sh --input default-viscoelastic.params
```

Compile without running:

```bash
bash runSimulation.sh --input default.params --compile-only
```

Run the planar case:

```bash
bash runSimulation.sh --case simulationCases/TaylorCulickPlanar.c \
  --input default-planar.params
```

The executable receives case.params as argv[1]; this keeps each output
directory self-describing and restartable.

`--outdir DIR` puts the run directory somewhere other than
`simulationCases/c<CaseNo>/`, which is what you want when run data must stay
out of the checkout. `--openmp` adds `-fopenmp`, so the run honours
`OMP_NUM_THREADS`.

### Always build parallel: the serial build traps on the first step

Basilisk arms floating-point exception trapping in `grid/config.h` under

```c
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP
```

and `qcc` passes `-D_GNU_SOURCE` itself, so a plain serial build of these
cases has the trap armed and aborts with `SIGFPE` on the first timestep. An
OpenMP build defines `_OPENMP` and compiles `enable_fpe()` out, which is why
the trap has never shown up locally. Use `--openmp` or `--mpi`; the runner
warns if you ask for neither.

For MPI:

```bash
bash runSimulation.sh --case simulationCases/TaylorCulickPlanar.c \
  --input default-planar.params --mpi --np 48
```

`--mpi` compiles with `-D_MPI=1 -D_DEFAULT_SOURCE` and points `$CC99` at
`mpicc`. `-D_DEFAULT_SOURCE` still exposes the GNU features Basilisk needs
(`MAP_ANONYMOUS`, `MADV_DONTNEED`) while leaving `enable_fpe()` a no-op:
`qcc` strips `-D_GNU_SOURCE` from `$CC99` for MPI builds (`qcc.c:415`), so
the resulting binary has the same FP behaviour as the trusted OpenMP one.
Counting the symbol confirms the parity rather than assuming it:

```
nm <openmp-build>              | grep -c feenableexcept   ->  0
nm <mpi-build -D_DEFAULT_SOURCE> | grep -c feenableexcept ->  0
nm <mpi-build -D_GNU_SOURCE>     | grep -c feenableexcept ->  1
```

Never pass `-D_GNU_SOURCE` to an MPI build. `qcc` refuses to combine MPI with
OpenMP (*"OpenMP cannot be used with MPI (yet)"*), so an MPI build has no
`_OPENMP` to disarm the trap.

Omit `--np` to compile and launch a single rank; on a scheduler, use
`--compile-only` and submit the binary yourself.

## Three silent failures the cases now guard against

All three produced a run that completed normally and looked plausible. They
are written up in the case sources next to the code that avoids them; this is
the short version.

**`dtmax=` in a parameter file binds `DT`, not `dtmax`.**
`navier-stokes/centered.h` contains `event set_dtmax (i++,last) dtmax = DT;`,
so `dtmax` is re-read from the global `DT` at the top of every timestep and
an assignment in `main()` survives step 0 only. `DT` defaults to `HUGE`
(`utils.h:8`), so from `i = 1` the requested cap was simply gone. The trap is
that it is resolution-dependent: a coarse run survives and looks converged
while the same script blows up at higher `MAXlevel`, which is backwards from
how a refinement study should behave. Both cases now assign `DT`.

**`refine()` tests cell centres, so the initial sheet could be smeared.**
When `Ldomain/2^MINlevel` exceeds about `h0`, the first row of cells has its
centre outside the band the pre-refinement is meant to resolve, nothing is
refined, and `fraction()` writes one smeared value across the whole sheet.
Both cases now test the cell's near edge (`y - Delta < ...`) instead, which
makes the pre-refinement independent of `MINlevel`, and both check the
initial volume against its analytic value and abort if it is off by more
than 5%. Measured relative error on a healthy run is 3e-3 (planar) and
1.4e-4 (axisymmetric).

**A periodic event with a zero-initialised increment fires once.**
`qcc` classifies `t += tsnap` as an increment or a condition by evaluating it
twice at registration time, inside `_init_solver()`, which runs before
`main()`. With `tsnap` still `0` the value does not advance, so the clause is
classified as a *condition*, and that classification is never revisited --
the snapshot event fired at `t = 0` and never again. Both cases now give
`tsnap` a non-zero file-scope initialiser; the parameter-file value still
wins, because `init_event()` re-runs at `iter == 0`.

The capillary timestep constraint is **not** in this list. `tension.h`'s
`stability` event was instrumented against this solver stack and does apply
`dt_sigma = sqrt(rho_m*Delta_min^3/(pi*sigma))` correctly, using the actual
finest interface cell. The cases print `dt_sigma` at `MAXlevel` in their
header as the worst case a run can reach, so an achieved `dt` can be checked
against it, but they do not re-impose it.

## Planar case

`simulationCases/TaylorCulickPlanar.c` retracts a semi-infinite sheet of
**full** thickness `h0` that is symmetric about the midplane `y = 0`. The free
edge is closed by a semicircular cap of radius `h0/2` centred on the midplane
at `x = xtip0 + h0/2`, so the interface meets the midplane at `x = xtip0` and
retracts towards `+x`. The bottom boundary keeps Basilisk's default symmetry
condition, which is exactly the midplane mirror condition; left, right and top
are open.

It is the same code as the axisymmetric case apart from the geometry: no
`axi.h`, a planar initial condition and boundary conditions, no `AThTh`
conformation component (the solver already guards it with `#if AXI`), and no
`2*pi*y` weight in the kinetic-energy integral. The elastic and viscoelastic
capability is unchanged -- `G1`, `lambda1`, `G2`, `lambda2` behave exactly as
in the axisymmetric case, and `G1 = 0` gives the Newtonian limit by parameter
rather than by a separate stripped source file.

`Ldomain` must be large enough that the sheet stays effectively semi-infinite:
the edge travels at most `sqrt(2)*t`, so keep `Ldomain` well above
`xtip0 + sqrt(2)*tmax`.

### Non-dimensionalisation

Lengths are scaled with the **full** thickness `h0`, densities with the liquid
density, and stresses with `sigma/h0`. With `rho1 = sigma = h0 = 1`:

- `mu1` is the Ohnesorge number `Oh = mu/sqrt(rho*sigma*h0)`;
- the Taylor--Culick speed is `V_TC = sqrt(2*sigma/(rho*h0)) = sqrt(2)`;
- the capillary time is `sqrt(rho*h0^3/sigma) = 1`.

Savva & Bush (*JFM* **626**, 2009) put the **half**-thickness in their
Ohnesorge number,

```
Oh_SB = mu/sqrt(2*h0*rho*sigma) = Oh/sqrt(2),
```

so a case quoted at `Oh_SB` is run here with `mu1 = sqrt(2)*Oh_SB` -- the
factor of `sqrt(2)` is easy to lose. Their viscous time is
`tau_vis = mu*h0/(2*sigma) = mu1/2` and their reduced time is
`t* = t/tau_vis`. Both clocks are written to the tip file.

### Tip diagnostic

The planar case measures the retraction in-code and writes `c<CaseNo>-tip.dat`
with columns `t`, `tstar`, `x_tip`, `x_tip_vof`, `x_tip_global`:

- `x_tip` is the smallest `x` over the reconstructed VOF facets inside the
  midplane band `y < h0/10`, i.e. the interface position on the midplane;
- `x_tip_vof` is the independent estimate `integral (1-f) dx` along the bottom
  row of cells;
- `x_tip_global` is the same minimum taken over the whole interface.

All three are sub-cell accurate. `x_tip` and `x_tip_vof` agree while the
midplane is crossed exactly once, so a growing gap between them flags rim
pinch-off or an entrained bubble rather than a genuine tip motion.

Differentiate offline:

```bash
python3 postProcess/tip_to_csv.py simulationCases/c1000
```

which writes `tip_velocity.csv` with `t,tstar,x_tip,v_tip,v_over_VTC`. The
speed is a local first-order least-squares slope over `--window` samples,
which is far less noisy than a two-point difference at the sub-cell scale of
the VOF reconstruction.

## Post-processing from snapshots

`event writing_files` already `dump()`s a restartable snapshot every `tsnap`
into `intermediate/`. `postProcess/get_facets.c`, `get_tip.c` and
`get_fields.c` are small standalone Basilisk programs that `restore()` one of
those snapshots directly, rather than depending on anything the run wrote
about itself, so they can also re-derive the tip diagnostic and extract
fields the running case never logged. `get_tip.c`'s recomputed `x_tip` has
been checked bit-for-bit against the in-code `c<CaseNo>-tip.dat` log across a
full 81-snapshot run; a mismatch would flag a bug in one of the two
independent paths.

`get_fields.c` also computes `PHI`, the viscous dissipation rate per unit
volume `2*mu*E:E` (`E` the strain-rate tensor), the planar analogue of the
axisymmetric invariant in `comphy-lab/DropsAtLubis`'s `getData.c` and the
same reconstruction `Taylor-Culick-FEM/postProcess/animate_planar.py` uses
for its co-moving FEM cross-check. It is compiled once plain (`f`, `u.x`,
`u.y`, `KAPPA`, `PHI`) and once with `-DVISCOELASTIC=1` (adds
`A11 A12 A22 T11 T12 T22` from `log-conform-viscoelastic-scalar-2D.h`) --
compile the variant matching the case a snapshot came from; requesting a
field absent from that snapshot reads back as exactly zero rather than
erroring (`restore()`'s documented behaviour for a requested-but-undumped
field name).

Run the full pipeline over a case directory:

```bash
bash postProcess/run_postprocess.sh --case-dir <run-dir> [--viscoelastic] --video
```

This compiles the three programs with `qcc -disable-dimensions` (a plain
post-processing pass over restored data, not a physics solve, so Basilisk's
dimensional-consistency checker is not relevant here), then for every
`intermediate/snapshot-*` writes `<run-dir>/postprocess/{tip_snapshots.csv,
facets/, fields/, wide_fields/}`. Two field windows are extracted per
snapshot: a comoving one centred on the tip (`--window`/`--window-y`, in
units of `h0`; clamped to `[0, Ldomain]` since early in a run the tip sits
close to `x=0`), and a wide one fixed at `[0, --wide-fraction * Ldomain]`
for every snapshot (not tip-tracked -- it shows the domain's own span, not
a moving window on the rim).

`--video` renders PNG frames and an mp4 via `plot_fields.py`: two panels per
frame, both mirrored about the midplane `y=0` (the simulated half-domain is
symmetric, so mirroring shows the physical full-sheet cross-section) --
velocity magnitude (`Blues`, white at `|u|=0`) over the wide window, and
`PHI` (`hot_r`, `LogNorm`) over the comoving window, styled after
`Taylor-Culick-FEM/postProcess/animate_planar.py`. It needs matplotlib,
ffmpeg and a working system LaTeX (`text.usetex`) -- matplotlib, numpy and
ffmpeg are in the `elastic-tc-postprocess` conda env (`conda create -n
elastic-tc-postprocess --override-channels -c conda-forge python=3.11 numpy
matplotlib ffmpeg`), but LaTeX must come from the system install (`sudo apt
install texlive-latex-extra texlive-fonts-recommended dvipng cm-super`):
conda-forge's own `texlive-core` does not ship `latex.ltx` at all, so do not
`conda install` it into this env -- it would shadow the working
`/usr/bin/latex` on `PATH` with a package that cannot build `latex.fmt`.
Pass `--mathtext` to `plot_fields.py` to fall back to matplotlib's
no-subprocess renderer if this is ever parallelised (LaTeX spawns one
process per unique string and deadlocks under `multiprocessing`).

## Phase rheology

two-phaseVE.h maps the parameters below to the scalar solver:

- f = 1: rho1, mu1, G1, lambda1
- f = 0: rho2, mu2, G2, lambda2

Use a finite positive lambda for a viscoelastic liquid. Use lambda = 1e30
for the upstream purely elastic limit. Both phases can therefore be assigned
independently, including an elastic solid and a finite-relaxation liquid.
- `TOLelastic` is the phase-fraction cutoff used at mixed cells; it defaults to
  `1e-2` and should remain small compared with one. The runtime range is
  `0 <= TOLelastic < 0.5`, which keeps at least one phase active in every
  mixed cell.

## Parameter sweep

Preview the checked-in two-case sweep:

```bash
bash runParameterSweep.sh --dry-run
```

CASE_START..CASE_END must equal the Cartesian product generated by all
SWEEP_* entries. Each generated case is executed through runSimulation.sh,
so the single-case contract remains the only compilation path.

## Validation boundary

The local checks establish parser, runner, and compilation contracts. They do
not by themselves verify Taylor--Culick convergence or validate the physical
model against independent data. Those claims require separate refinement and
external-comparison evidence.
