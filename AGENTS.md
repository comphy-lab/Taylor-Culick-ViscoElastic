# Taylor-Culick-ViscoElastic repository guidance

This is a CoMPhy Basilisk project. Keep the canonical structure visible:

- simulationCases/ contains the active Basilisk entry point and generated case
  directories.
- src-local/ contains project-specific headers and the single runtime
  parameter API.
- Root runSimulation.sh and runParameterSweep.sh are the only supported entry
  points for case execution.

Use qcc resolved from PATH or BASILISK; never commit a machine-local compiler
path. The active constitutive solver is
src-local/log-conform-viscoelastic-scalar-2D.h, copied from the pinned
MultiRheoFlow source recorded in README.md. Do not reintroduce the retired
root-level _v0/_v1 headers.

Parameter files use key=value. C code reads them through src-local/params.h;
shell code reads and updates them through scripts/params.sh. Keep CaseNo
deterministic and preserve generated case directories as local run artifacts.

Classify evidence honestly: compile and runner checks are software tests;
refinement against an independent mathematical result is verification; an
external experiment or independently generated dataset is validation.
