# CHAMP CI tests

Every folder in this directory is one integration test: the CHAMP input
files plus a declarative manifest called **`test.json`** that says which
program to run, on how many MPI ranks, and which quantities must agree
with **reference values within statistical error bars**.

There are no shell scripts involved. A single cross-platform driver,
[`champ_test_runner.py`](champ_test_runner.py) (Python ≥ 3.6, standard
library only), stages each test into a scratch directory inside the build
tree, launches the MPI runs, extracts the results and applies the checks.
CMake discovers the manifests automatically and registers one `ctest`
entry per case.

## Running tests

```bash
cd build
ctest -N                     # list every test
ctest -j 4                   # run everything, 4 cores in parallel
ctest -L VMC                 # only tests labelled VMC
ctest -L DMC -j 4            # only DMC tests
ctest -LE "TREXIO|QMCKL"     # everything except TREXIO/QMCKL tests
ctest -R VMC-H2              # tests whose name matches a regex
ctest --print-labels         # show the label vocabulary
ctest --rerun-failed --output-on-failure
```

Each case advertises how many cores it needs (the `PROCESSORS` property),
so `ctest -j N` schedules tests without oversubscribing the machine.

Tests that need an optional feature (TREXIO, QMCkl) declare it in their
manifest; if the build was configured without that feature the test is
skipped at configure time and listed in the CMake output.

All run artifacts land in `build/tests/CI_test/scratch/<test>/<case>/` —
the source tree is never written to. To inspect a failure, look in that
directory: it contains the staged inputs, the CHAMP output file and the
`error` file of the failing run.

### Choosing the MPI launcher

The launcher is taken from CMake's FindMPI and can be overridden at
configure time, e.g. on a cluster or to allow oversubscription:

```bash
cmake -DCHAMP_TEST_MPIEXEC=/usr/bin/srun \
      -DCHAMP_TEST_NUMPROC_FLAG=-n \
      -DCHAMP_TEST_MPIEXEC_PREFLAGS="--oversubscribe" ...
```

## The pass criterion

QMC results are stochastic, so checks are statistical. A value check
passes when

```
|obtained − reference| ≤ nsigma · sqrt(error_ref² + error_obtained²)
```

where `error_ref` is the error bar stored in the manifest,
`error_obtained` is the error bar reported by the run itself, and
`nsigma` defaults to **2** (a deliberate ~5% false-failure rate per check
under ideal statistics; raise `nsigma` for noisy quantities rather than
inflating the reference error bar).

## Adding a new test

1. Create a folder with the CHAMP input files (input `.inp`, wave
   function, `pool/`, etc.). Everything **committed to git** in that
   folder is what gets staged at test time, so commit exactly the files
   the run needs.
2. Generate the manifest — no need to write JSON by hand. The wizard
   pre-fills everything it can detect (input file, vmc/dmc, TREXIO/QMCkl
   requirements, labels from the folder name) and writes a schema-valid
   `test.json` with placeholder reference values:

   ```bash
   tests/CI_test/new_test.py VMC-H2O-my-new-test
   ```

   Press Enter to accept a suggested default, or pass `--yes` (plus any
   flags, see `--help`) to skip the questions entirely. Editing the JSON
   by hand is only needed for the advanced features described below
   (multi-run pipelines, forces checks, file ops, ...).
3. Re-run `cmake` once so the new test is registered, then run it:

   ```bash
   cd build && cmake . && ctest -R my-new-test
   ```

   The check fails (the references are placeholders) but the output
   already shows the measured value.
4. Let the runner fill in the measured numbers for you:

   ```bash
   tests/CI_test/champ_test_runner.py suggest \
       --manifest tests/CI_test/My-New-Test/test.json \
       --case np2 \
       --scratch build/tests/CI_test/scratch/My-New-Test/np2
   ```

   Paste the printed `checks` block into the manifest. Ideally the
   reference comes from an independent long run (smaller error bar than
   the test run), but the suggested values are statistically consistent
   by construction.
5. Validate and re-run: `ctest -R my-new-test --output-on-failure`.

A note on reproducibility: the reference is a *distribution*, not a bit
pattern. The same test compiled with a different compiler or run with a
different rank count samples differently, which is why every reference
carries an error bar and why per-rank-count cases (`np1`, `np2`) have
separate references.

## Manifest reference

```jsonc
{
  "description": "Water VMC full optimization (TREXIO)",
  "labels": ["VMC", "TREXIO", "H2O", "optimization"],   // ctest -L selectors
  "requires": ["trexio"],          // optional: trexio and/or qmckl
  "enabled": true,                 // optional: false parks the test
  "timeout": 1500,                 // optional: seconds per case (default 1500)
  "cases": [
    {
      "name": "np2",               // unique; letters, digits, . _ -
      "comment": "free text",
      "runs": [
        {
          "program": "vmc",        // "vmc" or "dmc"
          "input": "vmc.inp",
          "output": "vmc_core_2.out",
          "nproc": 2,
          "error_output": "error", // optional (default "error")
          "before": [              // optional declarative file ops
            {"op": "remove", "paths": ["mc_configs"]}
          ],
          "after": [
            {"op": "concat", "sources": ["mc_configs_new*"], "dest": "mc_configs"},
            {"op": "remove", "paths": ["mc_configs_new*"]}
          ]
        }
      ],
      "checks": [ ... ]
    }
  ]
}
```

A *case* is one ctest entry. A case may chain several runs (e.g. VMC
generating walkers for DMC, or a checkpoint/restart pair); the checks are
evaluated after all runs finished. JSON does not allow comments — use the
`comment` fields (any check, case or run accepts one).

### Check types

**`value`** — compare a scalar printed as `<match> = <value> +- <error>`:

```json
{"type": "value", "file": "vmc_core_2.out", "match": "total E",
 "value": -17.228156, "error": 0.0038242}
```

- `match` is a whitespace-insensitive token sequence anchored at the start
  of the line (DMC example: `"total energy ( 100)"`).
- `occurrence` selects which matching line: `"last"` (default, i.e. the
  final optimization iteration), `"first"`, or an integer (`1` = first,
  `-1` = last, `-2` = second to last — used for multi-state outputs).
- `nsigma` (default 2.0) widens/narrows the acceptance window.
- `policy: "warn"` reports a mismatch without failing the test (used for
  statistically noisy methods).

**`difference`** — a derived quantity from two occurrences of the same
scalar, e.g. an excitation energy in eV:

```json
{"type": "difference", "file": "vmc_core_2.out", "match": "total E",
 "minuend": -1, "subtrahend": -2, "scale": 27.2114,
 "value": 4.5101235, "error": 0.3400853,
 "comment": "S1-S0 excitation energy in eV"}
```

The run-side error bar is the two extracted error bars combined in
quadrature and scaled.

**`forces`** — element-wise comparison of a force table
(`index fx fy fz σx σy σz` per atom) against a committed reference file,
each component within `nsigma` of the combined error:

```json
{"type": "forces", "file": "force_analytic", "reference": "force_reference"}
```

**`file_exists`** — assert the run produced a file (e.g. a checkpoint):

```json
{"type": "file_exists", "file": "restart_vmc.hdf5"}
```

**`values_equal`** — assert two runs printed the *identical* value (e.g.
the TREXIO text vs HDF5 backend comparison):

```json
{"type": "values_equal",
 "files": ["a_core_1.out", "b_core_1.out"], "match": "total E"}
```

### Labels

Labels drive test selection (`ctest -L`) locally and in the GitHub
workflows. Conventions: the method (`VMC`, `DMC`), the system (`H2O`,
`C4H6`, ...), the wave-function source (`TREXIO`, `QMCKL`), and the kind
of test (`energy`, `optimization`, `forces`, `restart`, ...). Keep using
existing labels where they fit — `ctest --print-labels` shows the
current vocabulary.

### Validating manifests

Manifests are strictly validated (unknown keys are rejected, so typos
cannot silently disable a check). Validation runs automatically during
`cmake`; to run it by hand:

```bash
tests/CI_test/champ_test_runner.py list tests/CI_test/*/test.json
```

The runner's own logic is covered by `ctest -R runner-selftest` (or
`champ_test_runner.py selftest`), which needs no CHAMP binaries.
