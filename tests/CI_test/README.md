# CHAMP CI tests

Each folder in this directory is one integration test: CHAMP input files
plus a declarative manifest, `test.json`, specifying which program to
run, on how many MPI ranks, and which quantities must agree with
reference values within statistical error bars.

[`champ_test_runner.py`](champ_test_runner.py) (Python ≥ 3.6, standard
library only) stages each test into a scratch directory, launches the
MPI runs and applies the checks. CMake discovers the manifests and
registers one ctest entry per case.

## Running a test directly

```bash
tests/CI_test/champ_test_runner.py VMC-H2          # run all H2 cases
tests/CI_test/champ_test_runner.py VMC-H2 --case energy-np2
tests/CI_test/champ_test_runner.py all             # run every test
tests/CI_test/champ_test_runner.py                 # list available tests
cd tests/CI_test/VMC-H2 && ../champ_test_runner.py # run from inside
```

`all` runs the cases sequentially; for parallel execution of the full
suite use `ctest -j N` (below).

Defaults: binaries `bin/vmc.mov1` / `bin/dmc.mov1` at the repository
root, MPI launcher `mpirun`/`mpiexec` from `PATH`, outputs in
`tests/CI_test/scratch/<test>/<case>/` (gitignored). Override with
`--vmc`, `--dmc`, `--mpiexec` or the environment variables `CHAMP_VMC`,
`CHAMP_DMC`, `CHAMP_MPIRUN`, e.g. to test a specific build:

```bash
tests/CI_test/champ_test_runner.py VMC-H2 --vmc /path/to/debug/vmc.mov1
```

## Running the suite with ctest

```bash
cd build
ctest -N                     # list every test
ctest -j 4                   # run everything, 4 cores in parallel
ctest -L quick               # fast smoke subset
ctest -L VMC                 # only tests labelled VMC
ctest -LE "TREXIO|QMCKL"     # everything except TREXIO/QMCKL tests
ctest -R VMC-H2              # tests whose name matches a regex
ctest --print-labels         # show the label vocabulary
ctest --rerun-failed --output-on-failure
```

ctest runs the binaries of the build directory it is invoked in: a
`-DCMAKE_BUILD_TYPE=Debug` build tests the Debug executables. To test a
different binary, configure with
`-DCHAMP_TEST_VMC_BINARY=/path/to/vmc.mov1` (and `CHAMP_TEST_DMC_BINARY`).
All builds place executables in `<source>/bin`, so build directories of
different build type overwrite each other there; build and test one mode
at a time, or use the override.

Each case declares its core count (`PROCESSORS` property), so `ctest -j N`
schedules tests without oversubscribing. Tests requiring an optional
feature (TREXIO, QMCkl) are skipped at configure time when the build
lacks it, and listed in the CMake output.

Run artifacts go to `build/tests/CI_test/scratch/<test>/<case>/`; the
source tree is never written to. To inspect a failure, look there: the
directory contains the staged inputs, the CHAMP output file and the
`error` file.

### MPI launcher

Taken from FindMPI; override at configure time, e.g. for `srun` or to
allow oversubscription:

```bash
cmake -DCHAMP_TEST_MPIEXEC=/usr/bin/srun \
      -DCHAMP_TEST_NUMPROC_FLAG=-n \
      -DCHAMP_TEST_MPIEXEC_PREFLAGS="--oversubscribe" ...
```

## Pass criterion

QMC results are stochastic, so checks are statistical. A value check
passes when

```
|obtained − reference| ≤ nsigma · sqrt(error_ref² + error_obtained²)
```

where `error_ref` is the error bar stored in the manifest,
`error_obtained` is the error bar reported by the run, and `nsigma` is
2 (≈5% false-failure rate per check under ideal statistics). The window
is fixed: checks known to fluctuate beyond it on some toolchains are
marked `policy: "warn"` so the deviation is reported without failing
the test.

The reference is a distribution, not a bit pattern: a different
compiler, different flags or a different rank count samples differently.
This is why every reference carries an error bar and why per-rank-count
cases (`np1`, `np2`) have separate references.

## Seeing the numbers of a run

Every case writes a JSON report of what it measured -- obtained value,
error bar, the reference it was compared against and the deviation in
sigma -- for the checks that passed as well as the ones that failed.
ctest writes them to `build/tests/CI_test/reports/`; `collect` turns the
whole run into one table:

```bash
ctest --test-dir build --output-on-failure --no-tests=error -j "$(nproc)"
tests/CI_test/champ_test_runner.py collect --reports build/tests/CI_test/reports
```

```
test/case               check           obtained                 reference                sigma  verdict
----------------------  --------------  -----------------------  -----------------------  -----  -------
VMC-H2/energy-np1       total E         -1.0046458 +- 0.0084583  -1.0046458 +- 0.0084583  0.00   PASS
VMC-H2/corsamp-np1      total E         -0.9781120 +- 0.0102988  -0.9923642 +- 0.0102988  0.98   PASS
```

Interactive runs collect the same reports with `--report-dir`:

```bash
tests/CI_test/champ_test_runner.py all --report-dir reports
tests/CI_test/champ_test_runner.py collect --reports reports
```

### From a GitHub Actions run

ctest prints a test's output only when it fails, so a green run -- or a
`policy: "warn"` check that drifted without failing -- would otherwise
show no numbers at all. Every workflow that runs the CI tests therefore
calls `.github/actions/champ-test-values` after ctest with
`if: always()`. It **reports only**; it never edits the manifests. The
table appears in the job summary, and the job uploads an artifact
`champ-reference-values-<toolchain>`:

| file | contents |
| ---- | -------- |
| `champ-obtained-values.txt` | the table above |
| `champ-suggested-references.txt` | the same numbers as paste-ready `checks` blocks, per case |
| `champ-reports/*.report.json` | the raw per-case reports |

### Runs that abort

`fatal_error` (`src/module/m_error.f90`) writes `Fatal error: <msg>` to the
output and error files and then calls `mpi_abort(..., 0, ...)` -- **with
abort code 0**. An aborted run therefore exits with status 0 and cannot be
told apart from a successful one by exit code alone, while the output file
it leaves behind is partial: the last `total E` is an intermediate
optimisation step, not a result.

The runner scans both files for that marker after every run. When it is
present the case fails with the abort message and the tails of both files,
and **no value from the run is compared** -- otherwise a dead calculation
is reported as a plausible-looking sigma deviation against the reference,
which sends you looking for a physics change that is not there.

## Refreshing the references

The references are the numbers a converged run produces; an algorithm
change (a different SR update, a new sampling scheme) moves them, and
they then have to be rebased **once**, deliberately, by hand -- not by
CI, which would otherwise silently re-baseline whatever it measured on
whichever runner ran last.

Do it from the reports of one chosen run. Those may come from a GitHub
Actions job, which is the case when the tests pass locally but not on the
runners:

```bash
gh run download <run-id> -n champ-reference-values-intel -D ci-values
tests/CI_test/champ_test_runner.py collect --reports ci-values/champ-reports \
    --update-manifests --dry-run       # what would change
tests/CI_test/champ_test_runner.py collect --reports ci-values/champ-reports \
    --update-manifests                 # write it
git diff tests/CI_test                 # review before committing
```

or from a local run:

```bash
tests/CI_test/champ_test_runner.py all --report-dir reports
tests/CI_test/champ_test_runner.py collect --reports reports --update-manifests
```

Only the `value` and `error` of each check are rewritten; comments,
policies, key order and formatting stay as they are, so the diff is
exactly the numbers that moved. A check whose `match` or `type` no longer
agrees with what the run measured is reported and left alone. Reports
carry the absolute manifest path of the machine that produced them and
fall back to the test folder name in this checkout (or in `--base DIR`),
so reports downloaded from a runner update the right files.

To transcribe values by hand instead, `--suggest` prints the paste-ready
`checks` blocks (this is what `champ-suggested-references.txt` holds).

Two things to weigh before committing a rebase: references are
toolchain-dependent samples, so adopting one runner's values can push
another toolchain outside the 2-sigma window -- check the other jobs'
tables. And a test-length run has a wide error bar; a reference from a
long, well-converged run is worth more. Checks that are genuinely noisy
belong under `policy: "warn"` rather than in a repeated rebase.

## Adding a new test

1. Create a folder with the CHAMP input files (`.inp`, wave function,
   `pool/`, ...). Only files **committed to git** are staged at test
   time, so commit exactly what the run needs.
2. Generate the manifest:

   ```bash
   tests/CI_test/new_test.py VMC-H2O-my-new-test
   ```

   Prompts are pre-filled from the folder contents (input file, program,
   TREXIO/QMCkl requirements, labels); `--yes` accepts all defaults.
   Hand-editing is only needed for the advanced features below
   (multi-run pipelines, forces checks, file ops).
3. Run it:

   ```bash
   tests/CI_test/champ_test_runner.py My-New-Test
   ```

   The check fails (the references are placeholders) and reports the
   measured value.
4. Fill in the measured reference values:

   ```bash
   tests/CI_test/champ_test_runner.py suggest \
       --manifest tests/CI_test/My-New-Test/test.json --case np2
   ```

   Paste the printed `checks` block into the manifest. Preferably the
   reference comes from an independent long run (smaller error bar than
   the test run); the suggested values are statistically consistent by
   construction.
5. Re-run step 3; ctest registers the test at the next `cmake` run.

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

A case is one ctest entry. A case may chain several runs (VMC generating
walkers for DMC, a checkpoint/restart pair); checks are evaluated after
all runs finish. JSON has no comments — use the `comment` fields (every
check, case and run accepts one).

### Check types

**`value`** — compare a scalar printed as `<match> = <value> +- <error>`:

```json
{"type": "value", "file": "vmc_core_2.out", "match": "total E",
 "value": -17.228156, "error": 0.0038242}
```

- `match`: whitespace-insensitive token sequence anchored at the start
  of the line (DMC example: `"total energy ( 100)"`).
- `occurrence`: which matching line; `"last"` (default, i.e. the final
  optimization iteration), `"first"`, or an integer (`1` = first,
  `-1` = last, `-2` = second to last — used for multi-state outputs).
- `policy: "warn"`: report a mismatch without failing the test
  (the per-test equivalent of continue-on-error).

**`difference`** — derived quantity from two occurrences of the same
scalar, e.g. an excitation energy in eV:

```json
{"type": "difference", "file": "vmc_core_2.out", "match": "total E",
 "minuend": -1, "subtrahend": -2, "scale": 27.2114,
 "value": 4.5101235, "error": 0.3400853}
```

The run-side error bar is the two extracted error bars combined in
quadrature and scaled.

**`forces`** — element-wise comparison of a force table
(`index fx fy fz σx σy σz` per atom) against a committed reference file,
each component within `nsigma` of the combined error:

```json
{"type": "forces", "file": "force_analytic", "reference": "force_reference"}
```

**`file_exists`** — assert the run produced a file:

```json
{"type": "file_exists", "file": "restart_vmc.hdf5"}
```

**`values_equal`** — assert two runs printed the identical value
(TREXIO text vs HDF5 backend comparison):

```json
{"type": "values_equal",
 "files": ["a_core_1.out", "b_core_1.out"], "match": "total E"}
```

### Labels

Labels drive test selection (`ctest -L`). Conventions: method (`VMC`,
`DMC`), system (`H2O`, `C4H6`, ...), wave-function source (`TREXIO`,
`QMCKL`), kind of test (`energy`, `optimization`, `forces`, `restart`,
...), and `quick` for sub-minute smoke tests (run per toolchain/build
mode by the matrix workflow). Reuse existing labels where they fit;
`ctest --print-labels` shows the vocabulary.

### Validating manifests

Manifests are strictly validated; unknown keys are rejected, so a typo
cannot silently disable a check. Validation runs during `cmake`, or by
hand:

```bash
tests/CI_test/champ_test_runner.py list tests/CI_test/*/test.json
```

The runner's own logic is covered by `ctest -R runner-selftest`
(equivalently `champ_test_runner.py selftest`); no CHAMP binaries
needed.

## Subcommands

| subcommand | purpose |
| ---------- | ------- |
| *(none)* | run test folders interactively (`champ_test_runner.py VMC-H2`) |
| `list` | validate manifests and enumerate cases (used by CMake) |
| `run` | run one case in a given scratch dir (used by ctest); `--report FILE` dumps the obtained values |
| `collect` | aggregate reports into a table; `--update-manifests` rebases the references on them (a manual step, never run by CI) |
| `suggest` | paste-ready `checks` block for one case from its scratch dir |
| `selftest` | the runner's own unit tests |
