#!/usr/bin/env python3
"""CHAMP CI test runner.

Usage:

    tests/CI_test/champ_test_runner.py VMC-H2            # one test folder
    tests/CI_test/champ_test_runner.py VMC-H2 --case np2 # one case
    tests/CI_test/champ_test_runner.py all               # every test
    tests/CI_test/champ_test_runner.py                   # inside a test
                                                         # folder: run it;
                                                         # elsewhere: list
                                                         # available tests

Defaults: binaries bin/vmc.mov1 and bin/dmc.mov1 at the repository root,
MPI launcher mpirun/mpiexec from PATH, outputs in
tests/CI_test/scratch/<test>/<case>/.  Overrides: --vmc/--dmc/--mpiexec
or the environment variables CHAMP_VMC, CHAMP_DMC, CHAMP_MPIRUN.

Each test folder carries a declarative ``test.json`` manifest describing
one or more *cases* (an MPI run -- or a pipeline of runs -- followed by
*checks* against reference values with statistical error bars).  ctest
uses the plumbing subcommands:

  list     enumerate and validate manifests at configure time; prints one
           tab-separated registration line per enabled case
  run      execute a single case inside a given scratch directory
           (used as the ctest COMMAND)
  suggest  after a run, print a ready-to-paste "checks" JSON snippet
           with the values measured in the scratch directory
  selftest exercise the parsing/comparison logic without CHAMP binaries

Pass criterion for stochastic values:

    |obtained - reference| <= nsigma * sqrt(err_ref^2 + err_obtained^2)

with nsigma = 2.  Both the reference error bar and the error bar
reported by the run enter the combined uncertainty; references are
toolchain-dependent samples of a distribution, not bit patterns.
"""

from __future__ import print_function

import argparse
import fnmatch
import glob as globmod
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

# ----------------------------------------------------------------------
# Manifest schema
# ----------------------------------------------------------------------

MANIFEST_NAME = "test.json"

KNOWN_REQUIRES = ("trexio", "qmckl")
PROGRAMS = ("vmc", "dmc")
DEFAULT_TIMEOUT = 1500          # seconds, mirrors the ctest default
DEFAULT_NSIGMA = 2.0
DEFAULT_ERROR_FILE = "error"

TOP_KEYS = {
    "description": str,
    "labels": list,
    "requires": list,
    "enabled": bool,
    "timeout": int,
    "cases": list,
}
CASE_KEYS = {
    "name": str,
    "comment": str,
    "runs": list,
    "checks": list,
}
RUN_KEYS = {
    "program": str,
    "input": str,
    "output": str,
    "nproc": int,
    "error_output": str,
    "before": list,
    "after": list,
}
CHECK_KEYS_COMMON = {
    "type": str,
    "comment": str,
    "policy": str,
}
CHECK_KEYS = {
    "value": {"file": str, "match": str, "value": float, "error": float,
              "occurrence": (str, int)},
    "difference": {"file": str, "match": str, "minuend": int,
                   "subtrahend": int, "scale": float, "value": float,
                   "error": float},
    "forces": {"file": str, "reference": str},
    "file_exists": {"file": str},
    "values_equal": {"files": list, "match": str, "occurrence": (str, int)},
}
OP_KEYS = {
    "remove": {"paths": list},
    "concat": {"sources": list, "dest": str},
}

# Files never copied into the scratch directory when git is unavailable
# (with git, exactly the tracked files are staged, which is authoritative).
FALLBACK_EXCLUDES = [
    "*.out", "*.log", "error", "error_*", "restart_*",
    "mc_configs_start", "mc_configs_new*", "temp_*",
    "det_optimal*", "*_optimal.*", "geo_optimal*", "force_analytic",
    "*.sh", "*.inactive", "*.debug", "*.update_this_test",
    ".DS_Store", "__pycache__",
]


class ManifestError(Exception):
    """Raised for any schema violation; message names the offending key."""


def _fail(where, msg):
    raise ManifestError("%s: %s" % (where, msg))


def _type_ok(value, typ):
    if isinstance(value, bool):
        return typ is bool
    if typ is float:
        return isinstance(value, (int, float))
    return isinstance(value, typ)


def _check_keys(obj, allowed, where):
    for key in obj:
        if key not in allowed:
            _fail(where, "unknown key '%s' (allowed: %s)"
                  % (key, ", ".join(sorted(allowed))))
    for key, typ in allowed.items():
        if key in obj:
            types = typ if isinstance(typ, tuple) else (typ,)
            if not any(_type_ok(obj[key], t) for t in types):
                names = "/".join(t.__name__ for t in types)
                _fail(where, "key '%s' must be of type %s" % (key, names))


def _require(obj, key, where):
    if key not in obj:
        _fail(where, "missing required key '%s'" % key)
    return obj[key]


def _safe_relpath(path, where):
    if not isinstance(path, str) or not path:
        _fail(where, "path must be a non-empty string")
    norm = path.replace("\\", "/")
    if norm.startswith("/") or ".." in norm.split("/") or ":" in norm:
        _fail(where, "path '%s' must be relative to the test folder" % path)
    return path


def _validate_ops(ops, where):
    for i, op in enumerate(ops):
        w = "%s.op[%d]" % (where, i)
        if not isinstance(op, dict):
            _fail(w, "must be an object")
        kind = _require(op, "op", w)
        if kind not in OP_KEYS:
            _fail(w, "unknown op '%s' (allowed: %s)"
                  % (kind, ", ".join(sorted(OP_KEYS))))
        allowed = dict(OP_KEYS[kind])
        allowed["op"] = str
        _check_keys(op, allowed, w)
        for key in OP_KEYS[kind]:
            _require(op, key, w)
        for key in ("paths", "sources"):
            if key in op:
                for p in op[key]:
                    _safe_relpath(p, w)
        if "dest" in op:
            _safe_relpath(op["dest"], w)


def _validate_check(check, where):
    if not isinstance(check, dict):
        _fail(where, "must be an object")
    ctype = _require(check, "type", where)
    if ctype not in CHECK_KEYS:
        _fail(where, "unknown check type '%s' (allowed: %s)"
              % (ctype, ", ".join(sorted(CHECK_KEYS))))
    allowed = dict(CHECK_KEYS_COMMON)
    allowed.update(CHECK_KEYS[ctype])
    _check_keys(check, allowed, where)
    if check.get("policy", "fail") not in ("fail", "warn"):
        _fail(where, "policy must be 'fail' or 'warn'")
    occurrence = check.get("occurrence", "last")
    if isinstance(occurrence, str):
        if occurrence not in ("first", "last"):
            _fail(where, "occurrence must be 'first', 'last' or a non-zero "
                         "integer (1 = first, -1 = last, -2 = second last...)")
    elif occurrence == 0:
        _fail(where, "occurrence 0 is meaningless; use 1 (first) or -1 (last)")
    if ctype == "value":
        for key in ("file", "match", "value", "error"):
            _require(check, key, where)
        _safe_relpath(check["file"], where)
        if float(check["error"]) < 0:
            _fail(where, "'error' (the reference error bar) must be >= 0")
    elif ctype == "difference":
        for key in ("file", "match", "minuend", "subtrahend", "value", "error"):
            _require(check, key, where)
        _safe_relpath(check["file"], where)
        if check["minuend"] == 0 or check["subtrahend"] == 0:
            _fail(where, "minuend/subtrahend are occurrence indices and "
                         "must be non-zero (1 = first, -1 = last)")
        if float(check["error"]) < 0:
            _fail(where, "'error' (the reference error bar) must be >= 0")
        if "scale" in check and float(check["scale"]) == 0:
            _fail(where, "'scale' must be non-zero")
    elif ctype == "forces":
        for key in ("file", "reference"):
            _safe_relpath(_require(check, key, where), where)
    elif ctype == "file_exists":
        _safe_relpath(_require(check, "file", where), where)
    elif ctype == "values_equal":
        files = _require(check, "files", where)
        _require(check, "match", where)
        if not isinstance(files, list) or len(files) != 2:
            _fail(where, "'files' must be a list of exactly two output files")
        for f in files:
            _safe_relpath(f, where)


def _validate_run(run, where):
    if not isinstance(run, dict):
        _fail(where, "must be an object")
    _check_keys(run, RUN_KEYS, where)
    program = _require(run, "program", where)
    if program not in PROGRAMS:
        _fail(where, "program must be one of: %s" % ", ".join(PROGRAMS))
    _safe_relpath(_require(run, "input", where), where)
    _safe_relpath(_require(run, "output", where), where)
    nproc = _require(run, "nproc", where)
    if not isinstance(nproc, int) or isinstance(nproc, bool) or nproc < 1:
        _fail(where, "nproc must be an integer >= 1")
    if "error_output" in run:
        _safe_relpath(run["error_output"], where)
    for key in ("before", "after"):
        if key in run:
            _validate_ops(run[key], "%s.%s" % (where, key))


def load_manifest(path):
    """Load and strictly validate a manifest; returns the parsed dict."""
    where = os.path.relpath(path)
    try:
        with open(path, "r") as fh:
            data = json.load(fh)
    except ValueError as exc:
        raise ManifestError("%s: invalid JSON: %s" % (where, exc))
    except OSError as exc:
        raise ManifestError("%s: %s" % (where, exc))
    if not isinstance(data, dict):
        _fail(where, "top level must be a JSON object")
    _check_keys(data, TOP_KEYS, where)

    labels = _require(data, "labels", where)
    if not labels or not all(isinstance(l, str) and l for l in labels):
        _fail(where, "'labels' must be a non-empty list of strings")
    for label in labels:
        if any(c in label for c in ";,\t \"'"):
            _fail(where, "label '%s' must not contain spaces, quotes, "
                         "commas, semicolons or tabs" % label)

    for req in data.get("requires", []):
        if req not in KNOWN_REQUIRES:
            _fail(where, "unknown requirement '%s' (known: %s)"
                  % (req, ", ".join(KNOWN_REQUIRES)))

    if data.get("timeout", DEFAULT_TIMEOUT) <= 0:
        _fail(where, "'timeout' must be a positive number of seconds")

    cases = _require(data, "cases", where)
    if not cases and data.get("enabled", True):
        _fail(where, "'cases' must not be empty for an enabled test")
    seen = set()
    for i, case in enumerate(cases):
        cwhere = "%s: cases[%d]" % (where, i)
        if not isinstance(case, dict):
            _fail(cwhere, "must be an object")
        _check_keys(case, CASE_KEYS, cwhere)
        name = _require(case, "name", cwhere)
        ok_chars = all(c.isalnum() or c in "._-" for c in name)
        if not name or not ok_chars:
            _fail(cwhere, "case name '%s' may only contain letters, digits, "
                          "'.', '_' and '-'" % name)
        if name in seen:
            _fail(cwhere, "duplicate case name '%s'" % name)
        seen.add(name)
        runs = _require(case, "runs", cwhere)
        if not runs:
            _fail(cwhere, "'runs' must not be empty")
        for j, run in enumerate(runs):
            _validate_run(run, "%s.runs[%d]" % (cwhere, j))
        checks = _require(case, "checks", cwhere)
        if not checks:
            _fail(cwhere, "'checks' must not be empty")
        for j, check in enumerate(checks):
            _validate_check(check, "%s.checks[%d]" % (cwhere, j))
    return data


def find_case(manifest, name):
    for case in manifest["cases"]:
        if case["name"] == name:
            return case
    raise ManifestError("case '%s' not found in manifest" % name)


# ----------------------------------------------------------------------
# Value extraction from CHAMP output
# ----------------------------------------------------------------------

class ExtractionError(Exception):
    pass


def _parse_match_line(raw, match, want):
    rest = raw.split()[len(want):]
    if not rest or rest[0] != "=":
        raise ExtractionError(
            "matched line lacks '= <value>' after '%s': %r" % (match, raw.rstrip()))
    try:
        value = float(rest[1])
    except (IndexError, ValueError):
        raise ExtractionError(
            "cannot parse value after '%s =' in line: %r" % (match, raw.rstrip()))
    if len(rest) < 4 or rest[2] != "+-":
        raise ExtractionError(
            "no '+- <error>' after the value in line: %r" % raw.rstrip())
    try:
        error = float(rest[3])
    except ValueError:
        raise ExtractionError("cannot parse error bar in line: %r" % raw.rstrip())
    return value, error


def extract_value(lines, match, occurrence="last"):
    """Return (value, error) from the occurrence-th line starting with the
    whitespace-delimited token sequence ``match`` and shaped like

        <match tokens> = <value> +- <error> ...

    ``occurrence`` is 'first', 'last', or a 1-based index (negative counts
    from the end: -1 = last, -2 = second last).  This reproduces the
    semantics of the historical tools/compare_value.py (line anchored at
    the first token, statistical error read after '+-').
    """
    want = match.split()
    if not want:
        raise ExtractionError("empty match pattern")
    matches = [raw for raw in lines if raw.split()[:len(want)] == want]
    if not matches:
        raise ExtractionError("no line starting with '%s' found" % match)
    index = {"first": 1, "last": -1}.get(occurrence, occurrence)
    pos = index - 1 if index > 0 else index
    try:
        raw = matches[pos]
    except IndexError:
        raise ExtractionError(
            "requested occurrence %s of '%s' but only %d line(s) matched"
            % (occurrence, match, len(matches)))
    return _parse_match_line(raw, match, want)


def extract_from_file(path, match, occurrence="last"):
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        raise ExtractionError("cannot read '%s': %s" % (path, exc))
    return extract_value(lines, match, occurrence)


def read_forces_table(path):
    """Parse a CHAMP forces file: one row per atom,
    ``index fx fy fz sigx sigy sigz``; returns list of (values, sigmas)."""
    rows = []
    try:
        with open(path, "r") as fh:
            for ln, raw in enumerate(fh, 1):
                if not raw.strip():
                    continue
                cols = raw.split()
                if len(cols) < 7:
                    raise ExtractionError(
                        "%s:%d: expected 7 columns (index, 3 forces, "
                        "3 error bars), got %d" % (path, ln, len(cols)))
                try:
                    nums = [float(c) for c in cols[1:7]]
                except ValueError:
                    raise ExtractionError("%s:%d: non-numeric entry" % (path, ln))
                rows.append((nums[0:3], nums[3:6]))
    except OSError as exc:
        raise ExtractionError("cannot read '%s': %s" % (path, exc))
    if not rows:
        raise ExtractionError("'%s' contains no force rows" % path)
    return rows


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------

class CheckResult(object):
    def __init__(self, title, passed, details, policy="fail"):
        self.title = title
        self.passed = passed
        self.details = details          # list of printable lines
        self.policy = policy

    @property
    def verdict(self):
        if self.passed:
            return "PASS"
        return "WARN" if self.policy == "warn" else "FAIL"


def _fmt(x):
    return "%.7f" % x


def check_value(scratch, check):
    policy = check.get("policy", "fail")
    nsigma = DEFAULT_NSIGMA
    occurrence = check.get("occurrence", "last")
    title = "value '%s' (%s) in %s" % (check["match"], occurrence, check["file"])
    ref_value = float(check["value"])
    ref_error = float(check["error"])
    try:
        value, error = extract_from_file(
            os.path.join(scratch, check["file"]), check["match"], occurrence)
    except ExtractionError as exc:
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy)
    combined = math.sqrt(ref_error * ref_error + error * error)
    allowed = nsigma * combined
    diff = abs(value - ref_value)
    passed = diff <= allowed
    sigmas = (diff / combined) if combined > 0 else float("inf") if diff else 0.0
    details = [
        "obtained  : %s +- %s" % (_fmt(value), _fmt(error)),
        "reference : %s +- %s" % (_fmt(ref_value), _fmt(ref_error)),
        "|diff| = %s  allowed = %s (%.1f sigma combined)  -> %.2f sigma"
        % (_fmt(diff), _fmt(allowed), nsigma, sigmas),
    ]
    return CheckResult(title, passed, details, policy)


def check_difference(scratch, check):
    """Difference of two occurrences of the same quantity (e.g. an
    excitation energy: last minus second-to-last 'total E', scaled to eV).
    The run-side error bar is the two extracted error bars added in
    quadrature, scaled by |scale|."""
    policy = check.get("policy", "fail")
    nsigma = DEFAULT_NSIGMA
    scale = float(check.get("scale", 1.0))
    title = "difference of '%s' (occurrence %d minus %d, scale %g) in %s" % (
        check["match"], check["minuend"], check["subtrahend"], scale, check["file"])
    ref_value = float(check["value"])
    ref_error = float(check["error"])
    path = os.path.join(scratch, check["file"])
    try:
        va, ea = extract_from_file(path, check["match"], check["minuend"])
        vb, eb = extract_from_file(path, check["match"], check["subtrahend"])
    except ExtractionError as exc:
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy)
    value = (va - vb) * scale
    error = math.sqrt(ea * ea + eb * eb) * abs(scale)
    combined = math.sqrt(ref_error * ref_error + error * error)
    allowed = nsigma * combined
    diff = abs(value - ref_value)
    passed = diff <= allowed
    sigmas = (diff / combined) if combined > 0 else float("inf") if diff else 0.0
    details = [
        "terms     : %s (occ %d)  -  %s (occ %d)"
        % (_fmt(va), check["minuend"], _fmt(vb), check["subtrahend"]),
        "obtained  : %s +- %s" % (_fmt(value), _fmt(error)),
        "reference : %s +- %s" % (_fmt(ref_value), _fmt(ref_error)),
        "|diff| = %s  allowed = %s (%.1f sigma combined)  -> %.2f sigma"
        % (_fmt(diff), _fmt(allowed), nsigma, sigmas),
    ]
    return CheckResult(title, passed, details, policy)


def check_forces(scratch, check):
    policy = check.get("policy", "fail")
    nsigma = DEFAULT_NSIGMA
    title = "forces %s vs reference %s" % (check["file"], check["reference"])
    try:
        obtained = read_forces_table(os.path.join(scratch, check["file"]))
        reference = read_forces_table(os.path.join(scratch, check["reference"]))
    except ExtractionError as exc:
        return CheckResult(title, False, [str(exc)], policy)
    if len(obtained) != len(reference):
        return CheckResult(title, False,
                           ["row count differs: obtained %d, reference %d"
                            % (len(obtained), len(reference))], policy)
    worst = (-1.0, None)            # (sigma distance, description)
    failures = 0
    for atom, ((vals, sigs), (rvals, rsigs)) in enumerate(zip(obtained, reference), 1):
        for k in range(3):
            diff = abs(vals[k] - rvals[k])
            combined = math.sqrt(sigs[k] ** 2 + rsigs[k] ** 2)
            allowed = nsigma * combined
            if diff > allowed:
                failures += 1
            dist = (diff / combined) if combined > 0 else float("inf") if diff else 0.0
            if dist > worst[0]:
                worst = (dist, "atom %d %s: |%s - %s| = %s, allowed %s"
                         % (atom, "xyz"[k], _fmt(vals[k]), _fmt(rvals[k]),
                            _fmt(diff), _fmt(allowed)))
    components = 3 * len(obtained)
    details = ["%d/%d components within %.1f sigma"
               % (components - failures, components, nsigma)]
    if worst[1]:
        details.append("largest deviation: %s (%.2f sigma)" % (worst[1], worst[0]))
    return CheckResult(title, failures == 0, details, policy)


def check_file_exists(scratch, check):
    policy = check.get("policy", "fail")
    path = os.path.join(scratch, check["file"])
    title = "file exists: %s" % check["file"]
    if os.path.isfile(path):
        return CheckResult(title, True,
                           ["found (%d bytes)" % os.path.getsize(path)], policy)
    return CheckResult(title, False, ["file was not produced"], policy)


def check_values_equal(scratch, check):
    policy = check.get("policy", "fail")
    occurrence = check.get("occurrence", "last")
    f1, f2 = check["files"]
    title = "values of '%s' equal in %s and %s" % (check["match"], f1, f2)
    try:
        v1, e1 = extract_from_file(os.path.join(scratch, f1), check["match"], occurrence)
        v2, e2 = extract_from_file(os.path.join(scratch, f2), check["match"], occurrence)
    except ExtractionError as exc:
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy)
    details = ["%s : %s +- %s" % (f1, _fmt(v1), _fmt(e1)),
               "%s : %s +- %s" % (f2, _fmt(v2), _fmt(e2))]
    return CheckResult(title, v1 == v2, details, policy)


CHECKERS = {
    "value": check_value,
    "difference": check_difference,
    "forces": check_forces,
    "file_exists": check_file_exists,
    "values_equal": check_values_equal,
}


# ----------------------------------------------------------------------
# Staging (source folder -> scratch directory)
# ----------------------------------------------------------------------

def _git_tracked_files(folder):
    """Return tracked paths (relative) or None if git/work-tree unavailable."""
    git = shutil.which("git")
    if git is None:
        return None
    try:
        proc = subprocess.run(
            [git, "-C", folder, "ls-files", "-z"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False)
    except OSError:
        return None
    if proc.returncode != 0:
        return None
    entries = [e for e in proc.stdout.decode("utf-8", "replace").split("\0") if e]
    return entries or None


def _copy_resolved(src, dst):
    """Copy a file or directory, resolving symlinks (several test inputs are
    symlinks into sibling test folders; a verbatim link would dangle)."""
    real = os.path.realpath(src)
    if not os.path.exists(real):
        raise RuntimeError("input '%s' is a broken link (target '%s' missing)"
                           % (src, real))
    if not os.path.isdir(os.path.dirname(dst)):
        os.makedirs(os.path.dirname(dst))
    if os.path.isdir(real):
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(real, dst)
    else:
        shutil.copy2(real, dst)


def _fallback_files(folder):
    """Without git: walk the folder, skipping known generated artifacts."""
    staged = []
    for root, dirs, files in os.walk(folder):
        rel_root = os.path.relpath(root, folder)
        dirs[:] = [d for d in dirs
                   if not any(fnmatch.fnmatch(d, p) for p in FALLBACK_EXCLUDES)]
        for name in files:
            rel = name if rel_root == "." else os.path.join(rel_root, name)
            if any(fnmatch.fnmatch(name, p) for p in FALLBACK_EXCLUDES):
                continue
            staged.append(rel)
        # also keep symlinked directories (os.walk does not descend by default)
        for name in list(dirs):
            full = os.path.join(root, name)
            if os.path.islink(full):
                rel = name if rel_root == "." else os.path.join(rel_root, name)
                staged.append(rel)
                dirs.remove(name)
    return staged


def stage(source, scratch):
    """(Re)create ``scratch`` with a pristine copy of the test inputs."""
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    entries = _git_tracked_files(source)
    origin = "git-tracked files"
    if entries is None:
        entries = _fallback_files(source)
        origin = "directory walk (git unavailable; generated files excluded)"
    count = 0
    for rel in entries:
        src = os.path.join(source, rel)
        if not os.path.exists(src) and not os.path.islink(src):
            # tracked but deleted locally; staging from a dirty tree
            print("[stage] warning: tracked file missing locally: %s" % rel)
            continue
        _copy_resolved(src, os.path.join(scratch, rel))
        count += 1
    print("[stage] %d entries staged from %s (%s)"
          % (count, source, origin))


# ----------------------------------------------------------------------
# Declarative file operations between runs
# ----------------------------------------------------------------------

def apply_ops(scratch, ops, phase):
    for op in ops:
        kind = op["op"]
        if kind == "remove":
            for pattern in op["paths"]:
                for path in globmod.glob(os.path.join(scratch, pattern)):
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                    print("[%s] removed %s" % (phase, os.path.relpath(path, scratch)))
        elif kind == "concat":
            sources = []
            for pattern in op["sources"]:
                sources.extend(sorted(globmod.glob(os.path.join(scratch, pattern))))
            if not sources:
                raise RuntimeError("concat: no files match %r" % (op["sources"],))
            dest = os.path.join(scratch, op["dest"])
            with open(dest, "wb") as out:
                for src in sources:
                    with open(src, "rb") as fh:
                        shutil.copyfileobj(fh, out)
            print("[%s] concatenated %d files into %s"
                  % (phase, len(sources), op["dest"]))


# ----------------------------------------------------------------------
# Defaults for interactive use
# ----------------------------------------------------------------------

RUNNER_DIR = os.path.dirname(os.path.abspath(__file__))      # tests/CI_test
REPO_ROOT = os.path.dirname(os.path.dirname(RUNNER_DIR))
DEFAULT_SCRATCH_ROOT = os.path.join(RUNNER_DIR, "scratch")   # gitignored


def default_binary(program):
    """bin/vmc.mov1 or bin/dmc.mov1 at the repository root;
    CHAMP_VMC/CHAMP_DMC override."""
    env = os.environ.get("CHAMP_" + program.upper())
    return env or os.path.join(REPO_ROOT, "bin", program + ".mov1")


def default_mpiexec():
    return (os.environ.get("CHAMP_MPIRUN")
            or shutil.which("mpirun") or shutil.which("mpiexec"))


def available_tests(base=None):
    """Sorted folder names under tests/CI_test that carry a manifest."""
    base = base or RUNNER_DIR
    found = []
    try:
        entries = sorted(os.listdir(base))
    except OSError:
        return found
    for name in entries:
        if os.path.isfile(os.path.join(base, name, MANIFEST_NAME)):
            found.append(name)
    return found


def resolve_test_dir(token, base=None):
    """Map a command-line token to a test folder: an explicit path, or a
    folder name under tests/CI_test."""
    base = base or RUNNER_DIR
    candidates = [token, os.path.join(base, token)]
    for cand in candidates:
        if os.path.isfile(os.path.join(cand, MANIFEST_NAME)):
            return os.path.abspath(cand)
    import difflib
    names = available_tests(base)
    close = difflib.get_close_matches(os.path.basename(token.rstrip("/\\")),
                                      names, n=3, cutoff=0.4)
    hint = ("; did you mean: " + ", ".join(close)) if close else \
           ("; available: " + ", ".join(names) if names else "")
    raise RuntimeError("no %s found for '%s'%s" % (MANIFEST_NAME, token, hint))


def default_scratch(folder_name, case_name):
    return os.path.join(DEFAULT_SCRATCH_ROOT, folder_name, case_name)


def expand_test_tokens(tokens, base=None):
    """Replace the special token 'all' with every test folder name.
    Returns (ordered unique tokens, set of names that came from 'all');
    disabled manifests are skipped only when pulled in via 'all'."""
    base = base or RUNNER_DIR
    out = []
    from_all = set()
    for token in tokens:
        if token == "all" and not os.path.isfile(
                os.path.join(base, "all", MANIFEST_NAME)):
            for name in available_tests(base):
                if name not in out:
                    out.append(name)
                    from_all.add(name)
        elif token not in out:
            out.append(token)
    return out, from_all


# ----------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------

def _tail(path, n=40):
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return []
    return [l.rstrip("\n") for l in lines[-n:]]


def build_command(run, args):
    binary = args.vmc if run["program"] == "vmc" else args.dmc
    if not binary:
        raise RuntimeError("no --%s binary was provided for a '%s' run"
                           % (run["program"], run["program"]))
    binary = os.path.abspath(binary)
    if not os.path.isfile(binary):
        raise RuntimeError("binary not found: %s" % binary)
    nproc = run["nproc"]
    base = [binary, "-i", run["input"], "-o", run["output"],
            "-e", run.get("error_output", DEFAULT_ERROR_FILE)]
    mpiexec = args.mpiexec
    if not mpiexec or mpiexec.endswith("-NOTFOUND"):
        if nproc == 1:
            return base      # MPI singleton run
        raise RuntimeError(
            "this case needs %d MPI ranks but no MPI launcher is configured "
            "(set CHAMP_TEST_MPIEXEC when configuring CMake)" % nproc)
    cmd = [mpiexec, args.numproc_flag, str(nproc)]
    cmd += args.mpiexec_preflags
    cmd += base
    return cmd


def run_case(args):
    manifest_path = os.path.abspath(args.manifest)
    source = os.path.dirname(manifest_path)
    manifest = load_manifest(manifest_path)
    case = find_case(manifest, args.case)
    scratch = os.path.abspath(args.scratch)

    print("=" * 72)
    print("CHAMP test  : %s / %s" % (os.path.basename(source), case["name"]))
    if manifest.get("description"):
        print("description : %s" % manifest["description"])
    if case.get("comment"):
        print("comment     : %s" % case["comment"])
    print("scratch dir : %s" % scratch)
    print("=" * 72)

    stage(source, scratch)

    for i, run in enumerate(case["runs"], 1):
        apply_ops(scratch, run.get("before", []), "before")
        input_path = os.path.join(scratch, run["input"])
        if not os.path.isfile(input_path):
            print("[run %d/%d] FAILED: input file '%s' not present in scratch dir"
                  % (i, len(case["runs"]), run["input"]))
            return 1
        cmd = build_command(run, args)
        print("[run %d/%d] %s" % (i, len(case["runs"]), " ".join(cmd)))
        sys.stdout.flush()
        # Sub-second runs can hit the still-held HDF5 file lock of the
        # previous restart dump (h5fcreate_f fails); single-writer test
        # runs do not need locking.  An explicit user setting wins.
        env = dict(os.environ)
        env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, cwd=scratch, env=env)
        except OSError as exc:
            print("[run %d/%d] FAILED to launch: %s" % (i, len(case["runs"]), exc))
            return 1
        elapsed = time.time() - t0
        print("[run %d/%d] exit code %d after %.1f s"
              % (i, len(case["runs"]), proc.returncode, elapsed))
        if proc.returncode != 0:
            print("---- tail of %s ----" % run["output"])
            for line in _tail(os.path.join(scratch, run["output"])):
                print("  " + line)
            err_name = run.get("error_output", DEFAULT_ERROR_FILE)
            print("---- tail of %s ----" % err_name)
            for line in _tail(os.path.join(scratch, err_name)):
                print("  " + line)
            print("RESULT: FAIL (program exited with %d)" % proc.returncode)
            return 1
        apply_ops(scratch, run.get("after", []), "after")

    print("-" * 72)
    results = []
    for i, check in enumerate(case["checks"], 1):
        result = CHECKERS[check["type"]](scratch, check)
        results.append(result)
        print("[check %d/%d] %s" % (i, len(case["checks"]), result.title))
        if check.get("comment"):
            print("    note: %s" % check["comment"])
        for line in result.details:
            print("    " + line)
        print("    %s" % result.verdict)

    hard_failures = [r for r in results if not r.passed and r.policy == "fail"]
    warnings = [r for r in results if not r.passed and r.policy == "warn"]
    print("-" * 72)
    print("RESULT: %s (%d/%d checks passed%s)"
          % ("FAIL" if hard_failures else "PASS",
             len(results) - len(hard_failures) - len(warnings), len(results),
             ", %d warn-only" % len(warnings) if warnings else ""))
    return 1 if hard_failures else 0


# ----------------------------------------------------------------------
# list mode (configure-time enumeration for CMake)
# ----------------------------------------------------------------------

def list_cases(args):
    status = 0
    for manifest_path in args.manifests:
        try:
            manifest = load_manifest(manifest_path)
        except ManifestError as exc:
            print("ERROR: %s" % exc, file=sys.stderr)
            status = 2
            continue
        folder = os.path.basename(os.path.dirname(os.path.abspath(manifest_path)))
        if not manifest.get("enabled", True):
            print("DISABLED\t%s" % folder, file=sys.stderr)
            continue
        requires = ",".join(manifest.get("requires", [])) or "-"
        labels = ",".join(manifest["labels"])
        timeout = manifest.get("timeout", DEFAULT_TIMEOUT)
        for case in manifest["cases"]:
            nproc = max(run["nproc"] for run in case["runs"])
            print("\t".join([os.path.abspath(manifest_path), folder,
                             case["name"], str(nproc), str(timeout),
                             requires, labels]))
    return status


# ----------------------------------------------------------------------
# suggest mode (fill in reference values from a previous run)
# ----------------------------------------------------------------------

def suggest(args):
    manifest = load_manifest(os.path.abspath(args.manifest))
    case = find_case(manifest, args.case)
    scratch = os.path.abspath(args.scratch)
    if not os.path.isdir(scratch):
        print("scratch directory '%s' does not exist; run the case first"
              % scratch, file=sys.stderr)
        return 2
    suggestion = []
    for check in case["checks"]:
        new = dict(check)
        try:
            path = os.path.join(scratch, check["file"]) if "file" in check else None
            if check["type"] == "value":
                new["value"], new["error"] = extract_from_file(
                    path, check["match"], check.get("occurrence", "last"))
            elif check["type"] == "difference":
                va, ea = extract_from_file(path, check["match"], check["minuend"])
                vb, eb = extract_from_file(path, check["match"], check["subtrahend"])
                scale = float(check.get("scale", 1.0))
                new["value"] = (va - vb) * scale
                new["error"] = math.sqrt(ea * ea + eb * eb) * abs(scale)
        except ExtractionError as exc:
            new["comment"] = "SUGGEST FAILED: %s" % exc
        suggestion.append(new)
    print(json.dumps({"checks": suggestion}, indent=2))
    return 0


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------

VMC_SAMPLE = """
 some header
total E =         -1.0156198 +-  0.0048495  0.68582  0.34470
total E =         -1.0046458 +-  0.0084583  0.68582  0.34470
"""

DMC_SAMPLE = """
total energy (   0) =   -26.2817058 +-  0.0128233  1.98470
total energy ( 100) =   -26.2823109 +-  0.0131892  1.99962
"""


def selftest(_args):
    import unittest

    class Extraction(unittest.TestCase):
        def test_last_and_first(self):
            lines = VMC_SAMPLE.splitlines()
            self.assertEqual(extract_value(lines, "total E", "last"),
                             (-1.0046458, 0.0084583))
            self.assertEqual(extract_value(lines, "total E", "first"),
                             (-1.0156198, 0.0048495))

        def test_dmc_keyword(self):
            lines = DMC_SAMPLE.splitlines()
            self.assertEqual(extract_value(lines, "total energy ( 100)"),
                             (-26.2823109, 0.0131892))

        def test_integer_occurrence(self):
            lines = VMC_SAMPLE.splitlines()
            self.assertEqual(extract_value(lines, "total E", 1),
                             extract_value(lines, "total E", "first"))
            self.assertEqual(extract_value(lines, "total E", -2),
                             (-1.0156198, 0.0048495))
            with self.assertRaises(ExtractionError):
                extract_value(lines, "total E", -3)

        def test_missing(self):
            with self.assertRaises(ExtractionError):
                extract_value(["nothing here"], "total E")

        def test_anchored_at_line_start(self):
            with self.assertRaises(ExtractionError):
                extract_value(["the total E = -1.0 +- 0.1"], "total E")

    class ValueCheck(unittest.TestCase):
        def setUp(self):
            self.dir = tempfile.mkdtemp()
            with open(os.path.join(self.dir, "o.out"), "w") as fh:
                fh.write(VMC_SAMPLE)

        def tearDown(self):
            shutil.rmtree(self.dir)

        def _check(self, ref, err, **kw):
            check = {"type": "value", "file": "o.out", "match": "total E",
                     "value": ref, "error": err}
            check.update(kw)
            return check_value(self.dir, check)

        def test_within_two_sigma(self):
            self.assertTrue(self._check(-1.0046458, 0.0084583).passed)
            # 2 sigma of combined error: sqrt(0.0084583^2+0.0084583^2)*2
            self.assertTrue(self._check(-1.0046458 + 0.023, 0.0084583).passed)
            self.assertFalse(self._check(-1.0046458 + 0.025, 0.0084583).passed)

        def test_warn_policy(self):
            res = self._check(-2.0, 0.0001, policy="warn")
            self.assertFalse(res.passed)
            self.assertEqual(res.verdict, "WARN")

        def test_zero_error_bars(self):
            res = self._check(-1.0046458, 0.0)
            self.assertTrue(res.passed)

        def test_difference_check(self):
            # gap = (last - second-last) * scale, errors in quadrature
            gap = (-1.0046458 - -1.0156198) * 27.2114
            gap_err = math.sqrt(0.0084583 ** 2 + 0.0048495 ** 2) * 27.2114
            check = {"type": "difference", "file": "o.out",
                     "match": "total E", "minuend": -1, "subtrahend": -2,
                     "scale": 27.2114, "value": gap, "error": gap_err}
            self.assertTrue(check_difference(self.dir, check).passed)
            check["value"] = gap + 10.0
            self.assertFalse(check_difference(self.dir, check).passed)

    class ForcesCheck(unittest.TestCase):
        def setUp(self):
            self.dir = tempfile.mkdtemp()
            self.rows = "1 0.001 0.002 -0.003 0.004 0.004 0.004\n" \
                        "2 -0.002 0.000 0.003 0.004 0.004 0.004\n"
            with open(os.path.join(self.dir, "ref"), "w") as fh:
                fh.write(self.rows)

        def tearDown(self):
            shutil.rmtree(self.dir)

        def test_identical_pass(self):
            with open(os.path.join(self.dir, "got"), "w") as fh:
                fh.write(self.rows)
            res = check_forces(self.dir, {"type": "forces", "file": "got",
                                          "reference": "ref"})
            self.assertTrue(res.passed)

        def test_outlier_fails(self):
            with open(os.path.join(self.dir, "got"), "w") as fh:
                fh.write("1 0.1 0.002 -0.003 0.004 0.004 0.004\n"
                         "2 -0.002 0.000 0.003 0.004 0.004 0.004\n")
            res = check_forces(self.dir, {"type": "forces", "file": "got",
                                          "reference": "ref"})
            self.assertFalse(res.passed)

    class Ops(unittest.TestCase):
        def setUp(self):
            self.dir = tempfile.mkdtemp()
            for n in ("mc_configs_new1", "mc_configs_new2"):
                with open(os.path.join(self.dir, n), "w") as fh:
                    fh.write(n + "\n")

        def tearDown(self):
            shutil.rmtree(self.dir)

        def test_concat_and_remove(self):
            apply_ops(self.dir, [
                {"op": "concat", "sources": ["mc_configs_new*"],
                 "dest": "mc_configs"},
                {"op": "remove", "paths": ["mc_configs_new*"]},
            ], "test")
            with open(os.path.join(self.dir, "mc_configs")) as fh:
                self.assertEqual(fh.read(), "mc_configs_new1\nmc_configs_new2\n")
            self.assertEqual(globmod.glob(os.path.join(self.dir, "mc_configs_new*")), [])

    class Schema(unittest.TestCase):
        def _load(self, payload):
            tmp = tempfile.mkdtemp()
            try:
                path = os.path.join(tmp, MANIFEST_NAME)
                with open(path, "w") as fh:
                    json.dump(payload, fh)
                return load_manifest(path)
            finally:
                shutil.rmtree(tmp)

        def _good(self):
            return {
                "description": "x", "labels": ["VMC"],
                "cases": [{"name": "np1",
                           "runs": [{"program": "vmc", "input": "a.inp",
                                     "output": "a.out", "nproc": 1}],
                           "checks": [{"type": "value", "file": "a.out",
                                       "match": "total E", "value": -1.0,
                                       "error": 0.01}]}],
            }

        def test_good_manifest(self):
            self.assertEqual(self._load(self._good())["labels"], ["VMC"])

        def test_unknown_key_rejected(self):
            bad = self._good()
            bad["cases"][0]["checks"][0]["refrence"] = -1.0
            with self.assertRaises(ManifestError):
                self._load(bad)

        def test_bad_requires_rejected(self):
            bad = self._good()
            bad["requires"] = ["sparkles"]
            with self.assertRaises(ManifestError):
                self._load(bad)

        def test_absolute_path_rejected(self):
            bad = self._good()
            bad["cases"][0]["runs"][0]["input"] = "/etc/passwd"
            with self.assertRaises(ManifestError):
                self._load(bad)

        def test_duplicate_case_rejected(self):
            bad = self._good()
            bad["cases"].append(json.loads(json.dumps(bad["cases"][0])))
            with self.assertRaises(ManifestError):
                self._load(bad)

    class InteractiveDefaults(unittest.TestCase):
        def setUp(self):
            self.base = tempfile.mkdtemp()
            os.makedirs(os.path.join(self.base, "VMC-Fake"))
            with open(os.path.join(self.base, "VMC-Fake", MANIFEST_NAME), "w") as fh:
                fh.write("{}")

        def tearDown(self):
            shutil.rmtree(self.base)

        def test_resolve_by_name_and_path(self):
            expect = os.path.join(self.base, "VMC-Fake")
            self.assertEqual(resolve_test_dir("VMC-Fake", self.base), expect)
            self.assertEqual(resolve_test_dir(expect, self.base), expect)

        def test_resolve_suggests_close_match(self):
            try:
                resolve_test_dir("VMC-Fak", self.base)
            except RuntimeError as exc:
                self.assertIn("VMC-Fake", str(exc))
            else:
                self.fail("expected RuntimeError")

        def test_all_token_expansion(self):
            tokens, from_all = expand_test_tokens(["all"], self.base)
            self.assertEqual(tokens, ["VMC-Fake"])
            self.assertEqual(from_all, {"VMC-Fake"})
            # explicit names are kept and deduplicated against 'all'
            tokens, from_all = expand_test_tokens(["VMC-Fake", "all"], self.base)
            self.assertEqual(tokens, ["VMC-Fake"])
            self.assertEqual(from_all, set())

        def test_available_and_binary_defaults(self):
            self.assertEqual(available_tests(self.base), ["VMC-Fake"])
            self.assertTrue(default_binary("vmc").endswith(
                os.path.join("bin", "vmc.mov1")))
            os.environ["CHAMP_VMC"] = "/custom/vmc"
            try:
                self.assertEqual(default_binary("vmc"), "/custom/vmc")
            finally:
                del os.environ["CHAMP_VMC"]

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    for cls in (Extraction, ValueCheck, ForcesCheck, Ops, Schema,
                InteractiveDefaults):
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    return 0 if runner.run(suite).wasSuccessful() else 1


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def run_tests(args):
    """Interactive mode: run all (or selected) cases of one or more test
    folders with default binaries, launcher and scratch locations."""
    if not args.tests:
        # bare invocation: run the test folder we are standing in, if any
        if os.path.isfile(os.path.join(os.getcwd(), MANIFEST_NAME)):
            args.tests = [os.getcwd()]
        else:
            names = available_tests()
            prog = os.path.basename(sys.argv[0])
            print("Available tests (tests/CI_test):\n")
            for name in names:
                print("    %s" % name)
            print("\nRun one with:    %s <name>" % prog)
            print("Run all with:    %s all" % prog)
            print("More options:    %s --help" % prog)
            return 0

    wanted_cases = []
    for spec in args.case or []:
        wanted_cases.extend(c for c in spec.split(",") if c)

    tokens, from_all = expand_test_tokens(args.tests)
    mpiexec = args.mpiexec or default_mpiexec()
    results = []                      # (folder, case, status-string, ok)
    for token in tokens:
        source = resolve_test_dir(token)
        folder = os.path.basename(source)
        manifest = load_manifest(os.path.join(source, MANIFEST_NAME))
        if not manifest.get("enabled", True):
            if token in from_all:
                print("note: skipping '%s' (disabled in its manifest)" % folder)
                continue
            print("note: '%s' is disabled in its manifest (parked test); "
                  "running anyway because it was named explicitly" % folder)
        for req in manifest.get("requires", []):
            print("note: %s requires a binary built with %s support"
                  % (folder, req.upper()))

        case_names = [c["name"] for c in manifest["cases"]]
        selected = case_names if not wanted_cases else \
            [c for c in case_names if c in wanted_cases]
        unknown = [c for c in wanted_cases
                   if c not in case_names and len(tokens) == 1]
        if unknown:
            raise RuntimeError("no case named %s in %s (available: %s)"
                               % (", ".join(unknown), folder,
                                  ", ".join(case_names)))

        for case_name in selected:
            run_args = argparse.Namespace(
                manifest=os.path.join(source, MANIFEST_NAME),
                case=case_name,
                scratch=default_scratch(folder, case_name),
                vmc=args.vmc or default_binary("vmc"),
                dmc=args.dmc or default_binary("dmc"),
                mpiexec=mpiexec,
                numproc_flag=args.numproc_flag,
                mpiexec_preflags=args.mpiexec_preflags.split())
            status = run_case(run_args)
            results.append((folder, case_name, status == 0))
            if status != 0:
                print("hint: outputs kept in %s" % run_args.scratch)
                print("hint: to refresh the reference values:  "
                      "%s suggest --manifest %s --case %s"
                      % (sys.argv[0], os.path.relpath(run_args.manifest),
                         case_name))
            print()

    if not results:
        if wanted_cases:
            raise RuntimeError("no case named %s in the selected tests"
                               % ", ".join(wanted_cases))
        raise RuntimeError("nothing to run")
    if len(results) > 1:
        print("=" * 72)
        print("Summary")
        for folder, case_name, ok in results:
            print("    %-50s %s" % (folder + "/" + case_name,
                                    "PASS" if ok else "FAIL"))
        npass = sum(1 for r in results if r[2])
        print("    %d/%d cases passed" % (npass, len(results)))
    return 0 if all(r[2] for r in results) else 1


USAGE = """\
Run CHAMP integration tests (see tests/CI_test/README.md).

    %(prog)s VMC-H2                    run every case of that test
    %(prog)s VMC-H2 DMC-C4H6-...      several tests in a row
    %(prog)s VMC-H2 --case energy-np2  a single case
    %(prog)s all                       run every integration test
    %(prog)s                           inside a test folder: run it;
                                       elsewhere: list available tests

Outputs go to tests/CI_test/scratch/<test>/<case>/ (wiped per run).

Defaults (each can be overridden by a flag or environment variable):
    vmc/dmc binaries   bin/vmc.mov1, bin/dmc.mov1 at the repository root
                       (--vmc/--dmc, $CHAMP_VMC/$CHAMP_DMC)
    MPI launcher       mpirun or mpiexec from PATH
                       (--mpiexec, $CHAMP_MPIRUN)

Subcommands (each has its own --help; ctest uses `run`):
    list / run / suggest / selftest
"""


def build_cli():
    parser = argparse.ArgumentParser(
        prog="champ_test_runner.py",
        usage=argparse.SUPPRESS,
        description=USAGE % {"prog": "champ_test_runner.py"},
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", metavar="")

    p_test = sub.add_parser("test", help="run test folders (default mode)")
    p_test.add_argument("tests", nargs="*", metavar="TEST",
                        help="test folder name (under tests/CI_test), a "
                             "path, or 'all' for every test")
    p_test.add_argument("--case", action="append", metavar="NAME",
                        help="run only this case (repeatable, or comma-"
                             "separated); default: all cases")
    p_test.add_argument("--vmc", default=None,
                        help="vmc.mov1 to test (default: repo bin/vmc.mov1 "
                             "or $CHAMP_VMC)")
    p_test.add_argument("--dmc", default=None,
                        help="dmc.mov1 to test (default: repo bin/dmc.mov1 "
                             "or $CHAMP_DMC)")
    p_test.add_argument("--mpiexec", default=None,
                        help="MPI launcher (default: mpirun/mpiexec from "
                             "PATH or $CHAMP_MPIRUN)")
    p_test.add_argument("--numproc-flag", default="-n",
                        help="launcher flag for the rank count (default -n)")
    p_test.add_argument("--mpiexec-preflags", default="",
                        help="extra launcher flags, e.g. --oversubscribe")

    p_list = sub.add_parser("list", help="validate manifests and list cases")
    p_list.add_argument("manifests", nargs="+")

    p_run = sub.add_parser("run", help="run one case (ctest plumbing)")
    p_run.add_argument("--manifest", required=True)
    p_run.add_argument("--case", required=True)
    p_run.add_argument("--scratch", required=True,
                       help="scratch directory (wiped and recreated)")
    p_run.add_argument("--vmc", default=None, help="path to vmc.mov1")
    p_run.add_argument("--dmc", default=None, help="path to dmc.mov1")
    p_run.add_argument("--mpiexec", default=None)
    p_run.add_argument("--numproc-flag", default="-n")
    p_run.add_argument("--mpiexec-preflags", default="",
                       help="extra launcher flags (space separated)")

    p_sug = sub.add_parser("suggest",
                           help="print checks filled with measured values")
    p_sug.add_argument("--manifest", required=True)
    p_sug.add_argument("--case", required=True)
    p_sug.add_argument("--scratch", default=None,
                       help="scratch dir of the run to read (default: the "
                            "interactive default for this test/case)")

    sub.add_parser("selftest", help="run the runner's own unit tests")
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Anything that is not a recognised subcommand means "run tests":
    # `champ_test_runner.py VMC-H2` works without ceremony.
    known = ("test", "list", "run", "suggest", "selftest")
    if not argv or (argv[0] not in known
                    and argv[0] not in ("-h", "--help")):
        argv.insert(0, "test")

    parser = build_cli()
    args = parser.parse_args(argv)
    if args.mode == "run":
        args.mpiexec_preflags = args.mpiexec_preflags.split()
    if args.mode == "suggest" and args.scratch is None:
        folder = os.path.basename(os.path.dirname(os.path.abspath(args.manifest)))
        args.scratch = default_scratch(folder, args.case)

    try:
        if args.mode == "test":
            return run_tests(args)
        if args.mode == "list":
            return list_cases(args)
        if args.mode == "run":
            return run_case(args)
        if args.mode == "suggest":
            return suggest(args)
        if args.mode == "selftest":
            return selftest(args)
    except ManifestError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    return 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\ninterrupted")
        sys.exit(130)
