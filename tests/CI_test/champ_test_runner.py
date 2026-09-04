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
           (used as the ctest COMMAND); ``--report FILE`` additionally
           dumps every obtained value as JSON
  collect  aggregate those JSON reports into one table of obtained vs
           reference values and, with ``--update-manifests``, write the
           obtained values back into the test.json files
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

# src/module/m_error.f90 writes this prefix to BOTH the output and the error
# unit and then calls mpi_abort(..., 0, ...) -- an aborted run therefore exits
# with status 0 and is indistinguishable from a successful one by exit code
# alone.  Without this marker the runner would compare whatever partial
# numbers the dead run had already printed.
FATAL_MARKER = "Fatal error:"

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
    def __init__(self, title, passed, details, policy="fail", record=None):
        self.title = title
        self.passed = passed
        self.details = details          # list of printable lines
        self.policy = policy
        # Structured record of what was measured: goes into the JSON
        # report written by `run --report` and read back by `collect`.
        self.record = record or {}

    @property
    def verdict(self):
        if self.passed:
            return "PASS"
        return "WARN" if self.policy == "warn" else "FAIL"


def _fmt(x):
    return "%.7f" % x


def _round(x):
    """Round to the 7 decimals CHAMP prints, so a suggested reference is
    written into a manifest exactly as the run reported it."""
    return round(float(x), 7)


def _sigma_field(sigmas):
    """JSON-safe sigma distance (nan/inf are not valid JSON numbers)."""
    if sigmas != sigmas or sigmas in (float("inf"), float("-inf")):
        return None
    return round(sigmas, 3)


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
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy,
                           {"match": check["match"], "file": check["file"],
                            "occurrence": occurrence,
                            "reference": _round(ref_value),
                            "reference_error": _round(ref_error),
                            "extraction_error": str(exc)})
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
    record = {
        "match": check["match"], "file": check["file"],
        "occurrence": occurrence,
        "obtained": _round(value), "obtained_error": _round(error),
        "reference": _round(ref_value), "reference_error": _round(ref_error),
        "sigma": _sigma_field(sigmas), "allowed": _round(allowed),
        "suggest": {"value": _round(value), "error": _round(error)},
    }
    return CheckResult(title, passed, details, policy, record)


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
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy,
                           {"match": check["match"], "file": check["file"],
                            "minuend": check["minuend"],
                            "subtrahend": check["subtrahend"],
                            "reference": _round(ref_value),
                            "reference_error": _round(ref_error),
                            "extraction_error": str(exc)})
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
    record = {
        "match": check["match"], "file": check["file"],
        "minuend": check["minuend"], "subtrahend": check["subtrahend"],
        "scale": scale,
        "terms": [_round(va), _round(vb)],
        "obtained": _round(value), "obtained_error": _round(error),
        "reference": _round(ref_value), "reference_error": _round(ref_error),
        "sigma": _sigma_field(sigmas), "allowed": _round(allowed),
        "suggest": {"value": _round(value), "error": _round(error)},
    }
    return CheckResult(title, passed, details, policy, record)


def check_forces(scratch, check):
    policy = check.get("policy", "fail")
    nsigma = DEFAULT_NSIGMA
    title = "forces %s vs reference %s" % (check["file"], check["reference"])
    try:
        obtained = read_forces_table(os.path.join(scratch, check["file"]))
        reference = read_forces_table(os.path.join(scratch, check["reference"]))
    except ExtractionError as exc:
        return CheckResult(title, False, [str(exc)], policy,
                           {"file": check["file"],
                            "extraction_error": str(exc)})
    if len(obtained) != len(reference):
        return CheckResult(title, False,
                           ["row count differs: obtained %d, reference %d"
                            % (len(obtained), len(reference))], policy,
                           {"file": check["file"],
                            "extraction_error": "row count differs"})
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
    record = {
        "file": check["file"], "reference_file": check["reference"],
        "components": components, "failed_components": failures,
        "sigma": _sigma_field(worst[0]) if worst[1] else None,
        "worst_component": worst[1],
        # The reference of a forces check is a committed file, not a
        # number in the manifest, so there is nothing to suggest here:
        # refresh it by copying the run's force table over it.
        "obtained_file": os.path.join(scratch, check["file"]),
    }
    return CheckResult(title, failures == 0, details, policy, record)


def check_file_exists(scratch, check):
    policy = check.get("policy", "fail")
    path = os.path.join(scratch, check["file"])
    title = "file exists: %s" % check["file"]
    if os.path.isfile(path):
        size = os.path.getsize(path)
        return CheckResult(title, True, ["found (%d bytes)" % size], policy,
                           {"file": check["file"], "exists": True,
                            "size_bytes": size})
    return CheckResult(title, False, ["file was not produced"], policy,
                       {"file": check["file"], "exists": False})


def check_values_equal(scratch, check):
    policy = check.get("policy", "fail")
    occurrence = check.get("occurrence", "last")
    f1, f2 = check["files"]
    title = "values of '%s' equal in %s and %s" % (check["match"], f1, f2)
    try:
        v1, e1 = extract_from_file(os.path.join(scratch, f1), check["match"], occurrence)
        v2, e2 = extract_from_file(os.path.join(scratch, f2), check["match"], occurrence)
    except ExtractionError as exc:
        return CheckResult(title, False, ["extraction failed: %s" % exc], policy,
                           {"match": check["match"], "file": f1,
                            "occurrence": occurrence,
                            "extraction_error": str(exc)})
    details = ["%s : %s +- %s" % (f1, _fmt(v1), _fmt(e1)),
               "%s : %s +- %s" % (f2, _fmt(v2), _fmt(e2))]
    # Both files must agree exactly; there is no manifest-side reference
    # value to refresh, so the second file plays the role of "reference".
    record = {
        "match": check["match"], "file": f1, "occurrence": occurrence,
        "obtained": _round(v1), "obtained_error": _round(e1),
        "reference": _round(v2), "reference_error": _round(e2),
        "reference_file": f2,
        "sigma": 0.0 if v1 == v2 else None,
    }
    return CheckResult(title, v1 == v2, details, policy, record)


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
# Reports (obtained values, machine-readable and as a table)
# ----------------------------------------------------------------------

REPORT_SUFFIX = ".report.json"


def report_filename(folder, case_name):
    """One file per ctest case; ctest -j writes them concurrently."""
    safe = "%s__%s" % (folder, case_name)
    for bad in "/\\ ":
        safe = safe.replace(bad, "_")
    return safe + REPORT_SUFFIX


def write_report(path, report):
    if not path:
        return
    directory = os.path.dirname(os.path.abspath(path))
    if directory and not os.path.isdir(directory):
        try:
            os.makedirs(directory)
        except OSError:
            pass
    try:
        with open(path, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=False)
            fh.write("\n")
    except OSError as exc:
        print("warning: could not write report '%s': %s" % (path, exc))
        return
    print("report written: %s" % path)


def load_reports(directory):
    """Read every *.report.json in `directory`, sorted by test/case."""
    paths = sorted(globmod.glob(os.path.join(directory, "*" + REPORT_SUFFIX)))
    reports = []
    for path in paths:
        try:
            with open(path) as fh:
                reports.append(json.load(fh))
        except (OSError, ValueError) as exc:
            print("warning: skipping unreadable report '%s': %s" % (path, exc),
                  file=sys.stderr)
    reports.sort(key=lambda r: (r.get("test", ""), r.get("case", "")))
    return reports


def _pair(value, error):
    if value is None:
        return "-"
    if error is None:
        return _fmt(value)
    return "%s +- %s" % (_fmt(value), _fmt(error))


def _check_label(entry):
    name = entry.get("match") or entry.get("file") or entry.get("type", "?")
    occurrence = entry.get("occurrence", "last")
    if occurrence != "last":
        name = "%s [%s]" % (name, occurrence)
    if entry.get("type") == "difference":
        name = "%s (occ %s - %s)" % (name, entry.get("minuend"),
                                     entry.get("subtrahend"))
    return name


def report_rows(reports):
    """Flatten reports into table rows: case, check, obtained, reference,
    sigma, verdict."""
    rows = []
    for rep in reports:
        label = "%s/%s" % (rep.get("test", "?"), rep.get("case", "?"))
        if not rep.get("checks"):
            rows.append([label, rep.get("error", "no checks evaluated"),
                         "-", "-", "-", rep.get("result", "ERROR")])
            continue
        for entry in rep["checks"]:
            note = _check_label(entry)
            if entry.get("extraction_error"):
                note = "%s (%s)" % (note, entry["extraction_error"])
            sigma = entry.get("sigma")
            rows.append([
                label, note,
                _pair(entry.get("obtained"), entry.get("obtained_error")),
                _pair(entry.get("reference"), entry.get("reference_error")),
                "%.2f" % sigma if sigma is not None else "-",
                entry.get("verdict", "?"),
            ])
    return rows


def format_table(rows, headers):
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    # the last column is not padded, so lines do not end in blanks
    def line(cells):
        out = [cells[i].ljust(widths[i]) for i in range(len(cells) - 1)]
        out.append(cells[-1])
        return "  ".join(out).rstrip()
    text = [line(headers), line(["-" * w for w in widths])]
    text.extend(line(row) for row in rows)
    return "\n".join(text)


REPORT_HEADERS = ["test/case", "check", "obtained", "reference",
                  "sigma", "verdict"]


def render_reports(reports, title="OBTAINED VALUES"):
    rows = report_rows(reports)
    out = ["=" * 72, title, "=" * 72]
    if rows:
        out.append(format_table(rows, REPORT_HEADERS))
    else:
        out.append("(no checks were evaluated)")
    verdicts = [row[-1] for row in rows]
    out.append("")
    out.append("%d checks: %d PASS, %d WARN, %d FAIL/ERROR"
               % (len(verdicts), verdicts.count("PASS"), verdicts.count("WARN"),
                  sum(1 for v in verdicts if v not in ("PASS", "WARN"))))
    deviations = []
    for rep in reports:
        for entry in rep.get("checks", []):
            if entry.get("sigma") is not None:
                deviations.append((entry["sigma"],
                                   "%s/%s  %s" % (rep.get("test"), rep.get("case"),
                                                  _check_label(entry))))
    deviations.sort(reverse=True)
    if deviations and deviations[0][0] > 0:
        out.append("")
        out.append("largest deviations from the stored references:")
        for sigma, where in deviations[:10]:
            out.append("    %6.2f sigma   %s" % (sigma, where))
    return "\n".join(out)


def manifest_label(manifest_path):
    """Readable manifest path: relative to the cwd when that stays inside
    the tree (interactive use), else the repository-style path (ctest runs
    from the build directory, whose relpath is a stack of '..')."""
    rel = os.path.relpath(manifest_path)
    if not rel.startswith(os.pardir):
        return rel
    folder = os.path.basename(os.path.dirname(os.path.abspath(manifest_path)))
    return "/".join(["tests", "CI_test", folder, MANIFEST_NAME])


def suggested_checks(report):
    """The case's checks with the obtained numbers substituted, ready to
    paste into the manifest."""
    out = []
    for entry in report.get("checks", []):
        suggest = entry.get("suggest")
        if not suggest:
            continue
        out.append(dict(entry.get("check", {}), **suggest))
    return out


def update_manifests(reports, base=None, dry_run=False):
    """Write the obtained values back into the test.json manifests.

    Only `value`/`error` of checks whose type and match still agree with
    what the run measured are touched; everything else in the manifest
    (comments, policies, key order, formatting) is preserved."""
    base = base or os.path.dirname(os.path.abspath(__file__))
    by_manifest = {}
    for rep in reports:
        path = rep.get("manifest") or ""
        if not os.path.isfile(path):
            # reports produced on another machine (a downloaded CI
            # artifact): fall back to the test folder name
            path = os.path.join(base, rep.get("test", ""), MANIFEST_NAME)
        by_manifest.setdefault(os.path.abspath(path), []).append(rep)

    changed, skipped = [], []
    for path, reps in sorted(by_manifest.items()):
        if not os.path.isfile(path):
            skipped.append("%s: manifest not found" % path)
            continue
        with open(path) as fh:
            data = json.load(fh)
        cases = dict((c.get("name"), c) for c in data.get("cases", []))
        edits = 0
        for rep in reps:
            case = cases.get(rep.get("case"))
            if case is None:
                skipped.append("%s: no case '%s'" % (path, rep.get("case")))
                continue
            checks = case.get("checks", [])
            for entry in rep.get("checks", []):
                suggest = entry.get("suggest")
                index = entry.get("index", 0) - 1
                if not suggest or not 0 <= index < len(checks):
                    continue
                check = checks[index]
                if (check.get("type") != entry.get("type")
                        or check.get("match") != entry.get("match")):
                    skipped.append("%s/%s check %d: manifest changed since "
                                   "the run, not updated"
                                   % (rep.get("test"), rep.get("case"), index + 1))
                    continue
                if (check.get("value") == suggest["value"]
                        and check.get("error") == suggest["error"]):
                    continue
                check["value"] = suggest["value"]
                check["error"] = suggest["error"]
                edits += 1
        if edits and not dry_run:
            with open(path, "w") as fh:
                json.dump(data, fh, indent=2)
                fh.write("\n")
        if edits:
            changed.append((path, edits))
    return changed, skipped


# ----------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------

def scan_fatal(path):
    """First 'Fatal error:' line in `path`, or None.

    CHAMP aborts with MPI code 0 (see m_error.f90), so a run that died
    mid-calculation still reports exit status 0.  The output it leaves
    behind is partial: the last 'total E' is some intermediate
    optimisation step, not a result, and comparing it against a reference
    yields a meaningless sigma deviation instead of an error."""
    try:
        with open(path, "r", errors="replace") as fh:
            for raw in fh:
                if FATAL_MARKER in raw:
                    return raw.strip()
    except OSError:
        pass
    return None


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

    report = {
        "test": os.path.basename(source),
        "case": case["name"],
        "manifest": manifest_path,
        "scratch": scratch,
        "result": "ERROR",
        "runs": [],
        "checks": [],
    }
    report_path = getattr(args, "report", None)

    stage(source, scratch)

    for i, run in enumerate(case["runs"], 1):
        apply_ops(scratch, run.get("before", []), "before")
        input_path = os.path.join(scratch, run["input"])
        if not os.path.isfile(input_path):
            print("[run %d/%d] FAILED: input file '%s' not present in scratch dir"
                  % (i, len(case["runs"]), run["input"]))
            report["error"] = ("input file '%s' not present in scratch dir"
                               % run["input"])
            write_report(report_path, report)
            return 1
        try:
            cmd = build_command(run, args)
        except RuntimeError as exc:
            print("[run %d/%d] FAILED: %s" % (i, len(case["runs"]), exc))
            report["error"] = str(exc)
            write_report(report_path, report)
            return 1
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
            report["error"] = "failed to launch: %s" % exc
            write_report(report_path, report)
            return 1
        elapsed = time.time() - t0
        err_name = run.get("error_output", DEFAULT_ERROR_FILE)
        report["runs"].append({
            "program": run["program"], "input": run["input"],
            "output": run["output"], "nproc": run["nproc"],
            "returncode": proc.returncode, "seconds": round(elapsed, 1)})
        print("[run %d/%d] exit code %d after %.1f s"
              % (i, len(case["runs"]), proc.returncode, elapsed))
        # An aborted run exits 0, so the marker -- not the status -- is what
        # tells us the numbers in the output file are not a result.
        fatal = (scan_fatal(os.path.join(scratch, run["output"]))
                 or scan_fatal(os.path.join(scratch, err_name)))
        if proc.returncode == 0 and fatal:
            print("[run %d/%d] the program aborted despite exit code 0:"
                  % (i, len(case["runs"])))
            print("    %s" % fatal)
            print("    its output is partial, so no value from it is compared")
            for name in (run["output"], err_name):
                print("---- tail of %s ----" % name)
                for line in _tail(os.path.join(scratch, name), 20):
                    print("  " + line)
            report["runs"][-1]["fatal"] = fatal
            report["error"] = fatal
            print("RESULT: FAIL (%s)" % fatal)
            write_report(report_path, report)
            return 1
        if proc.returncode != 0:
            print("---- tail of %s ----" % run["output"])
            for line in _tail(os.path.join(scratch, run["output"])):
                print("  " + line)
            print("---- tail of %s ----" % err_name)
            for line in _tail(os.path.join(scratch, err_name)):
                print("  " + line)
            print("RESULT: FAIL (program exited with %d)" % proc.returncode)
            report["error"] = "program exited with %d" % proc.returncode
            write_report(report_path, report)
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
        entry = {"index": i, "type": check["type"], "title": result.title,
                 "policy": result.policy, "verdict": result.verdict}
        if check.get("comment"):
            entry["comment"] = check["comment"]
        entry.update(result.record)
        # the manifest check itself, so `collect` can rebuild a paste-ready
        # block without re-reading the manifest
        entry["check"] = check
        report["checks"].append(entry)

    hard_failures = [r for r in results if not r.passed and r.policy == "fail"]
    warnings = [r for r in results if not r.passed and r.policy == "warn"]
    print("-" * 72)
    print("RESULT: %s (%d/%d checks passed%s)"
          % ("FAIL" if hard_failures else "PASS",
             len(results) - len(hard_failures) - len(warnings), len(results),
             ", %d warn-only" % len(warnings) if warnings else ""))

    report["result"] = "FAIL" if hard_failures else "PASS"
    report["warnings"] = len(warnings)
    # Always echo the numbers: with `ctest --output-on-failure` this block
    # is what a failing job shows, and `collect` reproduces it for the
    # cases that passed.
    print(render_reports([report],
                         title="OBTAINED VALUES  %s / %s"
                               % (report["test"], report["case"])))
    # A case whose every check passed needs no new reference, so only the
    # table is printed for it; `collect --suggest` prints the blocks of a
    # whole run when references are being rebased deliberately.
    snippet = suggested_checks(report) if (hard_failures or warnings) else []
    if snippet:
        print("")
        print("checks block with the values measured by this run "
              "(paste into %s):" % manifest_label(manifest_path))
        print(json.dumps({"checks": snippet}, indent=2))
    write_report(report_path, report)
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
# collect mode (aggregate the reports of a whole ctest run)
# ----------------------------------------------------------------------

def collect(args):
    """Turn the per-case reports of a ctest run into one table of
    obtained vs reference values -- for every case, not only the failing
    ones -- and optionally write the obtained values into the manifests."""
    if not os.path.isdir(args.reports):
        print("no report directory '%s': no test ran, or the build predates "
              "--report" % args.reports)
        return 0
    reports = load_reports(args.reports)
    if not reports:
        print("no *%s files in '%s': no test ran" % (REPORT_SUFFIX, args.reports))
        return 0

    if args.format == "json":
        print(json.dumps(reports, indent=2))
    elif not args.no_table:
        print(render_reports(reports))
        print("")
        print("%d case(s): %s" % (len(reports), ", ".join(
            "%d %s" % (sum(1 for r in reports if r.get("result") == v), v)
            for v in ("PASS", "FAIL", "ERROR")
            if any(r.get("result") == v for r in reports))))
        print("refresh the stored references with:  %s collect --reports %s "
              "--update-manifests" % (os.path.basename(sys.argv[0]), args.reports))

    if args.suggest:
        print("")
        print("=" * 72)
        print("CHECKS BLOCKS WITH THE MEASURED VALUES")
        print("=" * 72)
        for rep in reports:
            snippet = suggested_checks(rep)
            if not snippet:
                continue
            print("")
            print("# %s / %s  ->  %s"
                  % (rep.get("test"), rep.get("case"),
                     manifest_label(rep.get("manifest")
                                    or os.path.join(rep.get("test", ""),
                                                    MANIFEST_NAME))))
            print(json.dumps({"checks": snippet}, indent=2))

    if args.update_manifests:
        changed, skipped = update_manifests(reports, base=args.base,
                                            dry_run=args.dry_run)
        print("")
        print("=" * 72)
        print("REFERENCE VALUES %s"
              % ("THAT WOULD BE WRITTEN (dry run)" if args.dry_run
                 else "WRITTEN INTO THE MANIFESTS"))
        print("=" * 72)
        for path, edits in changed:
            print("    %-60s %d value(s)" % (os.path.relpath(path, args.base
                                                             or os.getcwd()), edits))
        if not changed:
            print("    (none: the stored references already match this run)")
        for note in skipped:
            print("    skipped: %s" % note)
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

    class FatalAbort(unittest.TestCase):
        """CHAMP aborts with mpi_abort(...,0,...), so exit status alone
        cannot distinguish a dead run from a finished one."""

        def setUp(self):
            self.dir = tempfile.mkdtemp()

        def tearDown(self):
            shutil.rmtree(self.dir)

        def _write(self, name, text):
            path = os.path.join(self.dir, name)
            with open(path, "w") as fh:
                fh.write(text)
            return path

        def test_marker_found_in_partial_output(self):
            path = self._write("o.out", VMC_SAMPLE +
                               "MATINV: u(k,k)=0 with k=    10\n"
                               "Fatal error: MATINV: info ne 0 in dgetrf\n")
            self.assertEqual(scan_fatal(path),
                             "Fatal error: MATINV: info ne 0 in dgetrf")

        def test_clean_output_and_warnings_are_not_fatal(self):
            self.assertIsNone(scan_fatal(self._write("clean.out", VMC_SAMPLE)))
            # a successful run's error file carries warnings, not failures
            self.assertIsNone(scan_fatal(self._write(
                "error", "Warning:: No information about orbital symmetries\n"
                         "INPUT: multideterminant bloc MISSING\n")))

        def test_missing_file_is_not_fatal(self):
            self.assertIsNone(scan_fatal(os.path.join(self.dir, "absent")))

        def test_partial_output_still_parses_a_value(self):
            # the reason the marker matters: extraction happily returns the
            # last intermediate energy from a dead run
            path = self._write("o.out", VMC_SAMPLE +
                               "Fatal error: MATINV: info ne 0 in dgetrf\n")
            value, _ = extract_from_file(path, "total E")
            self.assertEqual(value, -1.0046458)
            self.assertIsNotNone(scan_fatal(path))

    class Reporting(unittest.TestCase):
        """The report -> table -> manifest path used to refresh references."""

        def setUp(self):
            self.dir = tempfile.mkdtemp()
            with open(os.path.join(self.dir, "o.out"), "w") as fh:
                fh.write(VMC_SAMPLE)
            self.manifest = {
                "description": "x", "labels": ["VMC"],
                "cases": [{"name": "np1",
                           "runs": [{"program": "vmc", "input": "a.inp",
                                     "output": "o.out", "nproc": 1}],
                           "checks": [{"type": "value", "file": "o.out",
                                       "match": "total E", "value": -2.0,
                                       "error": 0.01, "comment": "keep me"}]}],
            }
            folder = os.path.join(self.dir, "VMC-Fake")
            os.makedirs(folder)
            self.path = os.path.join(folder, MANIFEST_NAME)
            with open(self.path, "w") as fh:
                json.dump(self.manifest, fh, indent=2)

        def tearDown(self):
            shutil.rmtree(self.dir)

        def _report(self):
            check = self.manifest["cases"][0]["checks"][0]
            result = check_value(self.dir, check)
            entry = {"index": 1, "type": check["type"], "title": result.title,
                     "policy": result.policy, "verdict": result.verdict}
            entry.update(result.record)
            entry["check"] = check
            return {"test": "VMC-Fake", "case": "np1", "manifest": self.path,
                    "result": "FAIL", "runs": [], "checks": [entry]}

        def test_record_carries_the_obtained_value(self):
            entry = self._report()["checks"][0]
            self.assertEqual(entry["obtained"], -1.0046458)
            self.assertEqual(entry["obtained_error"], 0.0084583)
            self.assertEqual(entry["suggest"],
                             {"value": -1.0046458, "error": 0.0084583})

        def test_table_shows_obtained_and_reference(self):
            text = render_reports([self._report()])
            self.assertIn("total E", text)
            self.assertIn("-1.0046458", text)      # obtained
            self.assertIn("-2.0000000", text)      # reference
            self.assertIn("FAIL", text)

        def test_report_round_trip(self):
            path = os.path.join(self.dir, report_filename("VMC-Fake", "np1"))
            write_report(path, self._report())
            back = load_reports(self.dir)
            self.assertEqual(len(back), 1)
            self.assertEqual(back[0]["case"], "np1")
            self.assertEqual(back[0]["checks"][0]["obtained"], -1.0046458)

        def test_update_manifest_replaces_only_the_numbers(self):
            changed, skipped = update_manifests([self._report()])
            self.assertEqual([edits for _, edits in changed], [1])
            self.assertEqual(skipped, [])
            with open(self.path) as fh:
                check = json.load(fh)["cases"][0]["checks"][0]
            self.assertEqual(check["value"], -1.0046458)
            self.assertEqual(check["error"], 0.0084583)
            self.assertEqual(check["comment"], "keep me")
            self.assertEqual(check["type"], "value")
            # already up to date: nothing more to write
            self.assertEqual(update_manifests([self._report()])[0], [])

        def test_update_is_skipped_when_the_manifest_moved_on(self):
            report = self._report()
            report["checks"][0]["match"] = "some other quantity"
            changed, skipped = update_manifests([report])
            self.assertEqual(changed, [])
            self.assertTrue(skipped)

        def test_dry_run_does_not_write(self):
            changed, _ = update_manifests([self._report()], dry_run=True)
            self.assertEqual([edits for _, edits in changed], [1])
            with open(self.path) as fh:
                self.assertEqual(json.load(fh)["cases"][0]["checks"][0]["value"],
                                 -2.0)

        def test_suggested_checks_keep_the_manifest_keys(self):
            snippet = suggested_checks(self._report())
            self.assertEqual(snippet[0]["comment"], "keep me")
            self.assertEqual(snippet[0]["value"], -1.0046458)
            self.assertEqual(snippet[0]["match"], "total E")

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
    for cls in (Extraction, ValueCheck, ForcesCheck, Ops, FatalAbort,
                Reporting, Schema,
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
                mpiexec_preflags=args.mpiexec_preflags.split(),
                report=os.path.join(args.report_dir,
                                    report_filename(folder, case_name))
                       if args.report_dir else None)
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
    if args.report_dir:
        print("    reports in %s (aggregate them with: %s collect --reports %s)"
              % (args.report_dir, os.path.basename(sys.argv[0]), args.report_dir))
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
    list / run / collect / suggest / selftest

Collecting the numbers of a whole run (e.g. to refresh references after
an algorithm change):

    %(prog)s all --report-dir reports
    %(prog)s collect --reports reports --update-manifests
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
    p_test.add_argument("--report-dir", default=None, metavar="DIR",
                        help="also write one JSON report of obtained values "
                             "per case into DIR (read back by `collect`)")

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
    p_run.add_argument("--report", default=None, metavar="FILE",
                       help="write a JSON report of every obtained value "
                            "(read back by `collect`)")

    p_sug = sub.add_parser("suggest",
                           help="print checks filled with measured values")
    p_sug.add_argument("--manifest", required=True)
    p_sug.add_argument("--case", required=True)
    p_sug.add_argument("--scratch", default=None,
                       help="scratch dir of the run to read (default: the "
                            "interactive default for this test/case)")

    p_col = sub.add_parser(
        "collect",
        help="aggregate run reports: obtained vs reference values")
    p_col.add_argument("--reports", required=True, metavar="DIR",
                       help="directory holding the *%s files" % REPORT_SUFFIX)
    p_col.add_argument("--format", choices=("table", "json"), default="table",
                       help="table (default) or the raw merged JSON")
    p_col.add_argument("--no-table", action="store_true",
                       help="suppress the table (e.g. when only the "
                            "paste-ready blocks are wanted)")
    p_col.add_argument("--suggest", action="store_true",
                       help="also print a paste-ready 'checks' block per case")
    p_col.add_argument("--update-manifests", action="store_true",
                       help="write the obtained values into the test.json "
                            "files as the new reference values")
    p_col.add_argument("--dry-run", action="store_true",
                       help="with --update-manifests: report what would "
                            "change without writing")
    p_col.add_argument("--base", default=None, metavar="DIR",
                       help="tests/CI_test directory to update (default: the "
                            "one this script lives in); used when the reports "
                            "come from another machine")

    sub.add_parser("selftest", help="run the runner's own unit tests")
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Anything that is not a recognised subcommand means "run tests":
    # `champ_test_runner.py VMC-H2` works without ceremony.
    known = ("test", "list", "run", "collect", "suggest", "selftest")
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
        if args.mode == "collect":
            return collect(args)
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
