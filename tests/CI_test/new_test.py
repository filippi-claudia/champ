#!/usr/bin/env python3
"""Generate a test.json manifest for a CHAMP CI test.

Interactive wizard (just run it) or fully scriptable via flags:

    tests/CI_test/new_test.py VMC-H2O-my-test
    tests/CI_test/new_test.py VMC-H2O-my-test --yes \
        --program vmc --input vmc.inp --nproc 1 2 --labels VMC H2O

The generated manifest contains placeholder reference values (0.0); run
the test once and let the runner measure them for you:

    cd build && cmake . && ctest -R VMC-H2O-my-test
    tests/CI_test/champ_test_runner.py suggest --manifest <...>/test.json \
        --case np2 --scratch build/tests/CI_test/scratch/VMC-H2O-my-test/np2

then paste the printed "checks" block into test.json.

Standard library only, Python >= 3.6.  See README.md for the full schema
(multi-run pipelines, forces/difference/file_exists checks, file ops).
"""

from __future__ import print_function

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import champ_test_runner as runner          # noqa: E402  (schema validation)

DMC_DEFAULT_MATCH = "total energy ( 100)"
VMC_DEFAULT_MATCH = "total E"
PLACEHOLDER_NOTE = ("PLACEHOLDER reference - run the test once, then fill "
                    "in measured values with champ_test_runner.py suggest")


# ----------------------------------------------------------------------
# small prompt helpers (skipped entirely with --yes)
# ----------------------------------------------------------------------

def ask(prompt, default="", assume_yes=False):
    if assume_yes:
        return default
    suffix = " [%s]" % default if default != "" else ""
    try:
        value = input("  %-28s%s : " % (prompt, suffix)).strip()
    except EOFError:                      # piped/non-interactive input ran out
        print()
        return default
    return value or default


def ask_list(prompt, default, assume_yes=False):
    raw = ask(prompt, " ".join(str(d) for d in default), assume_yes)
    if isinstance(raw, list):
        return raw
    return [tok for tok in re.split(r"[,\s]+", raw.strip()) if tok]


def ask_choice(prompt, choices, default, assume_yes=False):
    while True:
        value = ask("%s (%s)" % (prompt, "/".join(choices)), default, assume_yes)
        if value in choices:
            return value
        print("    please answer one of: %s" % ", ".join(choices))
        if assume_yes:                    # avoid an infinite loop in scripts
            raise SystemExit("invalid value %r for %s" % (value, prompt))


# ----------------------------------------------------------------------
# context sniffing: derive good defaults from the test folder
# ----------------------------------------------------------------------

def find_inputs(folder):
    if not os.path.isdir(folder):
        return []
    return sorted(f for f in os.listdir(folder) if f.endswith(".inp"))


def guess_program(folder_name, input_name):
    text = (folder_name + " " + (input_name or "")).lower()
    return "dmc" if "dmc" in text else "vmc"


def guess_labels(folder_name, program, requires):
    """Folder names like VMC-TREXIO-H2O-DFT-optall already carry good tags."""
    labels = []
    for token in re.split(r"[-_+]+", folder_name):
        if not token or token.isdigit():
            continue
        if len(token) < 2 or token.lower() in ("dets", "csfs", "csf", "test"):
            continue
        if token not in labels:
            labels.append(token)
    for tag in [program.upper()] + [r.upper() for r in requires]:
        if tag not in labels:
            labels.insert(0, tag)
    return labels[:8]


def guess_requires(folder, folder_name, input_name):
    requires = []
    blob = folder_name.lower()
    input_path = os.path.join(folder, input_name) if input_name else None
    if input_path and os.path.isfile(input_path):
        try:
            with open(input_path, errors="replace") as fh:
                blob += " " + fh.read().lower()
        except OSError:
            pass
    if "trexio" in blob:
        requires.append("trexio")
    if "qmckl" in blob:
        requires.append("qmckl")
    return requires


# ----------------------------------------------------------------------
# manifest assembly
# ----------------------------------------------------------------------

def build_manifest(description, labels, requires, timeout, program,
                   input_file, output_stem, nprocs, match):
    cases = []
    for np in nprocs:
        output = "%s_core_%d.out" % (output_stem, np)
        cases.append({
            "name": "np%d" % np,
            "runs": [{
                "program": program,
                "input": input_file,
                "output": output,
                "nproc": np,
            }],
            "checks": [{
                "type": "value",
                "file": output,
                "match": match,
                "value": 0.0,
                "error": 0.0,
                "comment": PLACEHOLDER_NOTE,
            }],
        })
    manifest = {"description": description, "labels": labels}
    if requires:
        manifest["requires"] = requires
    if timeout != runner.DEFAULT_TIMEOUT:
        manifest["timeout"] = timeout
    manifest["cases"] = cases
    return manifest


def next_steps(folder, manifest_path, first_case):
    rel = os.path.relpath
    name = os.path.basename(folder)
    return """
Wrote %s (schema-valid).

Next steps:
  1. commit the input files in %s
     (only git-tracked files are staged when the test runs)
  2. run it (it FAILS now: the reference values are placeholders):
       tests/CI_test/champ_test_runner.py %s
  3. fill in the measured reference values and error bars:
       tests/CI_test/champ_test_runner.py suggest \\
           --manifest %s --case %s
     and paste the printed "checks" block into test.json
  4. re-run until happy:
       tests/CI_test/champ_test_runner.py %s
     (ctest registers the test automatically at the next cmake run)

For multi-run pipelines, forces/difference/file_exists/values_equal
checks, warn-only policies or file ops, see tests/CI_test/README.md.
""" % (rel(manifest_path), rel(folder), name,
       rel(manifest_path), first_case, name)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a test.json manifest for a CHAMP CI test "
                    "(interactive unless --yes is given)")
    parser.add_argument("folder",
                        help="test folder (name or path); created if missing")
    parser.add_argument("--description", default=None)
    parser.add_argument("--program", choices=runner.PROGRAMS, default=None)
    parser.add_argument("--input", default=None, help="CHAMP input file (.inp)")
    parser.add_argument("--output-stem", default=None,
                        help="output files become <stem>_core_<N>.out")
    parser.add_argument("--nproc", type=int, nargs="+", default=None,
                        help="MPI rank counts, one ctest case each (default: 1 2)")
    parser.add_argument("--labels", nargs="+", default=None)
    parser.add_argument("--requires", nargs="*", choices=runner.KNOWN_REQUIRES,
                        default=None)
    parser.add_argument("--match", default=None,
                        help="output line to check (default: '%s', dmc: '%s')"
                             % (VMC_DEFAULT_MATCH, DMC_DEFAULT_MATCH))
    parser.add_argument("--timeout", type=int, default=None,
                        help="seconds per case (default %d)" % runner.DEFAULT_TIMEOUT)
    parser.add_argument("--yes", action="store_true",
                        help="accept all defaults, no prompts")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing test.json")
    args = parser.parse_args(argv)

    folder = args.folder.rstrip("/\\")
    if os.path.dirname(folder) == "" and not os.path.isdir(folder):
        folder = os.path.join(HERE, folder)       # bare name -> tests/CI_test/
    folder = os.path.abspath(folder)
    name = os.path.basename(folder)
    manifest_path = os.path.join(folder, runner.MANIFEST_NAME)

    if os.path.exists(manifest_path) and not args.force:
        raise SystemExit("%s already exists (use --force to overwrite)"
                         % os.path.relpath(manifest_path))

    print("CHAMP CI test generator -- folder: %s" % os.path.relpath(folder))
    if not os.path.isdir(folder):
        print("  note: folder does not exist yet, it will be created")

    inputs = find_inputs(folder)
    if inputs and args.input is None:
        print("  input files found: %s" % ", ".join(inputs))

    default_input = args.input or (inputs[0] if inputs else "vmc.inp")
    input_file = ask("input file", default_input, args.yes) if args.input is None \
        else args.input
    if input_file not in inputs and not args.yes:
        print("    note: '%s' is not in the folder yet -- remember to add it"
              % input_file)

    default_prog = args.program or guess_program(name, input_file)
    program = args.program or ask_choice("program", list(runner.PROGRAMS),
                                         default_prog, args.yes)

    description = args.description or ask(
        "description", name.replace("-", " ").replace("_", " "), args.yes)

    default_req = guess_requires(folder, name, input_file)
    requires = args.requires if args.requires is not None else ask_list(
        "requires (trexio/qmckl)", default_req, args.yes)
    for req in requires:
        if req not in runner.KNOWN_REQUIRES:
            raise SystemExit("unknown requirement '%s' (known: %s)"
                             % (req, ", ".join(runner.KNOWN_REQUIRES)))

    default_labels = args.labels or guess_labels(name, program, requires)
    labels = args.labels or ask_list("labels (ctest -L selectors)",
                                     default_labels, args.yes)

    nprocs = args.nproc or [int(n) for n in
                            ask_list("MPI rank counts", [1, 2], args.yes)]

    default_stem = args.output_stem or os.path.splitext(input_file)[0]
    output_stem = args.output_stem or ask("output stem", default_stem, args.yes)

    default_match = args.match or (DMC_DEFAULT_MATCH if program == "dmc"
                                   else VMC_DEFAULT_MATCH)
    match = args.match or ask("quantity to check", default_match, args.yes)

    timeout = args.timeout or int(ask("timeout (seconds)",
                                      runner.DEFAULT_TIMEOUT, args.yes))

    manifest = build_manifest(description, labels, requires, timeout,
                              program, input_file, output_stem, nprocs, match)

    if not os.path.isdir(folder):
        os.makedirs(folder)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    # The generator must never produce an invalid manifest: validate with
    # the same code cmake/ctest will use, and fail loudly if it does not.
    try:
        runner.load_manifest(manifest_path)
    except runner.ManifestError as exc:
        os.remove(manifest_path)
        raise SystemExit("internal error, generated manifest is invalid: %s" % exc)

    print(next_steps(folder, manifest_path, "np%d" % nprocs[0]))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\naborted, nothing written")
        sys.exit(130)
