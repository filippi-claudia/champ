#!/usr/bin/python3

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict, deque
from urllib import error, request


DEFAULT_ECP_BASE_URL = (
    "https://raw.githubusercontent.com/filippi-claudia/champ/main/pool/BFD/ECP_gamess"
)
DEFAULT_GITHUB_TREE_API_URL = (
    "https://api.github.com/repos/filippi-claudia/champ/git/trees/main?recursive=1"
)

_GITHUB_TREE_CACHE = None


def normalize_element(label):
    if isinstance(label, bytes):
        label = label.decode("utf-8", errors="ignore")
    text = str(label).strip()
    if not text:
        raise ValueError(f"Invalid empty nucleus label: {label!r}")

    match = re.match(r"([A-Za-z]+)", text)
    if not match:
        raise ValueError(f"Could not parse element from nucleus label: {label!r}")

    symbol = match.group(1)
    return symbol[0].upper() + symbol[1:].lower()


def parse_bfd_ecp_text(ecp_text, element):
    lines = [line.strip() for line in ecp_text.splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"ECP file for {element} is empty")

    if "GEN" in lines[0].upper():
        return parse_bfd_ecp_gamess_text(lines, element)

    return parse_bfd_ecp_champ_text(lines, element)


def parse_bfd_ecp_champ_text(lines, element):
    if len(lines) < 3:
        raise ValueError(f"ECP file for {element} is too short")

    try:
        n_components = int(float(lines[1].split()[0]))
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Could not parse number of components for {element}") from exc

    local_l = n_components - 1
    pos = 2
    channels = {}

    for component_index in range(n_components):
        if pos >= len(lines):
            raise ValueError(f"Unexpected end of ECP file for {element}")

        try:
            n_terms = int(float(lines[pos].split()[0]))
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Could not parse number of terms for {element}") from exc
        pos += 1

        channel_l = local_l if component_index == 0 else (component_index - 1)
        terms = []

        for _ in range(n_terms):
            if pos >= len(lines):
                raise ValueError(f"Unexpected end of ECP terms for {element}")
            fields = lines[pos].split()
            pos += 1
            if len(fields) < 3:
                raise ValueError(f"Malformed ECP term in file for {element}: {fields}")

            coefficient = float(fields[0])
            power_file = int(round(float(fields[1])))
            exponent = float(fields[2])
            power_trexio = power_file - 2
            terms.append((power_trexio, coefficient, exponent))

        channels[channel_l] = terms

    return channels


def parse_bfd_ecp_gamess_text(lines, element):
    if len(lines) < 2:
        raise ValueError(f"GAMESS ECP file for {element} is too short")

    header_fields = lines[0].split()
    lmax = None
    if len(header_fields) >= 2:
        try:
            lmax = int(float(header_fields[-1]))
        except ValueError:
            lmax = None

    pos = 1
    blocks = []
    while pos < len(lines):
        try:
            n_terms = int(float(lines[pos].split()[0]))
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Could not parse number of terms in GAMESS ECP for {element} at line: {lines[pos]!r}"
            ) from exc
        pos += 1

        terms = []
        for _ in range(n_terms):
            if pos >= len(lines):
                raise ValueError(f"Unexpected end of GAMESS ECP terms for {element}")
            fields = lines[pos].split()
            pos += 1
            if len(fields) < 3:
                raise ValueError(f"Malformed GAMESS ECP term in file for {element}: {fields}")

            coefficient = float(fields[0])
            power_file = int(round(float(fields[1])))
            exponent = float(fields[2])
            power_trexio = power_file - 2
            terms.append((power_trexio, coefficient, exponent))

        blocks.append(terms)

    if not blocks:
        raise ValueError(f"No ECP terms found in GAMESS file for {element}")

    if lmax is None:
        lmax = len(blocks) - 1

    channels = {}
    if len(blocks) == 1:
        channels[0] = blocks[0]
        return channels

    if len(blocks) == lmax + 1:
        channels[lmax] = blocks[0]
        for l in range(lmax):
            channels[l] = blocks[l + 1]
        return channels

    channels[len(blocks) - 1] = blocks[0]
    for idx, terms in enumerate(blocks[1:]):
        channels[idx] = terms
    return channels


def fetch_text_from_url(url):
    try:
        with request.urlopen(url) as response:
            return response.read().decode("utf-8")
    except error.HTTPError:
        raise
    except error.URLError as exc:
        raise RuntimeError(f"Could not reach {url}: {exc.reason}") from exc


def fetch_github_repo_tree(api_url):
    global _GITHUB_TREE_CACHE
    if _GITHUB_TREE_CACHE is not None:
        return _GITHUB_TREE_CACHE

    payload = fetch_text_from_url(api_url)
    parsed = json.loads(payload)
    tree = parsed.get("tree")
    if not isinstance(tree, list):
        raise RuntimeError("Unexpected GitHub API response while listing CHAMP repository files")
    _GITHUB_TREE_CACHE = tree
    return tree


def discover_ecp_path_in_champ(element, github_tree_api_url):
    tree = fetch_github_repo_tree(github_tree_api_url)
    candidates = []

    preferred_1 = f"pool/BFD/ECP_gamess/{element}"
    preferred_2 = f"pool/BFD/ECP_champ/BFD.gauss_ecp.dat.{element}"
    preferred_3 = f"pool/BFD/BFD.gauss_ecp.dat.{element}"
    fallback_name_1 = f"BFD.gauss_ecp.dat.{element}"
    fallback_name_2 = f"ECP.gauss_ecp.dat.{element}"

    for node in tree:
        path = node.get("path", "")
        if not isinstance(path, str):
            continue
        if (
            path.endswith("/" + element)
            or path.endswith(fallback_name_1)
            or path.endswith(fallback_name_2)
        ):
            score = 1000
            if path == preferred_1:
                score = 0
            elif path == preferred_2:
                score = 10
            elif path == preferred_3:
                score = 20
            elif "/pool/" in path:
                score = 100
            elif "/tests/" in path:
                score = 200
            score += len(path)
            candidates.append((score, path))

    if not candidates:
        return None

    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def fetch_ecp_text(element, base_url, github_tree_api_url):
    direct_candidates = [
        f"{base_url.rstrip('/')}/{element}",
        f"{base_url.rstrip('/')}/BFD.gauss_ecp.dat.{element}",
        f"{base_url.rstrip('/')}/ECP.gauss_ecp.dat.{element}",
        f"https://raw.githubusercontent.com/filippi-claudia/champ/main/pool/BFD/ECP_gamess/{element}",
        f"https://raw.githubusercontent.com/filippi-claudia/champ/main/pool/BFD/BFD.gauss_ecp.dat.{element}",
    ]

    errors_seen = []
    for url in direct_candidates:
        try:
            return fetch_text_from_url(url)
        except error.HTTPError as exc:
            errors_seen.append(f"{url} -> HTTP {exc.code}")

    discovered_path = discover_ecp_path_in_champ(element, github_tree_api_url)
    if discovered_path is None:
        detail = "; ".join(errors_seen) if errors_seen else "no candidates tested"
        raise RuntimeError(
            f"Could not find ECP file for element {element} in CHAMP repository ({detail})"
        )

    discovered_url = (
        "https://raw.githubusercontent.com/filippi-claudia/champ/main/" + discovered_path
    )
    try:
        return fetch_text_from_url(discovered_url)
    except error.HTTPError as exc:
        raise RuntimeError(
            "Discovered ECP path but failed to download for element "
            f"{element}: {discovered_url} (HTTP {exc.code})"
        ) from exc


def read_ecp_text_from_dir(element, source_dir):
    candidates = [
        os.path.join(source_dir, element),
        os.path.join(source_dir, f"BFD.gauss_ecp.dat.{element}"),
        os.path.join(source_dir, f"ECP.gauss_ecp.dat.{element}"),
    ]
    path = next((candidate for candidate in candidates if os.path.isfile(candidate)), None)
    if path is None:
        raise RuntimeError(
            f"Missing ECP file for element {element}. Checked: {', '.join(candidates)}"
        )
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()


def load_trexio_data(trexio_file_name):
    try:
        import trexio
    except ImportError as exc:
        raise RuntimeError("Unable to import trexio. Please install TREXIO Python bindings.") from exc

    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Unable to import numpy. Please install numpy.") from exc

    file_handle = trexio.File(trexio_file_name, mode="r", back_end=trexio.TREXIO_HDF5)
    try:
        nucleus_labels_raw = trexio.read_nucleus_label(file_handle)
        ecp_nucleus_index = list(trexio.read_ecp_nucleus_index(file_handle))
        ecp_ang_mom = list(trexio.read_ecp_ang_mom(file_handle))
        ecp_power = list(trexio.read_ecp_power(file_handle))
        ecp_coefficient = np.array(trexio.read_ecp_coefficient(file_handle), dtype=float)
        ecp_exponent = np.array(trexio.read_ecp_exponent(file_handle), dtype=float)
    finally:
        file_handle.close()

    nucleus_labels = [normalize_element(label) for label in nucleus_labels_raw]

    if not (
        len(ecp_nucleus_index)
        == len(ecp_ang_mom)
        == len(ecp_power)
        == len(ecp_coefficient)
        == len(ecp_exponent)
    ):
        raise RuntimeError("TREXIO ECP arrays have inconsistent lengths.")

    return {
        "nucleus_labels": nucleus_labels,
        "ecp_nucleus_index": ecp_nucleus_index,
        "ecp_ang_mom": ecp_ang_mom,
        "ecp_power": ecp_power,
        "ecp_coefficient": ecp_coefficient,
        "ecp_exponent": ecp_exponent,
        "trexio_module": trexio,
        "numpy_module": np,
    }


def build_updated_arrays(data, ecp_data_by_element):
    np = data["numpy_module"]

    labels = data["nucleus_labels"]
    ecp_nucleus_index = data["ecp_nucleus_index"]
    ecp_ang_mom = data["ecp_ang_mom"]
    ecp_power = data["ecp_power"]

    coeff_new = np.array(data["ecp_coefficient"], dtype=float)
    expo_new = np.array(data["ecp_exponent"], dtype=float)

    positions_by_atom = defaultdict(list)
    for pos, atom_idx in enumerate(ecp_nucleus_index):
        positions_by_atom[int(atom_idx)].append(pos)

    for atom_idx, positions in positions_by_atom.items():
        element = labels[atom_idx]
        if element not in ecp_data_by_element:
            raise RuntimeError(f"No ECP data loaded for element {element}")

        channel_terms = ecp_data_by_element[element]

        existing_channel_counts = Counter(ecp_ang_mom[pos] for pos in positions)
        target_channel_counts = {ang: len(terms) for ang, terms in channel_terms.items()}
        if dict(existing_channel_counts) != target_channel_counts:
            raise RuntimeError(
                "Mismatch in ECP channel term counts for atom "
                f"{atom_idx} ({element}). Existing={dict(existing_channel_counts)}, "
                f"target={target_channel_counts}"
            )

        terms_by_channel_and_power = {
            ang: defaultdict(deque) for ang in channel_terms
        }
        for ang, terms in channel_terms.items():
            for power, coeff, expo in terms:
                terms_by_channel_and_power[ang][power].append((coeff, expo))

        for pos in positions:
            ang = ecp_ang_mom[pos]
            power = ecp_power[pos]
            channel_map = terms_by_channel_and_power.get(ang)
            if channel_map is None:
                raise RuntimeError(
                    f"No ECP channel l={ang} for atom {atom_idx} ({element})"
                )
            queue = channel_map.get(power)
            if queue is None or len(queue) == 0:
                available = sorted(channel_map.keys())
                raise RuntimeError(
                    "Could not map ECP term for atom "
                    f"{atom_idx} ({element}), l={ang}, power={power}. "
                    f"Available powers in source for l={ang}: {available}"
                )
            coeff, expo = queue.popleft()
            coeff_new[pos] = coeff
            expo_new[pos] = expo

        leftovers = []
        for ang, power_map in terms_by_channel_and_power.items():
            for power, queue in power_map.items():
                if queue:
                    leftovers.append((ang, power, len(queue)))
        if leftovers:
            raise RuntimeError(
                f"Unused source ECP terms for atom {atom_idx} ({element}): {leftovers}"
            )

    return coeff_new, expo_new


def write_updated_arrays(trexio_file_name, coeff_new, expo_new, trexio):
    file_handle = trexio.File(trexio_file_name, mode="u", back_end=trexio.TREXIO_HDF5)
    try:
        trexio.write_ecp_coefficient(file_handle, coeff_new)
        trexio.write_ecp_exponent(file_handle, expo_new)
    finally:
        file_handle.close()


def verify_arrays(trexio_file_name, coeff_ref, expo_ref, trexio, np):
    file_handle = trexio.File(trexio_file_name, mode="r", back_end=trexio.TREXIO_HDF5)
    try:
        coeff_check = np.array(trexio.read_ecp_coefficient(file_handle), dtype=float)
        expo_check = np.array(trexio.read_ecp_exponent(file_handle), dtype=float)
    finally:
        file_handle.close()

    coeff_ok = np.allclose(coeff_check, coeff_ref, rtol=0.0, atol=0.0)
    expo_ok = np.allclose(expo_check, expo_ref, rtol=0.0, atol=0.0)
    return coeff_ok and expo_ok


def run(args):
    if not os.path.isfile(args.trexio_file):
        raise RuntimeError(f"TREXIO file does not exist: {args.trexio_file}")

    data = load_trexio_data(args.trexio_file)

    labels = data["nucleus_labels"]
    atoms_with_ecp = sorted(set(data["ecp_nucleus_index"]))
    unique_elements = sorted({labels[idx] for idx in atoms_with_ecp})

    if args.verbose:
        print(f"Detected {len(labels)} atoms")
        print(f"Atoms with ECP terms: {len(atoms_with_ecp)}")
        print(f"Elements requiring ECP data: {', '.join(unique_elements)}")

    ecp_data_by_element = {}
    for element in unique_elements:
        if args.source_dir:
            text = read_ecp_text_from_dir(element, args.source_dir)
        else:
            text = fetch_ecp_text(element, args.base_url, args.github_tree_api_url)
        ecp_data_by_element[element] = parse_bfd_ecp_text(text, element)

    coeff_new, expo_new = build_updated_arrays(data, ecp_data_by_element)

    if args.dry_run:
        print("Dry run: computed updated ECP arrays but did not write to file.")
        return

    write_updated_arrays(
        args.trexio_file,
        coeff_new,
        expo_new,
        data["trexio_module"],
    )

    verified = verify_arrays(
        args.trexio_file,
        coeff_new,
        expo_new,
        data["trexio_module"],
        data["numpy_module"],
    )
    if not verified:
        raise RuntimeError("Verification failed: written ECP arrays differ from expected values")

    print("ECP coefficients and exponents updated successfully.")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Automatically replace TREXIO ECP coefficients/exponents with high-precision "
            "BFD values from CHAMP ECP libraries (default: pool/BFD/ECP_gamess) "
            "based on atom labels and ECP layout."
        )
    )
    parser.add_argument("trexio_file", help="Path to TREXIO HDF5 file (e.g. ala3.hdf5)")
    parser.add_argument(
        "--base-url",
        default=DEFAULT_ECP_BASE_URL,
        help=(
            "Base URL for BFD ECP files (default: CHAMP ECP_gamess raw path). "
            "Ignored when --source-dir is used."
        ),
    )
    parser.add_argument(
        "--source-dir",
        default=None,
        help=(
            "Optional local directory containing GAMESS-style files (e.g. O, C, H) "
            "or CHAMP-style files (BFD.gauss_ecp.dat.<Element>). "
            "Use this for offline runs."
        ),
    )
    parser.add_argument(
        "--github-tree-api-url",
        default=DEFAULT_GITHUB_TREE_API_URL,
        help=(
            "GitHub API URL used as fallback to discover ECP files when direct URLs fail."
        ),
    )
    parser.add_argument("--dry-run", action="store_true", help="Parse and map only; do not write")
    parser.add_argument("--verbose", action="store_true", help="Print progress information")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
