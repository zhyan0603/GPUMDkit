#!/usr/bin/env python3
import argparse
import os
import re
import tempfile

import dpdata
from dpdata.system import MultiSystems


def normalize_xyz(content):
    content = re.sub(
        r"(?<![\w])(Energy|Virial)(?==)",
        lambda m: m.group(1).lower(),
        content,
    )
    content = re.sub(
        r"Properties=\S+",
        lambda m: m.group(0).replace("forces", "force"),
        content,
    )
    return content


def detect_type_map(content):
    type_map = []
    for line in content.splitlines():
        if re.match(r"^\s*\d+\s*$", line):
            continue
        if "Properties=" in line or "Lattice=" in line:
            continue
        elem = line.split()[0]
        if elem not in type_map:
            type_map.append(elem)
    return type_map


def main():
    parser = argparse.ArgumentParser(
        description="Convert a labeled extended XYZ file to DeePMD-kit npy format"
    )
    parser.add_argument("xyz_file", help="input labeled extended XYZ file")
    parser.add_argument(
        "-o", "--output", default="deepmd_data", help="output directory (default: deepmd_data)"
    )
    parser.add_argument(
        "-t",
        "--type-map",
        nargs="*",
        default=None,
        help="element type order, e.g. -t Cl Li Y Mg Ca (default: auto-detected from file)",
    )
    parser.add_argument(
        "-f",
        "--from-format",
        default="extxyz",
        help="input format name (default: extxyz)",
    )
    args = parser.parse_args()

    with open(args.xyz_file) as f:
        content = f.read()
    content = normalize_xyz(content)

    type_map = args.type_map if args.type_map else detect_type_map(content)
    print("type_map:", type_map)

    fd, tmp = tempfile.mkstemp(suffix=".xyz")
    with os.fdopen(fd, "w") as f:
        f.write(content)

    try:
        ms = MultiSystems.from_file(tmp, fmt=args.from_format)
    finally:
        os.unlink(tmp)

    print("systems:", len(ms.systems))
    for key, s in ms.systems.items():
        print(
            "  {}: nframes={} natoms={} names={}".format(
                key, s.get_nframes(), s.get_natoms(), s.get_atom_names()
            )
        )
    s0 = list(ms.systems.values())[0]
    print("has virials:", "virials" in s0.data)

    ms.to_deepmd_npy(args.output, type_map=type_map)
    print("DONE ->", args.output)
    print(
        "NOTE: use the order in {}/<system>/type_map.raw as input.json type_map "
        "and for the ZBL pair table".format(args.output)
    )


if __name__ == "__main__":
    main()
