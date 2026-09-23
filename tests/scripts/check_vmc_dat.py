#! /usr/bin/env python3

"""Check weighted LocalEnergy means written to VMC per-step data files."""

import argparse
import math
import sys


def read_vmc_dat(filename, steps):
    with open(filename, encoding="utf-8") as vmc_file:
        header = vmc_file.readline().split()
        if len(header) < 3 or header[0] != "#":
            raise ValueError("missing vmc.dat header")
        columns = header[2:]
        try:
            energy_index = columns.index("LocalEnergy")
            weight_index = columns.index("BlockWeight")
        except ValueError as error:
            raise ValueError("vmc.dat lacks LocalEnergy or BlockWeight") from error

        values_by_step = []
        for line in vmc_file:
            fields = line.split()
            if not fields:
                continue
            values = fields[1:]
            if len(values) != len(columns):
                raise ValueError("malformed vmc.dat row")
            values_by_step.append((float(values[energy_index]), float(values[weight_index])))
    if len(values_by_step) < steps:
        raise ValueError("vmc.dat has fewer rows than one block")
    weighted_energy = sum(energy * weight for energy, weight in values_by_step[-steps:])
    total_weight = sum(weight for _, weight in values_by_step[-steps:])
    if total_weight == 0.0:
        raise ValueError("vmc.dat has zero total weight")
    return weighted_energy / total_weight


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix")
    parser.add_argument("checks", nargs="+", metavar="SERIES:ENERGY:TOLERANCE:STEPS",
                        help="one final-block LocalEnergy mean per VMC series")
    args = parser.parse_args()

    failed = False
    for check in args.checks:
        try:
            series, expected, tolerance, steps = check.split(":")
            actual = read_vmc_dat("{}.s{:03d}.vmc.dat".format(args.prefix, int(series)), int(steps))
            expected = float(expected)
            tolerance = float(tolerance)
        except (OSError, ValueError) as error:
            print("{}: {}".format(check, error))
            failed = True
            continue

        difference = abs(actual - expected)
        print("series {}: vmc.dat LocalEnergy = {:.10f}, expected = {:.10f}".format(series, actual, expected))
        if not math.isfinite(actual) or difference > tolerance:
            print("series {}: difference {:.3e} exceeds {:.3e}".format(series, difference, tolerance))
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
