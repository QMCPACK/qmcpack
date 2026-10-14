#!/usr/bin/env python3

import math
import sys
import xml.etree.ElementTree as ET

# Check optimized wavefunction parameters against a reference.
# Usage: check_h4_descent.py OUTPUT_XML REFERENCE_XML
# OUTPUT_XML is the optimized wavefunction produced by QMCPACK.
# REFERENCE_XML contains the expected optimized parameter values.

TOLERANCE = 1.0e-6
EXPECTED_PARAMETER_COUNT = 31


def read_parameters(filename):
    root = ET.parse(filename).getroot()
    parameters = {}

    for coefficients in root.findall("wavefunction/jastrow/correlation/coefficients"):
        parameters[coefficients.attrib["id"]] = [
            float(value) for value in coefficients.text.split()
        ]

    for csf in root.findall("wavefunction/determinantset/multideterminant/detlist/csf"):
        parameters[csf.attrib["id"]] = [float(csf.attrib["coeff"])]

    return parameters


def main():
    output = read_parameters(sys.argv[1])
    reference = read_parameters(sys.argv[2])
    reference_count = sum(len(values) for values in reference.values())

    if reference_count != EXPECTED_PARAMETER_COUNT:
        print(
            f"Expected {EXPECTED_PARAMETER_COUNT} reference parameters, found {reference_count}"
        )
        return 1

    for name, expected_values in reference.items():
        actual_values = output.get(name)
        if actual_values is None:
            print(f"Missing optimized parameter group: {name}")
            return 1
        if len(actual_values) != len(expected_values):
            print(
                f"Parameter count differs for {name}: {len(actual_values)} != {len(expected_values)}"
            )
            return 1

        for index, (actual, expected) in enumerate(zip(actual_values, expected_values)):
            if not math.isfinite(actual) or abs(actual - expected) > TOLERANCE:
                print(f"Parameter {name}[{index}] differs: {actual} != {expected}")
                return 1

    print(f"Validated {reference_count} optimized parameters")
    return 0


if __name__ == "__main__":
    sys.exit(main())
