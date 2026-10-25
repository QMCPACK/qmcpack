#!/usr/bin/env python3

# Compare variational-parameter names and values in a *.vp.h5 file with the corresponding *.opt.xml file.
# Report missing, ambiguous, non-finite, or numerically inconsistent parameters as test failures.

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET

import h5py

REL_TOLERANCE = 1.0e-9
ABS_TOLERANCE = 1.0e-12


# Read and validate the ordered parameter name/value pairs stored in variational-parameter HDF5 file.
def read_h5_parameters(filename):
    with h5py.File(filename, "r") as h5_file:
        names_dataset = h5_file["name_value_lists/parameter_names"]
        values_dataset = h5_file["name_value_lists/parameter_values"]

        if names_dataset.ndim != 1 or values_dataset.ndim != 1:
            raise ValueError(
                "parameter_names and parameter_values must be one-dimensional"
            )

        names = names_dataset.asstr()[:].tolist()
        values = values_dataset[:].tolist()

    if len(names) != len(values):
        raise ValueError(
            f"HDF5 parameter name/value counts differ: {len(names)} != {len(values)}"
        )
    if len(names) != len(set(names)):
        raise ValueError("HDF5 parameter names are not unique")

    return list(zip(names, values))


# Convert an XML element's whitespace-separated numeric text to a list of floats, or return an empty list.
def parse_numeric_text(element):
    if element.text is None:
        return []
    try:
        return [float(value) for value in element.text.split()]
    except ValueError:
        return []


# Resolve an HDF5 parameter name to its scalar value in the indexed XML elements.
def find_xml_value(elements_by_id, parameter_name):
    exact_values = []
    for element in elements_by_id.get(parameter_name, []):
        if "coeff" in element.attrib:
            exact_values.append(float(element.attrib["coeff"]))
        else:
            text_values = parse_numeric_text(element)
            if len(text_values) == 1:
                exact_values.append(text_values[0])

    if len(exact_values) == 1:
        return exact_values[0]
    if len(exact_values) > 1:
        raise ValueError(f"XML parameter {parameter_name} is ambiguous")

    indexed_name = re.fullmatch(r"(.+)_(\d+)", parameter_name)
    if indexed_name:
        base_name, index_text = indexed_name.groups()
        index = int(index_text)
        indexed_values = []
        for element in elements_by_id.get(base_name, []):
            text_values = parse_numeric_text(element)
            if index < len(text_values):
                indexed_values.append(text_values[index])

        if len(indexed_values) == 1:
            return indexed_values[0]
        if len(indexed_values) > 1:
            raise ValueError(f"XML parameter {parameter_name} is ambiguous")

    raise ValueError(f"HDF5 parameter {parameter_name} was not found in the XML file")


# Parse an optimized-wavefunction XML file and index all elements carrying an id attribute.
def read_xml_elements(filename):
    root = ET.parse(filename).getroot()
    elements_by_id = {}
    for element in root.iter():
        parameter_id = element.attrib.get("id")
        if parameter_id:
            elements_by_id.setdefault(parameter_id, []).append(element)
    return elements_by_id


# Parse command-line arguments, compare all HDF5/XML parameters, and return the test status.
def main():
    parser = argparse.ArgumentParser(
        description="Compare QMCPACK variational parameters in HDF5 and XML files."
    )
    parser.add_argument("h5_file", help="QMCPACK *.vp.h5 variational-parameter file")
    parser.add_argument(
        "xml_file", help="QMCPACK *.opt.xml optimized-wavefunction file"
    )
    args = parser.parse_args()

    try:
        h5_parameters = read_h5_parameters(args.h5_file)
        xml_elements = read_xml_elements(args.xml_file)

        for name, h5_value in h5_parameters:
            xml_value = find_xml_value(xml_elements, name)
            if not math.isfinite(h5_value) or not math.isfinite(xml_value):
                raise ValueError(f"Parameter {name} contains a non-finite value")
            if not math.isclose(
                h5_value, xml_value, rel_tol=REL_TOLERANCE, abs_tol=ABS_TOLERANCE
            ):
                raise ValueError(
                    f"Parameter {name} differs: HDF5 {h5_value} != XML {xml_value}"
                )
    except (KeyError, OSError, TypeError, ValueError, ET.ParseError) as error:
        print(f"Test status: fail\n{error}")
        return 1

    print(f"Validated {len(h5_parameters)} HDF5/XML variational parameters")
    print("Test status: pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
