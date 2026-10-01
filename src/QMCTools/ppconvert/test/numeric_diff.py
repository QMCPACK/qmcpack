#!/usr/bin/env python3

"""Compare text files with a configurable absolute tolerance for numeric values.

All non-numeric text must match exactly.

Example
-------
>>> python3 numeric_diff.py --absolute-tolerance 0.001 actual.xml reference.xml
"""

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Literal


NUMBER_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def tokenize(text: str) -> tuple[list[str], list[float]]:
    """Split text around numeric tokens and convert those tokens to floats."""
    text_parts = []
    numbers = []
    position = 0
    for match in NUMBER_PATTERN.finditer(text):
        text_parts.append(text[position : match.start()])
        numbers.append(float(match.group()))
        position = match.end()
    text_parts.append(text[position:])
    return text_parts, numbers


def compare_files(
    actual_path: Path, reference_path: Path, absolute_tolerance: float
) -> bool:
    """Return whether two files match within the numeric absolute tolerance."""
    actual_parts, actual_numbers = tokenize(actual_path.read_text())
    reference_parts, reference_numbers = tokenize(reference_path.read_text())

    if actual_parts != reference_parts:
        print(
            f"Non-numeric content differs: {actual_path} != {reference_path}",
            file=sys.stderr,
        )
        return False

    if len(actual_numbers) != len(reference_numbers):
        print(
            f"Numeric value counts differ: {len(actual_numbers)} != {len(reference_numbers)}",
            file=sys.stderr,
        )
        return False

    for index, (actual, reference) in enumerate(
        zip(actual_numbers, reference_numbers, strict=True), start=1
    ):
        difference = abs(actual - reference)
        if (
            not math.isfinite(actual)
            or not math.isfinite(reference)
            or difference > absolute_tolerance
        ):
            print(
                f"Numeric value {index} differs: actual={actual:.17g}, reference={reference:.17g}, "
                f"absolute difference={difference:.17g}, tolerance={absolute_tolerance:.17g}",
                file=sys.stderr,
            )
            return False

    return True


def main() -> Literal[0, 1, 2]:
    """Parse command-line arguments and return a process exit status."""
    parser = argparse.ArgumentParser(
        description="Compare text files with an absolute numeric tolerance."
    )
    parser.add_argument("actual", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--absolute-tolerance", type=float, required=True)
    args = parser.parse_args()

    if not math.isfinite(args.absolute_tolerance) or args.absolute_tolerance < 0.0:
        parser.error("--absolute-tolerance must be a finite, non-negative value")

    try:
        matches = compare_files(args.actual, args.reference, args.absolute_tolerance)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"Unable to compare files: {error}", file=sys.stderr)
        return 2

    return 0 if matches else 1


if __name__ == "__main__":
    sys.exit(main())
