#!/usr/bin/env python3
#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//////////////////////////////////////////////////////////////////////////////////////
"""Check that the real-checkpoint smoke test produced finite scalar data."""
import math
import pathlib
import sys

workdir = pathlib.Path(sys.argv[1])
files = list(workdir.glob("*.scalar.dat"))
if len(files) != 1:
    raise SystemExit(f"expected one scalar.dat output in {workdir}, found {len(files)}")
contents = files[0].read_text().splitlines()
header_index = next((i for i, line in enumerate(contents) if "LocalEnergy" in line), None)
if header_index is None:
    raise SystemExit("scalar output header lacks LocalEnergy")
header = contents[header_index].lstrip("# ").split()
rows = [line.split() for line in contents[header_index + 1:] if line.strip() and not line.lstrip().startswith("#")]
if len(rows) != 2:
    raise SystemExit(f"expected exactly two data rows, found {len(rows)}")
for row in rows:
    if len(row) != len(header):
        raise SystemExit("scalar output row width does not match header")
    try:
        values = [float(value) for value in row]
    except ValueError as error:
        raise SystemExit(f"nonnumeric scalar output: {error}") from error
    if not all(math.isfinite(value) for value in values):
        raise SystemExit("scalar output contains a non-finite value")
