#!/usr/bin/env python3
"""Check every ensemble group's statistics and the existing LiH ion-ion reference."""
import argparse
import math
from pathlib import Path

# Run from the ensemble test's output directory, for example:
#   python3 /path/to/check_ensemble.py --ranks 5 \
#       --prefix det_ensemble_blocks_vmc_msdj --ionion 0.995380119260
# --ranks is the total MPI rank count across all four ensemble groups (4 or 5).
# --prefix is the output project ID before .gNNN.s000.scalar.dat; the checker
# reads groups g000 through g003 relative to the current working directory.
# --ionion is the reference IonIon energy, checked in every block with an
# absolute tolerance of 1e-6.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--ranks', type=int, choices=(4, 5), required=True)
parser.add_argument('--prefix', required=True)
parser.add_argument('--ionion', type=float, required=True)
args = parser.parse_args()

for group, blocks in enumerate((1, 2, 3, 4)):
    path = Path(f'{args.prefix}.g{group:03d}.s000.scalar.dat')
    lines = path.read_text().splitlines()
    names = lines[0].split()[1:]
    rows = [dict(zip(names, map(float, line.split())))
            for line in lines[1:] if line.strip() and not line.startswith('#')]
    # FairDivideLow assigns any extra MPI rank to the last group.
    ranks = 1 + int(args.ranks == 5 and group == 3)
    expected_samples = blocks * 3 * ranks
    if len(rows) != blocks:
        raise SystemExit(f'{path}: expected {blocks} blocks, found {len(rows)}')
    for row in rows:
        if not all(math.isfinite(value) for value in row.values()):
            raise SystemExit(f'{path}: non-finite scalar data')
        if not math.isclose(row['IonIon'], args.ionion, rel_tol=0, abs_tol=1e-6):
            raise SystemExit(f'{path}: IonIon {row["IonIon"]} != {args.ionion}')
    samples = sum(row['BlockWeight'] for row in rows)
    if not math.isclose(samples, expected_samples, rel_tol=0, abs_tol=1e-8):
        raise SystemExit(f'{path}: samples {samples} != {expected_samples}')
    print(f'{path}: {blocks} blocks, {samples:g} samples, IonIon {args.ionion:g}: pass')
