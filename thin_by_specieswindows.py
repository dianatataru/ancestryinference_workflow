#!/usr/bin/python

import sys
import random
from collections import defaultdict

if len(sys.argv) < 3:
    print("Usage: python thin_equal_density.py aims.txt window_size > thinned.txt")
    sys.exit(1)

infile = sys.argv[1]
WINDOW = int(sys.argv[2])

random.seed(42)

# windows[(chrom, window_index)][species] = list of lines
windows = defaultdict(lambda: defaultdict(list))

species_counts_before = defaultdict(int)
species_counts_after = defaultdict(int)

with open(infile) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.strip().split()
        chrom = fields[0]
        pos = int(fields[1])
        species = fields[10]      # last column

        species_counts_before[species] += 1

        win = pos // WINDOW
        windows[(chrom, win)][species].append(line.strip())

output_lines = []

# Process each window
for key, species_dict in windows.items():
    # Count number of sites per species in this window
    counts = {sp: len(species_dict[sp]) for sp in species_dict}

    if not counts:
        continue

    # Minimum count across species in this window
    target = min(counts.values())

    # Randomly select target number of sites per species
    for sp, lines in species_dict.items():
        keep_n = min(target, len(lines))
        kept = random.sample(lines, keep_n)

        for site in kept:
            species_counts_after[sp] += 1
            output_lines.append(site)

# Output thinned dataset
for line in output_lines:
    print(line)

# Report stats to STDERR
sys.stderr.write("\n=== Species counts BEFORE thinning ===\n")
for sp, c in species_counts_before.items():
    sys.stderr.write(f"{sp}\t{c}\n")

sys.stderr.write("\n=== Species counts AFTER thinning ===\n")
for sp, c in species_counts_after.items():
    sys.stderr.write(f"{sp}\t{c}\n")
sys.stderr.write("\n")

