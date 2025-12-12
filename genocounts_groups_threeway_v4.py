#!/usr/bin/python

import sys

popfile = sys.argv[1]        # pop/species assignment
genofile = sys.argv[2]       # merged genotype table
min_prop_diff = 0.8          # threshold for AIM calling

# Load populations for species mapping
pop2species = {}
species_list = []

with open(popfile) as f:
    for line in f:
        pop, species = line.strip().split()
        pop2species[pop] = species
        if species not in species_list:
            species_list.append(species)

# Count number of AIMs for each species
highest_alt_counts = {sp: 0 for sp in species_list}

# Print header
cols = ["CHROM","POS","REF","ALT"]
for sp in species_list:
    cols += [f"{sp}_ref", f"{sp}_alt", f"{sp}_na", f"{sp}_nsamples"]
print("\t".join(cols))

# Process genotypes
with open(genofile) as f:
    header = f.readline().strip().split()
    sample_species = [pop2species.get(s, "NA") for s in header[4:]]
    
    for line in f:
        entries = line.strip().split()
        chrom, pos, ref, alt = entries[:4]
        genos = entries[4:]
        
        counts = {sp: {"ref":0,"alt":0,"na":0,"n":0} for sp in species_list}
        
        for g, sp in zip(genos, sample_species):
            counts[sp]["n"] += 1
            if g in ["0/0","0|0"]:
                counts[sp]["ref"] += 2
            elif g in ["0/1","1/0","0|1","1|0"]:
                counts[sp]["ref"] += 1
                counts[sp]["alt"] += 1
            elif g in ["1/1","1|1"]:
                counts[sp]["alt"] += 2
            else:  # ./.
                counts[sp]["na"] += 1
        

        # Compute ALT freq and only keep if >0.9 difference
        
        props = {}
        valid_site = True

        for sp in species_list:
            total = counts[sp]["ref"] + counts[sp]["alt"]
            if total == 0:
                valid_site = False   # Remove site if any species has zero counts
                break
            props[sp] = counts[sp]["alt"] / total

        if valid_site:
            # Count species with alt frequency >= 0.9
            high_alt_species = [sp for sp, p in props.items() if p >= 0.9]

            # Remove site if more than one species is high-alt
            if len(high_alt_species) == 1:
                max_prop = max(props.values())
                min_prop = min(props.values())
        
                # Keep site only if difference >= 0.9
                if (max_prop - min_prop) >= 0.9:
                    winner = high_alt_species[0]
                    highest_alt_counts[winner] += 1
            
                    out = [chrom, pos, ref, alt]
                    for sp in species_list:
                        c = counts[sp]
                        out += [str(c["ref"]), str(c["alt"]), str(c["na"]), str(c["n"])]
                    print("\t".join(out))


# Print number of AIMs for each species
print("Total high-ALT SNPs per species:", highest_alt_counts, file=sys.stderr)

#DONE
