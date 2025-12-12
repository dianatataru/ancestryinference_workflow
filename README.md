# ancestryinference_workflow
Alterations to the Schumer Lab's threeway ancestryinfer program (https://github.com/Schumerlab/ancestryinfer) to be run with three Mimulus species (M. guttatus, M. lacinatus, M. nasutus). Requires downloading the ancestryinfer program to run.

Written by: Diana Tataru
Created: June 11, 2025

## 0. Create AIMS: Run 0_genotype_filter_finalDT.sh

This script is adapted for three species from the following manuscript: Fluctuating reproductive isolation and stable ancestry structure in a fine-scaled mosaic of hybridizing Mimulus monkeyflowers" by Matthew Farnitano, Keith Karoly, and Andrea Sweigart, and github: https://github.com/mfarnitano/CAC_popgen/#blob/main/reference_panels/genotype_filter.sh  

These are the populations used for my AIMs:

|    Species    |  Pop |    SRA     | Longitude  | Latitude  |
|    -------    | ---  | ---------  | ---------  | --------  |
|  M. guttatus  |  TOL | Phyt 551   | -120.63315 | 37.969917 |
|  M. guttatus  |  YVO | Unpub.     | -119.74643 | 37.723367 |
|  M. guttatus  | MAR  | SRX030542  | -123.29445 | 43.4786   |
|  M. guttatus  | LMC  | SRX030680  | -123.083917| 38.863983 |
|  M. guttatus  | IM   | SRR398937  | -122.508783|	45.57571 |
|  M. guttatus  | AHQT | SRX142379  | -110.813   | 44.431    |
| M. laciniatus | OPN  | SRR23709136| -119.4852	 | 37.8107   |
| M. laciniatus	| WLF  | Unpub.	    | -119.59385 | 37.841533 |
| M. laciniatus	| TRT  | SRX19570592| -119.70535 |	37.7165  |
| M. laciniatus	| PER  | SRX19570591| -119.3687	 | 37.055767 |
| M. laciniatus	| HUL  | Unpub.	    | -119.150168| 37.2334366|
| M. nasutus	| SF   | SRR29155563| -121.0225	 | 45.264444 |
| M. nasutus	| DPRN | SRR1259273	| -120.344	 | 37.828    |
| M. nasutus	| KOOT | SRR1259272	| -115.983	 | 48.104    |
| M. nasutus	| NHN  | SRX525051	| -124.16	 | 49.273    |
| M. nasutus	| CACN | SRR1259271	| -121.3667	 | 45.71076  |

Located in ```project/dtataru/lac_nas_gut/AIMS```. This script creates a file with ancestry informative sites (AIMs) from 5 allopatric populations each of three species: *M. guttatus, M.laciniatus,* and *M. nasutus*. Input file for it is a joint genotyped vcf ```lacnasgut_jointgeno.vcf.gz```, created using the Ferris Lab GATK variant calling pipeline, with TOLv5 of *M. guttatus* as the reference. This is the current version of the script:

```
#!/bin/bash
#SBATCH --job-name=genotype_filter                     
#SBATCH --output=/project/dtataru/lac_nas_gut/AIMS/logs/genotype_filter_%j.out
#SBATCH --error=/project/dtataru/lac_nas_gut/AIMS/logs/genotype_filter_%j.err
#SBATCH --time=24:00:00 
#SBATCH -p single
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac

###SETUP
PROJECT=lacnasgut
BATCH_NAME=lacnasgut_panel

SCRIPTS_DIR=/project/dtataru/lac_nas_gut/AIMS
LOG_DIR=/project/dtataru/lac_nas_gut/AIMS/logs
WORKING_DIR=/project/dtataru/lac_nas_gut/AIMS
VCF_DIR=/project/dtataru/lac_nas_gut/4_ref/3_Genotyped_GVCFs
TMPDIR=/work/dtataru/AIMS

GENOME=/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa
REPEATMASK=/project/dtataru/lac_nas_gut/AIMS/MguttatusTOL_551_v5.0.repeatmasked_assembly_v5.0.gff3

NTHREADS=32
MEM=120

#inputs
RAW_VCF=${VCF_DIR}/lacnasgut_jointgeno.vcf.gz
GROUP=panel15
REFCODE=TOL551
PREFIX=${WORKING_DIR}/VCFs/${GROUP}.${REFCODE}

cd $WORKING_DIR

###MODULES
module load gatk/4.5.0.0

###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/chr_positions.list ]; then
       head -n14 ${GENOME}.fai | awk '{print $1 ":1-"$2}' > ${WORKING_DIR}/chr_positions.list
fi

printf "\n...extracting biallelic SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
       -select-type SNP --restrict-alleles-to BIALLELIC -O ${PREFIX}.SNPs.vcf.gz

printf "\n...extracting invariant sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
       -select-type NO_VARIATION -O ${PREFIX}.INVTs.vcf.gz

printf "\n...filtering SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" VariantFiltration -V ${PREFIX}.SNPs.vcf.gz -O ${PREFIX}.SNPs.f.vcf.gz \
       -filter "QD < 2.0" --filter-name "QD2" \
       -filter "QUAL < 40.0" --filter-name "QUAL40" \
       -filter "SOR > 3.0" --filter-name "SOR4" \
       -filter "FS > 60.0" --filter-name "FS60" \
       -filter "MQ < 40.0" --filter-name "MQ40" \
       -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
       -filter "ReadPosRankSum < -12.5" --filter-name "ReadPosRankSum-12.5" \
       -filter "ReadPosRankSum > 12.5" --filter-name "ReadPosRankSum12.5" \
       --verbosity ERROR

printf "\n...filtering INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" VariantFiltration -V ${PREFIX}.INVTs.vcf.gz -O ${PREFIX}.INVTs.f.vcf.gz \
       -filter "QD < 2.0" --filter-name "QD2" \
       -filter "SOR > 3.0" --filter-name "SOR4" \
       -filter "MQ < 40.0" --filter-name "MQ40" \
       --verbosity ERROR

printf "\n...sorting SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
       -I ${PREFIX}.SNPs.f.vcf.gz -use_jdk_inflater --TMP_DIR ${TMPDIR} \
       -SD ${GENOME}.dict -O ${PREFIX}.SNPs.fs.vcf.gz

printf "\n...sorting INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
       -I ${PREFIX}.INVTs.f.vcf.gz -SD ${GENOME}.dict -O ${PREFIX}.INVTs.fs.vcf.gz

printf "\n...selecting passing SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
       -V ${PREFIX}.SNPs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.SNPs.fsp.vcf.gz

printf "\n...selecting passing INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
       -V ${PREFIX}.INVTs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.INVTs.fsp.vcf.gz

module purge
module load vcftools/0.1.16
module load htslib/1.21

printf "\n...prepare bed file of repeat-masked regions\n" | tee >(cat >&2)
printf '#chrom\tchromStart\tchromEnd\n' > ${REPEATMASK}.bed
cut -f1,4,5 ${REPEATMASK} >> ${REPEATMASK}.bed

printf "\n...filter individual SNP genotypes by depth and GQ, and exclude repeat-masked regions\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.SNPs.fsp.vcf.gz -c --minGQ 15 --minDP 6 --maxDP 100 \
       --exclude-bed ${REPEATMASK}.bed \
       --recode --recode-INFO-all | bgzip -c > ${PREFIX}.SNPs.fspi.rm.vcf.gz

printf "\n...filter individual INVT genotypes by depth\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.INVTs.fsp.vcf.gz -c --minDP 6 --maxDP 100 \
       --exclude-bed ${REPEATMASK}.bed \
       --recode --recode-INFO-all | bgzip -c > ${PREFIX}.INVTs.fspi.rm.vcf.gz

tabix -p vcf ${PREFIX}.SNPs.fspi.rm.vcf.gz
tabix -p vcf ${PREFIX}.INVTs.fspi.rm.vcf.gz

module purge
module load gatk/4.5.0.0

printf "\n...merge SNP and INVT sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" MergeVcfs \
       -I ${PREFIX}.SNPs.fspi.rm.vcf.gz -I ${PREFIX}.INVTs.fspi.rm.vcf.gz \
       -O ${PREFIX}.merged.fspi.rm.vcf.gz

printf "\n...sort merged VCF\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
       -I ${PREFIX}.merged.fspi.rm.vcf.gz \
       -SD ${GENOME}.dict -O ${PREFIX}.merged.fspi.rm.vcf.gz

printf "\n...filter SNPs by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.SNPs.fspi.rm.vcf.gz \
       --max-nocall-number 7 --exclude-filtered -O ${PREFIX}.SNPs.fspi.rm.31called.vcf.gz

printf "\n...filter merged VCF by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.merged.fspi.rm.vcf.gz \
       --max-nocall-number 7 --exclude-filtered -O ${PREFIX}.merged.fspi.rm.31called.vcf.gz

printf "\n...making genotype tables\n" | tee >(cat >&2)
module purge
module load bcftools/1.18

printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.SNPs.fspi.rm.table
printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.merged.fspi.rm.table

bcftools query -l ${PREFIX}.SNPs.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.SNPs.fspi.rm.table
bcftools query -l ${PREFIX}.merged.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.merged.fspi.rm.table

printf '\n' >> ${PREFIX}.SNPs.fspi.rm.table
printf '\n' >> ${PREFIX}.merged.fspi.rm.table

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.SNPs.fspi.rm.vcf.gz >> ${PREFIX}.SNPs.fspi.rm.table
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.merged.fspi.rm.vcf.gz >> ${PREFIX}.merged.fspi.rm.table

### filtering, edits starting here made by Diana Tataru, December 2025 ###

module load python3

python ${SCRIPTS_DIR}/genocounts_groups_threeway_v4.py ${WORKING_DIR}/lacnasgut_pops.txt ${PREFIX}.SNPs.fspi.rm.table > ${PREFIX}.SNPs.fspi.rm.counts.v5.txt

awk -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,$5,$6,$9,$10,$13,$14}' \
    ${PREFIX}.SNPs.fspi.rm.counts.v5.txt \
    > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.txt

printf "\n...make a column with species classification\n" | tee >(cat >&2)
awk 'NR>1 {
    g=$6; n=$8; l=$10;

    max=g; species="guttatus";

    if (n>max) {max=n; species="nasutus"}
    if (l>max) {species="laciniatus"}

    print $0 "\t" species
}' panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v5.txt \
   > panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v5.species.txt

printf "\n...balance number of aims by species density in 100kb windows\n" | tee >(cat >&2)
python ${SCRIPTS_DIR}/thin_by_specieswindows.py ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.species.txt 100000 \
    >  ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensity.txt \
    2> thinning.stats

sort -k1,1V -k2,2n ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensity.txt > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt

printf "\n...edited to remove underscore in chromosome\n" | tee >(cat >&2)
awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $3, $4}' \
    ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.v3.txt

awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $5, $6, $7, $8, $9, $10}' \
    ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.v3.txt
```
These are the two scripts called above, genocounts_groups_threeway_v4.py:

```
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
```
and thin_by_specieswindows.py:

```
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
```

Additional files needed to run this:
Located in ```/project/dtataru/lac_nas_gut/AIMS```. This has multiple input files that need to be created and changed.   
  a.  output ```AIMs_panel15_final.AIMs_counts.txt``` & ```AIMs_panel15_final.AIMs.txt``` 
  b.  Reference genomes for each of three species, I'm using the following, listed in ```/project/dtataru/hybrids/ancestryinfer/reference_genomes/```:  
      1.  MguttatusTOL_551_v5.0.fa  
      2.  Mnasutusvar_SF_822_v2.0.fa  
      3.  WLF47.fasta (consensus genome made by me with high coverage unpub. sequencing data  

To visualize distribution of AIMs genomewide, run ```visualizeAIMs.R ```. This also outputs number of AIMS per chromosome. This is aims_v1:

| Chr01 | Chr02 | Chr03 | Chr04 | Chr05 | Chr06 | Chr07 | Chr08 | Chr09 | Chr10 | Chr11 | Chr12 | Chr13 | Chr14 |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
|  9260 | 15019 | 11090 | 17851 | 13041 | 16073 | 8709  | 20077 | 10296 | 15111 | 11289 | 14642 | 15133 | 25823 | 

number of AIMS per chromosome for aims_v2:

| Chr01 | Chr02 | Chr03 | Chr04 | Chr05 | Chr06 | Chr07 | Chr08 | Chr09 | Chr10 | Chr11 | Chr12 | Chr13 | Chr14 |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
|33172  |55669  |41019  |64583  |48556  |60589  | 33813 | 74756 | 41348 | 57280 | 43404 | 56877 | 54787 | 95319 | 

### Editing different versions of AIMs (currently using aims_v5)

For aims_v1, ended up with 698623 remaining_sites.txt, then after filtering to sites uniquely fixed to each species and only one AIM for every 100 bp, ended up with 203414 sites. Need to revisit this filtering step because I'm worries that those many sites that weren't fixed were mostly guttatus and then when I subset I ended up with very few laciniatus sites. I think it has to do with the way sites were subset or some way that ref/alt is being called. I think I want to try to run this with the full set up unthinned AIMS that I got, with 1267995 total sites. After running this it is still off.

I found a typo in my script which meant that nasutus counts were being called incorrectly. Also, it was not picking up whether the alt allele was the one identified in the AIMS, so I changed code to only match counts when the ref and alt from AIMs match those in the vcf. Lastly, I edited it so that when genotype call is: 0/0 = 2 0, 0/1 = 1 1, 1/1 = 0 2. I fixed these and reran (panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v2.filtered.thinned.txt), ended up with 800357 sites, and am now going to rerun them through steps 3 and 4, with 4 having different HMM priors. 

I reran it and everything came up looking like 50% laciniatus, 50% nasutus. I think it has to do with this filtering line of code in my 0_genotype_filter_finalDT.sh script:
```
awk '($5 != 0 || $7 != 0 || $9 != 0)' ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v2.txt > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v2.filtered.txt
```
The goal of this line is to remove sites where all sites have the alternate allele, by removing sites where all species REF counts are 0, but that is not actually what I want. I actually want to just keep sites where the alt is unique to the species. 

Additionally, these lines of code in genocounts_group_threeway.py are trying to filter but incorrectly because it doesn't account for the third species proportion:
```
# Check that each species has at least one called allele
        if all((counts[sp]["ref"] + counts[sp]["alt"]) > 0 for sp in species_list):
            # Calculate alt allele frequencies
            props = {sp: counts[sp]["alt"] / max(1, counts[sp]["ref"]+counts[sp]["alt"]) for sp in species_list}
            max_prop = max(props.values()) #highest alt proportion across species
            min_prop = min(props.values()) #lowest alt proporiton across species
            if (max_prop - min_prop) >= min_prop_diff:
                out = [chrom,pos,ref,alt]
                for sp in species_list:
                    c = counts[sp]
                    out += [str(c["ref"]), str(c["alt"]), str(c["na"]), str(c["n"])]
                print("\t".join(out))
```
We want to do something similar to banarjee et al. 2023, "We thinned to an approximately equivalent number of informative sites between all pairs of species. To do so, we retained all ancestry-informative sites that distinguished X. birchmanni and X. malinche, and every other site that distinguished X. variatus from either of these two species." To do so I am changing that code to:
```
    # Compute ALT proportions for every species
    props = {
        sp: counts[sp]["alt"] / (counts[sp]["ref"] + counts[sp]["alt"])
        for sp in species_list
    }

    # Count how many species have alt freq > 0.8
    high_alt_species = [sp for sp, p in props.items() if p >= 0.8]
	for sp in high_alt_species:
    	high_alt_counts[sp] += 1

    # keep site if only one species has alt proportion higher than 0.8
    if len(high_alt_species) == 1:
        out = [chrom, pos, ref, alt]
        for sp in species_list:
            c = counts[sp]
            out += [str(c["ref"]), str(c["alt"]), str(c["na"]), str(c["n"])]
        print("\t".join(out))
		print("Total high-ALT SNPs per species:", high_alt_counts, file=sys.stderr)
```

panel15.TOL551.SNPs.fspi.rm.counts.v3.txt: 2,746,758 total sites
Total high-ALT SNPs per species: {'guttatus': 658786, 'nasutus': 1076451, 'laciniatus': 1011520}
Output thinned: 1,104,019 total sites

Okay, so now I need to figure out how many of these thinned sites are max proportion ALT for each species, and thin to about equal numbers for each. It's still a little unclear to me whether it's better to filter sites before or after thinning. If I filtered before thinning then I would maybe have to thin again.

```
awk 'NR>1 {
    max=$6; winner="guttatus";
    if ($8>max) {max=$8; winner="nasutus"}
    if ($10>max) {max=$10; winner="laciniatus"}
    count[winner]++
}
END {
    for (sp in count) print sp, count[sp]
}' panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v3.thinned.txt
```
laciniatus 377015
nasutus 427369
guttatus 299634

Not horribly different between species. So, the main difference between these aims_v3 and aims_v2 are that these ones contain sites where some species have all 0s. I think I do actually want filter by sites where all species have at least one count in ref or alt, and call that aims_v4.

thinned.aims_v4: laciniatus 287938, nasutus 348803, guttatus 131002
unthinned.aims_v4: {'guttatus': 257533, 'nasutus': 794065, 'laciniatus': 695118}

I also want to do some species-aware balancing and then thinning to the unthinned.aims_v4. First I'll create a column to classify the species:

```
awk 'NR>1 {
    g=$6; n=$8; l=$10;

    max=g; species="guttatus";
    if (n>max) {max=n; species="nasutus"}
    if (l>max) {species="laciniatus"}

    print $0 "\t" species
}' panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.txt \
    > panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.species.txt
```
Then balance by species in balance_species.py with
```python balance_species.py /project/dtataru/lac_nas_gut/AIMS/VCFs/panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.species.txt /project/dtataru/lac_nas_gut/AIMS/VCFs/panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.txt 257000``` 
where 257000 is max_per_species I want to keep. Here is that script:

```
#!/usr/bin/python
import sys, random

###USAGE:python balance_species.py {AIM_counts}.species.txt {AIM_counts}.speciesbalanced.txt max_per_species

infile = sys.argv[1]
outfile = sys.argv[2]
max_per_species = int(sys.argv[3])  # e.g. 200000

species_rows = {"guttatus": [], "nasutus": [], "laciniatus": []}

with open(infile) as f:
    header = f.readline()
    for line in f:
        fields = line.strip().split("\t")
        species = fields[-1]
        species_rows[species].append(line)

random.seed(42)

balanced = []
for sp in species_rows:
    rows = species_rows[sp]
    if len(rows) > max_per_species:
        balanced += random.sample(rows, max_per_species)
    else:
        balanced += rows

with open(outfile, "w") as out:
    out.write(header)
    for row in balanced:
        out.write(row)

for sp in species_rows:
    print(sp, "kept", min(len(species_rows[sp]), max_per_species),
          "of", len(species_rows[sp]), file=sys.stderr)
```
output said: guttatus kept 257000 of 275746, nasutus kept 257000 of 804072, laciniatus kept 257000 of 666896. and then do  thinning:
```
 python /project/dtataru/lac_nas_gut/AIMS/thin_positions.py 100 panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.txt > panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.thinned.txt

```
Now lets look at the numbers:
```
awk 'NR>1 {
    max=$6; winner="guttatus";
    if ($8>max) {max=$8; winner="nasutus"}
    if ($10>max) {max=$10; winner="laciniatus"}
    count[winner]++
}
END {
    for (sp in count) print sp, count[sp]
}' panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.thinned.txt
```
total of 490769 sites (771001 before thinning to one AIM per 100 bp), laciniatus 167507, nasutus 171174, guttatus 152087

REDOING for aims_v5:
wrote the script thin_byspecieswindows.py which is supposed to replace the original balancing and thinning scripts. this sets 100 kb windows and finds the species with the least amount of aims in the window and randomly thins to that number of aims, ending up with equal density per species across the genome. run it using:
```
python /project/dtataru/lac_nas_gut/AIMS/thin_by_specieswindows.py panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.species.txt 100000 \
    >  panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v4.speciesdensity.txt \
    2> thinning.stats
```
thinning stats: guttatus 260444, nasutus 260448, laciniatus 260445, total 781,337

v5 of the aims is created with  genocounts_groups_threeway_v4.py and thin_by_specieswindows.py, and specifically only keeps sites where there is >90% difference in alt allele frequency of highest frequency and lowest frequency species (like in banarjee et al. 2023). This ends up being almost equal parts each species across SNPs.

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from each hybrid individual to all parental references independently. Input files are fastq files for all samples, four lanes per sample and two reads per lane (Total 1,232 arrays).  In the TMPDIR, this creates a folder for each sample, with all reads aligned for each lane run. Output is three .sam files for each sample, corresponding to each species' reference genome. In the next step, these separate runs will be merged. This script is ```map_array_DT.sh ```.

```
#!/bin/bash
#SBATCH --output=/project/dtataru/ancestryinfer/logs/map_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/map_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1            #: Number of Nodes
#SBATCH -n 1            #: Number of Tasks per Node
#SBATCH -A loni_ferrislac
#SBATCH --cpus-per-task=12
#SBATCH --array=1-4   # Job array when n is number of unique samples, % is array throttling (do 32 at a time, 8 samples)

### LOAD MODULES ###
#For this step, bwa and needed
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES ARRAY SIZE > 1000 ###
mapfile -t R1_FILE_LIST < <(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R1_001.fastq.gz | sort)
mapfile -t R2_FILE_LIST < <(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R2_001.fastq.gz | sort)

R1_FILE_LIST_SIZE=${#R1_FILE_LIST[@]}
R2_FILE_LIST_SIZE=${#R2_FILE_LIST[@]}

# check for matching read files
for v in R1_FILE_LIST_SIZE R2_FILE_LIST_SIZE; do echo "Using $v=${!v}"; done
if [ $R1_FILE_LIST_SIZE -ne $R2_FILE_LIST_SIZE ]; then
        echo "ERROR: Encountered unmatched read files: R1_FILE_LIST_SIZE=$R1_FILE_LIST_SIZE and R2_FILE_LIST_SIZE=$R2_FILE_LIST._SIZE Exiting."
        exit 1
fi

NUMBER_FILE_PAIRS_ALREADY_DONE=1000
#NUMBER_FILE_PAIRS_ALREADY_DONE=0
NUMBER_FILE_PAIRS_TODO=$(( $R1_FILE_LIST_SIZE-$NUMBER_FILE_PAIRS_ALREADY_DONE ))
MAX_FILE_PAIRS_PER_JOB=$(( ($NUMBER_FILE_PAIRS_TODO+1000-1)/1000 ))
for v in NUMBER_FILE_PAIRS_ALREADY_DONE NUMBER_FILE_PAIRS_TODO MAX_FILE_PAIRS_PER_JOB; do
        echo "Using $v=${!v}"
done
START_FILE_INDEX=$(( $NUMBER_FILE_PAIRS_ALREADY_DONE+1+($SLURM_ARRAY_TASK_ID-1)*$MAX_FILE_PAIRS_PER_JOB))
END_FILE_INDEX=$(( $START_FILE_INDEX+$MAX_FILE_PAIRS_PER_JOB-1 ))
for v in START_FILE_INDEX END_FILE_INDEX; do echo "Using $v=${!v}"; done
for FILE_INDEX in $(seq $START_FILE_INDEX $END_FILE_INDEX); do
        echo "FILE_INDEX=$FILE_INDEX"
        if [ $FILE_INDEX -gt $R1_FILE_LIST_SIZE ]; then
                continue
        fi

R1=${R1_FILE_LIST[$((FILE_INDEX - 1))]}
R2=${R2_FILE_LIST[$((FILE_INDEX - 1))]}

SAMPLE=$(echo $R1 | cut -d "/" -f 6 | cut -d "_" -f 1-3)
HEADER=$(echo $R1 | cut -d "/" -f 6 | cut -d "_" -f 1)

genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

echo "R1=$R1"
echo "R2=$R2"
echo "SAMPLE=$SAMPLE"
echo "HEADER=$HEADER"
echo "genome1=$genome1"
echo "genome2=$genome2"
echo "genome3=$genome3"
echo "tag=$tag"

### SET TMPDIR ###
TMPDIR="/work/dtataru/TMPDIR/${HEADER}"
mkdir -p "$TMPDIR"
cd "$TMPDIR" 

### MAPPING ###
echo "Mapping ${SAMPLE} to three parental genomes"

# Read group string
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:hyblib1\tPU:LSIslowmode"

# Run bwa mem for each parental genome
bwa mem -M -R "$RG" "$genome1" "$R1" "$R2" > "${SAMPLE}.par1.sam"
bwa mem -M -R "$RG" "$genome2" "$R1" "$R2" > "${SAMPLE}.par2.sam"
bwa mem -M -R "$RG" "$genome3" "$R1" "$R2" > "${SAMPLE}.par3.sam"

echo "Mapping complete for ${SAMPLE}"

done

```

It takes ~12 hours to run. So if I'm running 32 samples at a time (max allowed by HPC, also don't queue more than 100 at a time), that would be ~30 samples a day, and it would take ~10 days to do the alignment.

## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them. Input is three .sam files for each sample (one per species reference), output is three sorted.pass.unique.bam for each sample. All of this takes place in each samples directory in TMPDIR, total array # is 308.

This is the script ```run_samtools_to_jointfiltering_DT.sh```:

```
#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --output=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac
#SBATCH --array=2  # For some reason this starts at array 2

### LOAD MODULES ###
#module load samtools/1.19
#module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
P=$(find /work/dtataru/TMPDIR/ -type d | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(echo $P | cut -d "/" -f 5 | cut -d "_" -f 1)

genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

TMPDIR="/work/dtataru/TMPDIR/${SAMPLE}"
cd "$TMPDIR"

echo "Working in TMPDIR: $TMPDIR"
echo "Sample: $SAMPLE"

### PARAMETERS ###
QUALITY=29
MAX_ALIGN=0 #was set to 2000000 before, but trying out keeping all reads
RATE=0.005
THREADS=20

### MERGE LANES ###
echo "merging lanes"

for P in 1 2 3; do
    SAM_FILES=("${TMPDIR}/${SAMPLE}"*_L00?.par${P}.sam)
	 BAM_FILES=()

   for SAM in "${SAM_FILES[@]}"; do
        BAM="${SAM%.sam}.bam"
        echo "Converting $SAM to $BAM"
        samtools fixmate -O bam "$SAM" "$BAM"
        BAM_FILES+=("$BAM")
    done

	MERGED_BAM="${TMPDIR}/${SAMPLE}.merged.par${P}.bam"
    samtools merge -f -@ ${THREADS} "$MERGED_BAM" "${BAM_FILES[@]}"
done

### SAMTOOLS SORT AND FILTER ###
echo "Starting Samtools"

for P in 1 2 3; do
    SORTED="${TMPDIR}/${SAMPLE}.par${P}.sorted.bam"
    UNIQUE="${TMPDIR}/${SAMPLE}.par${P}.sorted.unique.bam"

    samtools sort -@ 12 -o $SORTED $MERGED_BAM
    samtools index $SORTED
    samtools view -b -q $QUALITY $SORTED > $UNIQUE
    samtools index $UNIQUE
done

echo "Samtools Done"

### JOINT FILTERING ###

echo "Start Joint Filtering"

samtools view -F 4 ${SAMPLE}.par1.sorted.unique.bam | cut -f1 > p1_pass
samtools view -F 4 ${SAMPLE}.par2.sorted.unique.bam | cut -f1 > p2_pass
samtools view -F 4 ${SAMPLE}.par3.sorted.unique.bam | cut -f1 > p3_pass

comm -12 <(sort p1_pass) <(sort p2_pass) | \
comm -12 - <(sort p3_pass) > pass_all

if [[ $MAX_ALIGN -gt 0 ]]; then
    echo "Subsampling to $MAX_ALIGN alignments"
    shuf -n $MAX_ALIGN pass_all -o pass_all
fi

for P in 1 2 3; do
    samtools view -N pass_all -b ${SAMPLE}.par${P}.sorted.unique.bam > ${SAMPLE}.par${P}.sorted.pass.unique.bam
	samtools index ${SAMPLE}.par${P}.sorted.pass.unique.bam
done

echo "Job Done"

```

### 3. Variant Calling 

Joint variant calling across all samples with bcftools, then genotype calls matching parent 1 (coordinate space) alleles at ancestry informative sites are counted from a joint samtools mpileup file. Note, this script calls ```vcf_to_counts_non-colinear_DTv3.pl```. Input files are each of the ${SAMPLE}.par${P}.sorted.pass.unique.bam files, which are merged, and output files are hybrid1 vcf counts and .bed files. Initially, I looped this through all parent files but I think it is only needed for par1. I have kept the format which can easily incorporate all parents if needed (just add for P in 1 2 3). The script will fail if there are empty or corrupt bams. Somehow, everytime I have run Step 2 empty hidden files have been created in the directory, so before running remove those:

```
rm /work/dtataru/TMPDIR/.par1.sorted.pass.unique.bam
rm /work/dtataru/TMPDIR/.par2.sorted.pass.unique.bam
rm /work/dtataru/TMPDIR/.par3.sorted.pass.unique.bam
```
Also, there were three "Bad bams", GBG10, HHH46, and SHG11, which I move to the BAD_BAMS folder to be able to merge the other bams. Then run the script called ```varcall_readcount_DT.sh```:

```
#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall__%j.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

P=$(find /work/dtataru/TMPDIR/ -type d | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(echo $P | cut -d "/" -f 6 | cut -d "." -f 1)

PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
AIM_COUNTS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.txt"
WORKDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
BAMDIR="/work/dtataru/TMPDIR/"
THREADS=20

### CHECK BAM FILES FOR CORRUPTION ###
#script will fail if there are corrupt bams
#echo "Checking BAM files for corruption or emptiness..."
#CORRUPT_COUNT=0

#for f in /work/dtataru/TMPDIR/*/*.sorted.pass.unique.bam; do
#    if ! samtools quickcheck -v "$f"; then
#        echo "Corrupt or empty file detected: $f"
#        CORRUPT_COUNT=$((CORRUPT_COUNT+1))
#    fi
#done

#if [ $CORRUPT_COUNT -gt 0 ]; then
#    echo "ERROR: Found $CORRUPT_COUNT corrupt or empty BAM files. Exiting job."
#    exit 1
#else
#    echo "All BAM files passed samtools quickcheck."
#fi

### MERGE ALL BAMS FOR VARIANT CALLING ###
#echo "Merge BAM files"
#cd "$BAMDIR"

#for P in 1; do
#	BAM_FILES=($(find "$BAMDIR" -type f -name "*.par${P}.sorted.pass.unique.bam" | sort))
#	MERGED="${WORKDIR}/hybrids1merged.par${P}.pass.unique.bam"
#	SORTED="${WORKDIR}/hybrids1merged.par${P}.sorted.pass.unique.bam"

#   	samtools merge -r -c -p -@ ${THREADS} "$MERGED" "${BAM_FILES[@]}"
#	samtools sort -@ 12 -o "$SORTED" "$MERGED"
#	samtools index "$SORTED"
#done
#echo "BAM files merged"

```
In this script, the I changed the mpileup step from this to what is current, because default max-depth is 250 and that is way too low for 308 samples at 23x. It timed out after 3 days of trying to make it, though. Now I am going to run a second script for variant calling that runs an array by chromosome called ```varcall_simple.sh```:

```
#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall_%j.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -A loni_ferrislac
#SBATCH --array=1-14       
#SBATCH --mem=32G

module load bcftools
module load samtools

### ASSIGN VARIABLES ###
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
WORKDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
BAM_FILE="${WORKDIR}/hybrids1merged.par1.sorted.pass.unique.bam"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FOCAL_CHROM_LIST")
OUTVCF="hybrids1.par1.maxdepth6000.${CHR}.vcf"
PATH_SCRIPTS="/project/dtataru/ancestryinfer"

### VARIANT CALLING ###
echo "start variant calling"

cd "$WORKDIR"

#exit if vcf already exists
if [[ -s "$OUTVCF" ]]; then
    echo "Output exists, skipping: $OUTVCF"
    exit 0
fi

bcftools mpileup -Ou -d 6000 -r "$CHR" -f "$genome1" "$BAM_FILE" | bcftools call -m -Ov -o "$OUTVCF"

### GENERATE HMM INPUT FILES ###
echo "start generating hmm input"

### AIMS TO COUNTS ###
INFILE_AIMS="${OUTVCF}.v2.aims"
COUNTS="${OUTVCF}.v2_counts"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.v2.txt"
AIM_COUNTS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.v2.txt"
AIMS_MOD="${AIMS}.mod"
AIMS_BED="${AIMS}.mod.bed"
COUNTS_BED="${COUNTS}.bed"

#modify aims
#awk -v OFS='\t' '{print $1"_"$2, $1, $2, $3, $4}' "$AIMS" > "$AIMS_MOD"

#modify vcf
awk '!/^#/ {print $1"_"$2"\t"$0}' "$OUTVCF" > "$VCF_MOD"

#merge aims and vcf
perl "${PATH_SCRIPTS}/combine_FAS_scriptome_snippet.pl" "$AIMS_MOD" "$VCF_MOD" "$INFILE_AIMS"
perl "${PATH_SCRIPTS}/vcf_to_counts_non-colinear_DTv5.pl"  "$INFILE_AIMS" "$COUNTS"
perl -F'_|\t' -lane 'print join("\t", $F[0], $F[1], @F[4..$#F])' "$COUNTS" > "$COUNTS_BED"

### COUNTS TO HMM INPUT ###
perl "${PATH_SCRIPTS}/vcf_counts_to_hmm_DT.pl" "$COUNTS_BED" "$AIM_COUNTS" 0.0000000387 > "${COUNTS}.hmmsites1"

echo "Job Done"
```

Number of sites for each Chromosome

|Chrom|# of AIMS|
|-----|---------|
|  1  |  8948   |
|  2  |  14386  |
|  3  |  10409  |
|  4  |  17211  |
|  5  |  12362  |
|  6  |  15483  |
|  7  |  8235   |
|  8  |  18949  |
|  9  |  10019  |
|  10 |  14372  |
|  11 |  10866  |
|  12 |  13923  |
|  13 |  14492  |
|  14 |  24575  |

total 194,230 sites

### 4. AncestryHMM

Run main ancestryhmm program. Because ancestryinfer uses read counts instead of genotype calls, there is a step in between these in that wrapper that thin to one AIM per read for each individual, to account for non-independence. Because we are using genotype calls (-g), this is not necessary. Calling this ```run_hmm_DT.sh ```:

```
#!/bin/bash
#SBATCH --job-name=ancestryhmm
#SBATCH --output=/project/dtataru/ancestryinfer/logs/ancestryhmm_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/ancestryhmm__%j.err
#SBATCH --time=05:00:00
#SBATCH -p workq
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac
#SBATCH --array=1             # one task = one sample

### LOAD MODULES ###
module load python/3.11.5-anaconda
module load boost/1.83.0/intel-2021.5.0
module load gcc/13.2.0
module load gsl/2.7.1/intel-2021.5.0
module load bcftools
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer
export PATH=$PATH:/project/dtataru/ancestryinfer/Ancestry_HMM/src
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH

### ASSIGN VARIABLES ###
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
INPUTDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
WORKDIR="/work/dtataru/HYBRIDS/HMM_OUTPUT"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"
cd $WORKDIR

### SAMPLE INFO ###
#only create this sample list once
#echo "Creating current.samples.list file"
#bcftools query -l hybrids1.par1.maxdepth6000.Chr-01.vcf > current.samples.list

SAMPLE_ID=$SLURM_ARRAY_TASK_ID
SAMPLE_NAME=$(sed -n "${SAMPLE_ID}p" "${INPUTDIR}/current.samples.list" | awk '{print $1}')
echo "Processing SAMPLE: $SAMPLE_NAME    SAMPLE_ID=$SAMPLE_ID"

### RUN HMM ###

### LOOP OVER CHROMS ###
while read CHROM; do

    echo "  → Processing chromosome $CHROM"

    INPUT="${INPUTDIR}/hybrids1.par1.maxdepth6000.${CHROM}.vcf_counts.hmmsites1"

    ### DETERMINE SAMPLE COLUMNS ###
    A_col=$((10 + (SAMPLE_ID-1)*2))
    a_col=$((11 + (SAMPLE_ID-1)*2))

    echo "     Sample columns: A=$A_col   a=$a_col"

    OUTFILE="${WORKDIR}/${SAMPLE_NAME}.${CHROM}.counts.hmmsites1"

    ### Extract first 9 columns + sample's two genotype columns ###
    awk -v A=$A_col -v a=$a_col -v OFS="\t" \
        '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$A,$a}' \
        "$INPUT" > "$OUTFILE"

    ### Write sample-specific sample list file ###
    echo -e "${SAMPLE_NAME}\t2" > ${SAMPLE_NAME}.${CHROM}.samples

    ### RUN ANCESTRY_HMM ###
    ancestry_hmm \
        -i "$OUTFILE" \
        -s "${SAMPLE_NAME}.${CHROM}.samples" \
        -a 3 0.33 0.33 0.34 \
        -p 0 -10000 0.33 \
        -p 1 -1000 0.33 \
        -p 2 -500 0.34 \
        -e 0.05 -g

    echo "Done with $SAMPLE_NAME on $CHROM"

done < "$FOCAL_CHROM_LIST"

echo "ALL COMPLETE for sample $SAMPLE_NAME"

```

## 5. PARSE TSVs FOR PLOTTING
Parsing ancestry from the output tsvs for downstream analysis and visualization. Output of this can be analyzed in something similar to  Banarjee et al. 2023's```Tlalica_three-way_hybrids/local_ancestry_calling/local_ancestry_plots.R```. This is called '''5_ProcessHMM_Output.sh'''

```
#!/bin/bash
#SBATCH --job-name=ProcessHMM_Output
#SBATCH --output=/project/dtataru/ancestryinfer/logs/ProcessHMM_Output_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/ProcessHMM_Output_%j.err
#SBATCH --time=03:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A loni_ferrislac                      

# Load Modules
module load R/3.2.4

# Set variables
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
INPUT_DIR="/work/dtataru/HYBRIDS/HMM_OUT"
WORKING_DIR="work/dtataru/HYBRIDS/HMM_POSTPROCESS"
POSTERIOR_THRESH=0.8
INTERVALS_RSCRIPT="${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"      
TAG=$(paste -sd "_" "$FOCAL_CHROM_LIST")

cd $WORKING_DIR

### Create Input Files ###

current_list="current.posterior.samples.list_${TAG}"
read_list="current.samples.read.list_${TAG}"

# loop through all .posterior files
for p in *.posterior; do
    sample="${p%.posterior}"
    echo "$sample" >> $current_list

    for hmmsites in ${INPUT_DIR}/${sample}.Chr-*.counts.hmmsites1; do
        if [[ -f "$hmmsites" ]]; then
            echo "$hmmsites" >> "$read_list"
        else
            echo "WARNING: no hmmsites1 files found for sample $sample" >&2
        fi
    done
done

#remove duplicates
sort -u "$current_list" -o "$current_list"
sort -u "$read_list" -o "$read_list"

### CONVERT TO TSV ###

perl ${PATH_SCRIPTS}/convert_rchmm_to_ancestry_tsv_3way_v2_DT.pl current.posterior.samples.list current.samples.read.list 1 ${FOCAL_CHROM_LIST}

echo "Transposing ancestry probabilities"
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv

### parse ancestry in transposed tsvs ###
echo "parse transposed tsv start"
POSTERIOR_THRESH=0.8
perl ${PATH_SCRIPTS}/parse_3way_tsv_to_genotypes_file.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv $POSTERIOR_THRESH > genotypes.txt

perl ${PATH_SCRIPTS}/parse_3way_tsv_ancestry.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv 0.8 > ancestry_proportions.txt

echo "all files parsed"

### identify intervals ###
echo "run intervals R script"
Rscript ${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R genotypes.txt ${PATH_SCRIPTS}
echo "identified intervals"
```
