#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall_%j.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -A loni_ferrislac
#SBATCH --array=1-14

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
#echo "start variant calling"

cd "$WORKDIR"

#exit if vcf already exists
#if [[ -s "$OUTVCF" ]]; then
#    echo "Output exists, skipping: $OUTVCF"
#    exit 0
#fi

#bcftools mpileup -Ou -d 6000 -r "$CHR" -f "$genome1" "$BAM_FILE" | bcftools call -m -Ov -o "$OUTVCF"

### GENERATE HMM INPUT FILES ###
echo "start generating hmm input"

### AIMS TO COUNTS ###
INFILE_AIMS="${OUTVCF}.v3.aims"
COUNTS="${OUTVCF}.v3_counts"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.v3.txt"
AIM_COUNTS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.v3.txt"
AIMS_MOD="${AIMS}.mod"
COUNTS_BED="${COUNTS}.bed"
VCF_MOD="${OUTVCF}.mod"

#modify aims
#awk -v OFS='\t' '{print $1"_"$2, $1, $2, $3, $4}' "$AIMS" > "$AIMS_MOD"

#modify vcf
#awk '!/^#/ {print $1"_"$2"\t"$0}' "$OUTVCF" > "$VCF_MOD"

#merge aims and vcf
perl "${PATH_SCRIPTS}/combine_FAS_scriptome_snippet.pl" "$AIMS_MOD" "$VCF_MOD" "$INFILE_AIMS"
perl "${PATH_SCRIPTS}/vcf_to_counts_non-colinear_DTv5.pl"  "$INFILE_AIMS" "$COUNTS"
perl -F'_|\t' -lane 'print join("\t", $F[0], $F[1], @F[4..$#F])' "$COUNTS" > "$COUNTS_BED"

### COUNTS TO HMM INPUT ###
#changed recombination rate to 0.0000000387 from 0.00000002 following Farnitano et al. 2025
perl "${PATH_SCRIPTS}/vcf_counts_to_hmm_DT.pl" "$COUNTS_BED" "$AIM_COUNTS" 0.0000000387 > "${COUNTS}.hmmsites1"

echo "Job Done"
