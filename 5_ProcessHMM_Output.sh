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
INPUT_DIR="/work/dtataru/HYBRIDS/HMM_OUTPUT"
WORKING_DIR="/work/dtataru/HYBRIDS/HMM_POSTPROCESS"
POSTERIOR_THRESH=0.8
INTERVALS_RSCRIPT="${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"
TAG=$(paste -sd "_" "$FOCAL_CHROM_LIST")

cd $WORKING_DIR

### Create Input Files ###

current_list="current.posterior.samples.list_${TAG}"

# loop through all .posterior files
for p in *.posterior; do
    sample="${p%.posterior}"
    echo "$sample" >> $current_list

    for hmmsites in ${INPUT_DIR}/${sample}.Chr-*.v2.counts.hmmsites1; do
        if [[ -f "$hmmsites" ]]; then
            echo "$hmmsites" >> "$read_list"
        else
            echo "WARNING: no hmmsites1 files found for sample $sample" >&2
        fi
    done
done

#remove duplicates
sort -u "$current_list" -o "$current_list"

### CONVERT TO TSV ###

perl ${PATH_SCRIPTS}/convert_rchmm_to_ancestry_tsv_3way_v1.pl current.posterior.samples.list current.posterior.samples.list 1 ${FOCAL_CHROM_LIST}

echo "Transposing ancestry probabilities"
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv

### parse ancestry in transposed tsvs ###
#echo "parse transposed tsv start"
POSTERIOR_THRESH=0.8
perl ${PATH_SCRIPTS}/parse_3way_tsv_to_genotypes_file.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv $POSTERIOR_THRESH > genotypes.txt

perl ${PATH_SCRIPTS}/parse_3way_tsv_ancestry.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv $POSTERIOR_THRESH > ancestry_proportions_post0.3.txt

perl ${PATH_SCRIPTS}/parse_3way_tsv_ancestry_BYCHROM.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv $POSTERIOR_THRESH

echo "all files parsed"

### identify intervals ###
#echo "run intervals R script"
#Rscript ${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R genotypes.txt ${PATH_SCRIPTS}
#echo "identified intervals"
