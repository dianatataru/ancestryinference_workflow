#!/bin/bash
#SBATCH --job-name=ancestryhmm
#SBATCH --output=/project/dtataru/ancestryinfer/logs/ancestryhmm_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/ancestryhmm__%j.err
#SBATCH --time=05:00:00
#SBATCH -p workq
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac
#SBATCH --array=251-305%40             # one task = one sample, 305 total samples

### LOAD MODULES ###
module load python/3.11.5-anaconda
module load boost/1.83.0/intel-2021.5.0
module load gcc/13.2.0
module load gsl/2.7.1/intel-2021.5.0
module load bcftools
module load r/4.3.2/gcc-9.3.0 
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
        -p 0 -1000 0.33 \
        -p 1 -1000 0.33 \
        -p 2 -1000 0.34 \
        -e 0.05 -g

    echo "Done with $SAMPLE_NAME on $CHROM"


 #rename so it doesn't get overwritten
    mv "${SAMPLE_NAME}.posterior" "${SAMPLE_NAME}.${CHROM}.posterior"

done < "$FOCAL_CHROM_LIST"

### MERGE CHROMOSOMES IN ORDER ###
echo "Merging chromosome posterior files into final ${SAMPLE_NAME}.posterior"

while read CHROM; do
    cat "${SAMPLE_NAME}.${CHROM}.posterior" >> "${SAMPLE_NAME}.posterior"
done < "$FOCAL_CHROM_LIST"

cp "${SAMPLE_NAME}.posterior" "/work/dtataru/HYBRIDS/HMM_POSTPROCESS/${SAMPLE_NAME}.posterior"


echo "ALL COMPLETE for sample $SAMPLE_NAME"
