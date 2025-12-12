#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --output=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -A loni_ferrislac
#SBATCH --array=2  # Starts at array 2, total is  309

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
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

