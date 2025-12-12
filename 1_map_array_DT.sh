#!/bin/bash
#SBATCH --output=/project/dtataru/ancestryinfer/logs/map_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/map_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1            #: Number of Nodes
#SBATCH -n 1            #: Number of Tasks per Node
#SBATCH -A loni_ferrislac
#SBATCH --cpus-per-task=12
#SBATCH --array=201-234%32   # Job array when n is number of unique samples, % is array throttling (do 32 at a time, 8 samples)

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
