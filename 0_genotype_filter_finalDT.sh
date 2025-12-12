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

#this script is adapted from the following manuscript: Fluctuating reproductive isolation and stable ancestry structure in a fine-scaled mosaic of #hybridizing Mimulus monkeyflowers" by Matthew Farnitano, Keith Karoly, and Andrea Sweigart, and github: https://github.com/mfarnitano/CAC_popgen/#blob/main/reference_panels/genotype_filter.sh

PROJECT=lacnasgut
BATCH_NAME=lacnasgut_panel

SCRIPTS_DIR=/project/dtataru/lac_nas_gut/AIMS
LOG_DIR=/project/dtataru/lac_nas_gut/AIMS/logs
WORKING_DIR=/project/dtataru/lac_nas_gut/AIMS
VCF_DIR=/project/dtataru/lac_nas_gut/4_ref/3_Genotyped_GVCFs
TMPDIR=/work/dtataru/AIMS

GENOME=/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa
REPEATMASK=/project/dtataru/lac_nas_gut/AIMS/MguttatusTOL_551_v5.0.repeatmasked_assembly_v5.0.gff3
#FASTQ_TABLE=${WORKING_DIR}/lacnasgut_panel_fastq_table_bigger.txt #two columns, directory and prefix

NTHREADS=32
MEM=120

#inputs
RAW_VCF=${VCF_DIR}/lacnasgut_jointgeno.vcf.gz
GROUP=panel15
REFCODE=TOL551
PREFIX=${WORKING_DIR}/VCFs/${GROUP}.${REFCODE}

###SETUP DIRS
cd $WORKING_DIR

###MODULES
#module load gatk/4.5.0.0

###create chr_list
#printf "\n...creating chromosome list\n" | tee >(cat >&2)
#if [ ! -f ${WORKING_DIR}/chr_positions.list ]; then
#	head -n14 ${GENOME}.fai | awk '{print $1 ":1-"$2}' > ${WORKING_DIR}/chr_positions.list
#fi

#printf "\n...extracting biallelic SNPs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
#	-select-type SNP --restrict-alleles-to BIALLELIC -O ${PREFIX}.SNPs.vcf.gz

#printf "\n...extracting invariant sites\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
#	-select-type NO_VARIATION -O ${PREFIX}.INVTs.vcf.gz

#printf "\n...filtering SNPs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" VariantFiltration -V ${PREFIX}.SNPs.vcf.gz -O ${PREFIX}.SNPs.f.vcf.gz \
#	-filter "QD < 2.0" --filter-name "QD2" \
#	-filter "QUAL < 40.0" --filter-name "QUAL40" \
#	-filter "SOR > 3.0" --filter-name "SOR4" \
#	-filter "FS > 60.0" --filter-name "FS60" \
#	-filter "MQ < 40.0" --filter-name "MQ40" \
#	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#	-filter "ReadPosRankSum < -12.5" --filter-name "ReadPosRankSum-12.5" \
#	-filter "ReadPosRankSum > 12.5" --filter-name "ReadPosRankSum12.5" \
#	--verbosity ERROR

#printf "\n...filtering INVTs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}g" VariantFiltration -V ${PREFIX}.INVTs.vcf.gz -O ${PREFIX}.INVTs.f.vcf.gz \
#	-filter "QD < 2.0" --filter-name "QD2" \
#	-filter "SOR > 3.0" --filter-name "SOR4" \
#	-filter "MQ < 40.0" --filter-name "MQ40" \
#	--verbosity ERROR

#printf "\n...sorting SNPs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SortVcf \
#	-I ${PREFIX}.SNPs.f.vcf.gz -use_jdk_inflater --TMP_DIR ${TMPDIR} \
#	-SD ${GENOME}.dict -O ${PREFIX}.SNPs.fs.vcf.gz

#printf "\n...sorting INVTs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SortVcf \
#	-I ${PREFIX}.INVTs.f.vcf.gz -SD ${GENOME}.dict -O ${PREFIX}.INVTs.fs.vcf.gz

#printf "\n...selecting passing SNPs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants \
#	-V ${PREFIX}.SNPs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.SNPs.fsp.vcf.gz

#printf "\n...selecting passing INVTs\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants \
#	-V ${PREFIX}.INVTs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.INVTs.fsp.vcf.gz

#module purge
#module load vcftools/0.1.16
#module load htslib/1.21

#printf "\n...prepare bed file of repeat-masked regions\n" | tee >(cat >&2)
#printf '#chrom\tchromStart\tchromEnd\n' > ${REPEATMASK}.bed
#cut -f1,4,5 ${REPEATMASK} >> ${REPEATMASK}.bed

#printf "\n...filter individual SNP genotypes by depth and GQ, and exclude repeat-masked regions\n" | tee >(cat >&2)
#vcftools --gzvcf ${PREFIX}.SNPs.fsp.vcf.gz -c --minGQ 15 --minDP 6 --maxDP 100 \
#	--exclude-bed ${REPEATMASK}.bed \
#	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.SNPs.fspi.rm.vcf.gz

#printf "\n...filter individual INVT genotypes by depth\n" | tee >(cat >&2)
#vcftools --gzvcf ${PREFIX}.INVTs.fsp.vcf.gz -c --minDP 6 --maxDP 100 \
#	--exclude-bed ${REPEATMASK}.bed \
#	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.INVTs.fspi.rm.vcf.gz

#tabix -p vcf ${PREFIX}.SNPs.fspi.rm.vcf.gz
#tabix -p vcf ${PREFIX}.INVTs.fspi.rm.vcf.gz

#module purge
#module load gatk/4.5.0.0

#printf "\n...merge SNP and INVT sites\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" MergeVcfs \
#	-I ${PREFIX}.SNPs.fspi.rm.vcf.gz -I ${PREFIX}.INVTs.fspi.rm.vcf.gz \
#	-O ${PREFIX}.merged.fspi.rm.vcf.gz

#printf "\n...sort merged VCF\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SortVcf \
#	-I ${PREFIX}.merged.fspi.rm.vcf.gz \
#	-SD ${GENOME}.dict -O ${PREFIX}.merged.fspi.rm.vcf.gz

#printf "\n...filter SNPs by mincalled\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.SNPs.fspi.rm.vcf.gz \
#	--max-nocall-number 7 --exclude-filtered -O ${PREFIX}.SNPs.fspi.rm.31called.vcf.gz

#printf "\n...filter merged VCF by mincalled\n" | tee >(cat >&2)
#gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.merged.fspi.rm.vcf.gz \
#	--max-nocall-number 7 --exclude-filtered -O ${PREFIX}.merged.fspi.rm.31called.vcf.gz

#printf "\n...completing all filtering steps\n" | tee >(cat >&2)

#printf "\n...making genotype tables\n" | tee >(cat >&2)
#module purge
#module load bcftools/1.18

#printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.SNPs.fspi.rm.table
#printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.merged.fspi.rm.table

#bcftools query -l ${PREFIX}.SNPs.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.SNPs.fspi.rm.table
#bcftools query -l ${PREFIX}.merged.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.merged.fspi.rm.table

#printf '\n' >> ${PREFIX}.SNPs.fspi.rm.table
#printf '\n' >> ${PREFIX}.merged.fspi.rm.table

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.SNPs.fspi.rm.vcf.gz >> ${PREFIX}.SNPs.fspi.rm.table
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.merged.fspi.rm.vcf.gz >> ${PREFIX}.merged.fspi.rm.table


##filtering, edits starting here made by Diana Tataru, December 2025

module load python3

#python ${SCRIPTS_DIR}/genocounts_groups_threeway_v4.py ${WORKING_DIR}/lacnasgut_pops.txt ${PREFIX}.SNPs.fspi.rm.table > ${PREFIX}.SNPs.fspi.rm.counts.v5.txt

#awk -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,$5,$6,$9,$10,$13,$14}' \
#    ${PREFIX}.SNPs.fspi.rm.counts.v5.txt \
#    > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.txt


#printf "\n...make a column with species classification\n" | tee >(cat >&2)
#awk 'NR>1 {
#    g=$6; n=$8; l=$10;

#    max=g; species="guttatus";
#    if (n>max) {max=n; species="nasutus"}
#    if (l>max) {species="laciniatus"}

#    print $0 "\t" species
#}' panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v5.txt \
#    > panel15.TOL551.SNPs.fspi.rm.AIMs_counts.v5.species.txt

### OLD BALANCE AND THIN###
#printf "\n...balance number of aims by species\n" | tee >(cat >&2)
#python  ${SCRIPTS_DIR}/balance_species.py ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v4.species.txt ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.txt 257000

#printf "\n...thin to one AIM per 100 bp\n" | tee >(cat >&2)
#python ${SCRIPTS_DIR}/thin_positions.py 100 ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.txt ${PREFIX}.SNPS.thinned_coords.txt

#awk 'NR==FNR {coords[$1":"$2]; next} ($1":"$2) in coords' \
#    ${PREFIX}.SNPS.thinned_coords.txt \
#    ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.txt \
#    > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v4.speciesbalanced.thinned.txt

### NEW BALANCE BY DENSITY IN 100 KB WINDOWS
#printf "\n...balance number of aims by species density in 100kb windows\n" | tee >(cat >&2)
#python ${SCRIPTS_DIR}/thin_by_specieswindows.py ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.species.txt 100000 \
#    >  ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensity.txt \
#    2> thinning.stats

#sort -k1,1V -k2,2n ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensity.txt > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt

printf "\n...edited to remove underscore in chromosome\n" | tee >(cat >&2)
awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $3, $4}' \
    ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.v3.txt

awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $5, $6, $7, $8, $9, $10}' \
    ${PREFIX}.SNPs.fspi.rm.AIMs_counts.v5.speciesdensitysorted.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.v3.txt



