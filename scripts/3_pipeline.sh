#!/bin/bash
#
# Pipeline for aligning SRA reads to the desired genome
# Original workflow by Danny Arends
# Adapted by Sean Anderson (anderss@wfu.edu)
# Feb 2024
#

################################################################################
################################################################################

# Read software path from command line
SOFTWARE="$1"

# set up environment (improve in the future)
module load compilers/gcc/10.2.0 python/3.8.13 R/4.2.1 utils/git/2.36.1 && unset JAVA_HOME
export PATH="${SOFTWARE}/bin:${SOFTWARE}/java/bin:$PATH"

################################################################################
################################################################################

ulimit -n 24000

# command line inputs
BASE="$2"
INPUT="$3"
OUTDIR="$4"
GENOMESTAR="$5"
ASSEMBLY="$6"
SNPS="$7"
CORES="$8"
OUTPUT="${OUTDIR}/${BASE}.aln"

# Make both the input and output directorires
# mkdir -p $INPUT
mkdir -p $OUTPUT

# Change into the input directory
# cd $INPUT

# STEP 0 - SRA Download and Compress
# fasterq-dump --threads 8 -p --split-files $BASE
# bgzip --threads 8 "${BASE}_1.fastq"
# bgzip --threads 8 "${BASE}_2.fastq"

# STEP 1 - READ Trimming
TRIM_FILES=(
            "${INPUT}/${BASE//_}_1.fq.gz"
            "${INPUT}/${BASE//_}_2.fq.gz"
            "${OUTPUT}/${BASE//_}_1.P.fq.gz"
            "${OUTPUT}/${BASE//_}_1.U.fq.gz"
            "${OUTPUT}/${BASE//_}_2.P.fq.gz"
            "${OUTPUT}/${BASE//_}_2.U.fq.gz"
          )

java -Xmx16g -jar ${SOFTWARE}/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE \
     -threads $CORES \
     ${TRIM_FILES[@]} \
     ILLUMINACLIP:${SOFTWARE}/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# STEP 1.1 - UNZIP for STAR
gunzip < ${TRIM_FILES[2]} > ${TRIM_FILES[2]%.gz}
gunzip < ${TRIM_FILES[4]} > ${TRIM_FILES[4]%.gz}

# STEP 2 - Alignment using STAR
STAR --runMode alignReads \
     --runThreadN $CORES \
     --limitBAMsortRAM 4000000000 \
     --readFilesIn ${TRIM_FILES[2]%.gz} ${TRIM_FILES[4]%.gz} \
     --genomeDir ${GENOMESTAR} \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "${OUTPUT}/${BASE//_}"

# STEP 2.1 - Create a samtools index
STAR_BAM=${OUTPUT}/${BASE//_}Aligned.sortedByCoord.out.bam
samtools index $STAR_BAM
# STEP 2.2 - Create mapping and coverage statistics
samtools flagstats $STAR_BAM
samtools coverage $STAR_BAM

#STEP 3 - Remove duplicate reads using picard tools
P_BAM="${OUTPUT}/${BASE//_}Aligned.sortedByCoord.RD.out.bam"
METRICS="${OUTPUT}/${BASE//_}_metrics.txt"
java -Xmx16g -jar ${SOFTWARE}/picard/picard.jar \
     MarkDuplicates \
     --REMOVE_DUPLICATES true \
     -I $STAR_BAM \
     -O $P_BAM \
     -M $METRICS

# STEP 3.1 - Create a samtools index
samtools index $P_BAM
# STEP 3.2 - Create mapping and coverage statistics
samtools flagstats $P_BAM
samtools coverage $P_BAM

# STEP 4 - Add read group (1) and sample run, library, and name
RG_BAM=${OUTPUT}/${BASE//_}Aligned.sortedByCoord.RD.RG.out.bam
java -Xmx16g -jar ${SOFTWARE}/picard/picard.jar \
     AddOrReplaceReadGroups \
     -I $P_BAM \
     -O $RG_BAM \
     -PL ILLUMINA \
     -PU run \
     -LB ${BASE//_} \
     -SM ${BASE//_}

# STEP 4.1 - Create a samtools index
samtools index $RG_BAM

# STEP 5 - GATK prep
# STEP 5.1 - GATK BaseRecalibrator
GATK_COV1="${OUTPUT}/${BASE//_}_cov1.txt"
java -Xmx16g -jar ${SOFTWARE}/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
     BaseRecalibrator \
     -R $ASSEMBLY \
     --known-sites $SNPS \
     -I $RG_BAM \
     -O $GATK_COV1

# STEP 5.2 - GATK ApplyBQSR
RECAL_BAM="${OUTPUT}/${BASE//_}Aligned.sortedByCoord.RD.RG.RC.out.bam"
java -Xmx16g -jar ${SOFTWARE}/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
     ApplyBQSR \
     -R $ASSEMBLY \
     -bqsr $GATK_COV1 \
     -I $RG_BAM \
     -O $RECAL_BAM

# STEP 5.3 - GATK BaseRecalibrator
GATK_COV2="${OUTPUT}/${BASE//_}_cov2.txt"
java -Xmx16g -jar ${SOFTWARE}/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
     BaseRecalibrator \
     -R $ASSEMBLY \
     --known-sites $SNPS \
     -I $RECAL_BAM \
     -O $GATK_COV2

# STEP 5.4 - GATK AnalyzeCovariates
java -Xmx16g -jar ${SOFTWARE}/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
     AnalyzeCovariates \
     -before $GATK_COV1 \
     -after $GATK_COV2 \
     -plots "${OUTPUT}/${BASE//_}AnalyzeCovariates.pdf"

# STEP 6 - Index the recalibrated bam files
samtools index $RECAL_BAM

# STEP 6.1 - Create mapping and coverage statistics
samtools flagstats $RECAL_BAM
samtools coverage $RECAL_BAM

# STEP 7 - Convert bam files into bw files
bamCoverage -b ${OUTPUT}/${BASE//_}Aligned.sortedByCoord.out.bam -o ${OUTPUT}/${BASE//_}coverage.bw

# remove STAR tmp directory
rm -rf "${OUTPUT}/${BASE//_}_STARtmp"

# remove uncompressed files that STAR used
rm -f "${TRIM_FILES[2]%.gz}"
rm -f "${TRIM_FILES[4]%.gz}"
rm -f "${TRIM_FILES[2]}"
rm -f "${TRIM_FILES[3]}"
rm -f "${TRIM_FILES[4]}"
rm -f "${TRIM_FILES[5]}"
