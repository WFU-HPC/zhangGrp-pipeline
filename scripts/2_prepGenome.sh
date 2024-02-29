#!/bin/bash
#
# Prepare genome files
# Original workflow by Danny Arends
# Adapted by Sean Anderson (anderss@wfu.edu)
# Feb 2024
#

################################################################################
################################################################################

# Read software path from command line
SOFTWARE="$1"

# set up environment (improve in the future)
module load compilers/gcc/10.2.0 utils/git/2.36.1 && unset JAVA_HOME
export PATH="${SOFTWARE}/bin:${SOFTWARE}/java/bin:$PATH"

################################################################################
################################################################################

# Command line options
FILE_ASM="$2"
FILE_GTF="$3"
FILE_VCF="$4"
STAR="$5"
CORES="8"

# Index the genome using samtools
samtools faidx $FILE_ASM

# Generate genome/transcriptome index using STAR
STAR --runThreadN $CORES \
     --runMode genomeGenerate \
     --genomeDir $STAR \
     --genomeSAindexNbases 10 \
     --sjdbGTFfile "${FILE_GTF%.gz}" \
     --genomeFastaFiles "${FILE_ASM%.gz}"

# Get the reference SNPs and index using tabix
gunzip < $FILE_VCF > ${FILE_VCF%.gz}
bgzip < ${FILE_VCF%.gz} > $FILE_VCF
tabix -p vcf $FILE_VCF # CONSIDER 0 INDEX?????????????

# Index the genome using picard
java -Xmx4g -jar ${SOFTWARE}/picard/picard.jar \
     CreateSequenceDictionary \
     -R $FILE_ASM
