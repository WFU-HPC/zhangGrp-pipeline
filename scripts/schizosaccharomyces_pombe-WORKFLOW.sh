#!/bin/bash
#
# Complete workflow for schizosaccharomyces_pombe
# Original workflow by Danny Arends
# Adapted by Sean Anderson (anderss@wfu.edu)
# Feb 2024
#

################################################################################
################################################################################

export RESEARCH="/deac/bio/zhangGrp"

################################################################################
# Step 0: Install your software (ONLY ONCE!)
################################################################################

bash ${RESEARCH}/scripts/0_installSoftware.sh ${RESEARCH}/software

################################################################################
# Step 1: Download the genome for schizosaccharomyces_pombe
################################################################################

bash ${RESEARCH}/scripts/1_buildGenome.sh \
        ${RESEARCH}/software \
        ftp.ensemblgenomes.org/pub/fungi/release-57/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.chromosome \
        ftp.ensemblgenomes.org/pub/fungi/release-57/gtf/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.57.gtf.gz \
        ftp.ensemblgenomes.org/pub/fungi/release-57/variation/vcf/schizosaccharomyces_pombe/schizosaccharomyces_pombe.vcf.gz \
        ${RESEARCH}/DATA/genome

################################################################################
# Step 2: Prepare the genome for schizosaccharomyces_pombe
################################################################################

bash ${RESEARCH}/scripts/2_prepGenome.sh \
        ${RESEARCH}/software \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.dna.primary_assembly.fa.gz \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.57.gtf.gz \
        ${RESEARCH}/DATA/genome/schizosaccharomyces_pombe.vcf.gz \
        ${RESEARCH}/DATA/genome/STAR

################################################################################
# Step 3: Create batch jobs for aligning schizosaccharomyces_pombe
################################################################################

bash ${RESEARCH}/scripts/3a_batchPipeline.sh \
        ${RESEARCH}/software     \
        ${RESEARCH}/scripts/3_pipeline.sh \
        ${RESEARCH}/DATA/rna_seq_7.26.2021/raw_data \
        ${RESEARCH}/DATA/genome/STAR \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.dna.primary_assembly.fa.gz \
        ${RESEARCH}/DATA/genome/schizosaccharomyces_pombe.vcf.gz \
        "small" \
        "16" \
        "20G" \
        "00-06:00:00" \
        ${RESEARCH}/DATA/schizosaccharomyces_pombe \
        ./schizosaccharomyces_pombe-BATCH.sh
