#!/bin/bash
#
# Create batch jobs from raw data using pipeline
# Original workflow by Danny Arends
# Adapted by Sean Anderson (anderss@wfu.edu)
# Feb 2024
#

################################################################################
################################################################################

SOFTWARE="${1}"
PIPELINE="${2}"
DATA_RAW="${3}"
DATA_STR="${4}"
DATA_PRI="${5}"
DATA_VCF="${6}"
SLURM_PARTITION="${7}"
SLURM_NTASKS="${8}"
SLURM_MEM="${9}"
SLURM_TIME="${10}"
OUTPUT="${11}"
BATCH_SCRIPT="${12}"

################################################################################
################################################################################

declare -a CASES=($(find $DATA_RAW -maxdepth 1 -mindepth 1 -type d | sort -V))
# CASES=(/deac/inf/adminGrp/anderss/scratch/zhangGrp/DATA/rna_seq_7.26.2021/raw_data/B1)

cat <<-EOD > $BATCH_SCRIPT
#!/bin/bash
for DIR in ${CASES[@]}; do
CASE="\$( basename \$DIR )"
sbatch --job-name="\$CASE" \\
       --partition=$SLURM_PARTITION \\
       --nodes=1 \\
       --ntasks-per-node=$SLURM_NTASKS \\
       --mem=$SLURM_MEM \\
       --time=$SLURM_TIME \\
       --mail-type=fail \\
       --mail-user="${USER}@wfu.edu" \\
       --output=slurm-%x-%j.o \\
<< EOF
#!/bin/bash
bash $PIPELINE $SOFTWARE \$CASE \$DIR /scratch/\\\${SLURM_JOB_ID}/output $DATA_STR $DATA_PRI $DATA_VCF $SLURM_NTASKS
mkdir -p $OUTPUT
mv /scratch/\\\${SLURM_JOB_ID}/output/* $OUTPUT
EOF
done
EOD

printf "\n\nYour final batch submissions have been generated; submit your jobs by running:\n\nbash %s\n" $(realpath $BATCH_SCRIPT)
