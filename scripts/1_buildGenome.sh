#!/bin/bash
#
# Download genome files
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

function getgen () {
    local url="$1"
    local base="$2"
    lftp -e "cls --quiet -1 ${base}*; exit" "${url}" 2>/dev/null
}

function roman () {
    local number="$1"
    roman_numerals=$(echo $number | tr a-z A-Z)

    # Test that it is valid
    [[ "${roman_numerals//[IVXLCDM]/}" == "" ]] || {
        # echo "Roman numerals $1 contains invalid characters"
        echo "$1"
        return 1
    }

    [[ $roman_numerals =~ ^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$ ]] || {
        echo "$roman_numerals isn't a standard number"
        return 2
    }

    # We want to replace all tokens to eventually have
    # all I's, remove new lines and count the characters.
    number=$(
        echo ${roman_numerals} |
        sed -e 's/CM/DCD/g'  -e 's/M/DD/g'  -e 's/CD/CCCC/g' \
            -e 's/D/CCCCC/g'  -e 's/XC/LXL/g'  -e 's/C/LL/g' \
            -e 's/XL/XXXX/g'  -e 's/L/XXXXX/g'  -e 's/IX/VIV/g' \
            -e 's/X/VV/g'  -e 's/IV/IIII/g'  -e 's/V/IIIII/g' |
        tr -d '\n' |
        wc -m
    )

    printf "%04d" ${number}
}

################################################################################
################################################################################

# read options from command line
URL1_ORIG="$2"
URL2_ORIG="$3"
URL3_ORIG="$4"
GENOME="$5"

# do some parsing
URL1=$(dirname $URL1_ORIG)
URL2=$(dirname $URL2_ORIG)
URL3=$(dirname $URL3_ORIG)
BASE1=$(basename $URL1_ORIG)
BASE2=$(basename $URL2_ORIG)
BASE3=$(basename $URL3_ORIG)

# Make Directory
mkdir -p ${GENOME}

# Get listing of files
FILES=($(getgen "$URL1" "$BASE1"))

# Download
printf "Downloading Genome files:\n\n"
for FILE in ${FILES[@]}; do
    NEWNUM=$(roman $(echo ${FILE%.fa.gz} | awk -F"${BASE1}" '{print $2}' | sed -e 's/\.//g'))
    printf "Downloading file ${URL1}/${FILE} ...\n"
    wget --quiet --no-check-certificate "${URL1}/${FILE}" -O "${GENOME}/${BASE1}.${NEWNUM}.fa.gz"
done

# ASK USER IF OBTAINED FILES ARE CORRECT, and MAYBE CONFIRM LOCATION?

# New assembly file name
ASSEMBLY="${GENOME}/${BASE1%.chromosome}.primary_assembly.fa"

# If assembly file exists, remove it!
if [ -f "$ASSEMBLY" ]; then
    printf "\nRemoving previous assembly file at\n${ASSEMBLY}\nand generating a new one.\n"
    rm $ASSEMBLY
fi

# Extract and merge into a combined fast file
printf "Concatenating genome files into new fast file.\n\n"
for NEWFILE in "${GENOME}/${BASE1}.*.fa.gz"; do
    zcat $NEWFILE >> $ASSEMBLY
done

# Compress the fasta file using bgzip (keep original)
bgzip -fk $ASSEMBLY

# DOES THE ORIGINAL FILE NEED TO BE DELETED????

# Delete the chromosomes
for NEWFILE in "${GENOME}/${BASE1}.*.fa.gz"; do
    rm $NEWFILE
done

# Get the reference transcriptome and GTF for STAR
printf "Downloading file ${URL2}/${BASE2}...\n"
wget --quiet --no-check-certificate ${URL2}/${BASE2} -O ${GENOME}/${BASE2}
gunzip < ${GENOME}/${BASE2} > ${GENOME}/${BASE2%.gz}
printf "Downloading file ${URL3}/${BASE3}...\n"
wget --quiet --no-check-certificate ${URL3}/${BASE3} -O ${GENOME}/${BASE3}


printf "\n\nFinal resulting files are:\n"
ls ${GENOME}/${BASE1%.chromosome}* ${GENOME}/${BASE2}* ${GENOME}/${BASE3}*
