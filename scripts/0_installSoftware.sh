#!/bin/bash
#
# Workflow software installation
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

function sw_openjdk {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    wget --no-check-certificate https://download.java.net/openjdk/jdk17/ri/openjdk-17+35_linux-x64_bin.tar.gz                                               -O /tmp/openjdk.tar.gz
    # wget --no-check-certificate https://download.java.net/java/GA/jdk20.0.2/6e380f22cbe7469fa75fb448bd903d8e/9/GPL/openjdk-20.0.2_linux-x64_bin.tar.gz      -O /tmp/openjdk.tar.gz
    # wget --no-check-certificate https://download.java.net/java/GA/jdk21.0.2/f2283984656d49d69e91c558476027ac/13/GPL/openjdk-21.0.2_linux-x64_bin.tar.gz     -O /tmp/openjdk.tar.gz
    
    # untar and prepare with build directory outside of source tree
    tar -xf /tmp/openjdk.tar.gz -C /tmp
    mv /tmp/jdk-17 $target

    # remove spurious leftovers
    rm -rf /tmp/openjdk.tar.gz
}

function sw_gatk {
    local target="$1"

    # [ -d "$target" ] && rm -rf $target

    # download files
    wget --no-check-certificate https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip -O /tmp/gatk-4.2.6.1.zip

    # untar and prepare with build directory outside of source tree
    unzip /tmp/gatk-4.2.6.1.zip -d $target

    # remove spurious leftovers
    rm -rf /tmp/gatk-4.2.6.1.zip
}

function sw_sratoolkit {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    wget --no-check-certificate https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz -O /tmp/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz

    # untar and prepare with build directory outside of source tree
    mkdir -p $target
    tar -xf /tmp/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz -C $target
    # ${SOFTWARE}/sratoolkit/usr/local/ncbi/sra-tools/bin/vdb-config

    # remove spurious leftovers
    rm -rf /tmp/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz
}

function sw_igv {
    local target="$1"

    # [ -d "$target" ] && rm -rf $target

    local versionu="$(echo ${version} | awk -F. '{print $1"."$2}')"
    # download files
    wget --no-check-certificate https://data.broadinstitute.org/igv/projects/downloads/2.14/IGV_Linux_2.14.1_WithJava.zip -O /tmp/IGV_Linux_2.14.1_WithJava.zip

    # untar and prepare with build directory outside of source tree
    unzip /tmp/IGV_Linux_2.14.1_WithJava.zip -d $target

    # remove spurious leftovers
    rm -rf /tmp/IGV_Linux_2.14.1_WithJava.zip
}

function sw_trimmomatic {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    git clone https://github.com/usadellab/Trimmomatic.git $target

    # install
    cd $target
    ant
}

function sw_star {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    git clone https://github.com/alexdobin/STAR.git $target

    # install
    cd ${target}/source
    make -j4
}

function sw_picard {
    local target="$1"

    [ -d "$target" ] && rm -rf $target
    mkdir -p $target

    # download files
    # git clone https://github.com/broadinstitute/picard.git $target
    wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O ${target}/picard.jar

    # install
    # cd $target
    # ./gradlew shadowJar
}

function sw_htslib {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    git clone https://github.com/samtools/htslib.git $target

    # install
    cd $target
    git submodule update --init --recursive
    autoreconf -i
    ./configure
    make -j4
}

function sw_samtools {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    git clone https://github.com/samtools/samtools.git $target

    # install
    cd $target
    autoheader
    autoconf -Wno-syntax
    ./configure
    make -j4
}

function sw_bcftools {
    local target="$1"

    [ -d "$target" ] && rm -rf $target

    # download files
    git clone https://github.com/samtools/bcftools.git $target

    # install
    cd $target
    autoheader
    autoconf -Wno-syntax
    ./configure
    make -j4
}

################################################################################
################################################################################

# Make a couple of directories
mkdir -p ${SOFTWARE}

# Install each package
sw_openjdk      ${SOFTWARE}/java        # Install Java
sw_gatk         ${SOFTWARE}             # Install GATK
sw_igv          ${SOFTWARE}             # Install IGV
sw_sratoolkit   ${SOFTWARE}/sratoolkit  # Install SRA
sw_trimmomatic  ${SOFTWARE}/Trimmomatic # Install Trimmomatic
sw_star         ${SOFTWARE}/STAR        # Install STAR
sw_picard       ${SOFTWARE}/picard      # Install picard
sw_htslib       ${SOFTWARE}/htslib      # Install HTSlib
sw_samtools     ${SOFTWARE}/samtools    # Install samtools
sw_bcftools     ${SOFTWARE}/bcftools    # Install bcftools

# Make symbolic links
mkdir -p ${SOFTWARE}/bin
ln -s ${SOFTWARE}/STAR/source/STAR                                      ${SOFTWARE}/bin/STAR
ln -s ${SOFTWARE}/htslib/bgzip                                          ${SOFTWARE}/bin/bgzip
ln -s ${SOFTWARE}/samtools/samtools                                     ${SOFTWARE}/bin/samtools
ln -s ${SOFTWARE}/bcftools/bcftools                                     ${SOFTWARE}/bin/bcftools
ln -s ${SOFTWARE}/sratoolkit/usr/local/ncbi/sra-tools/bin/fasterq-dump  ${SOFTWARE}/bin/fasterq-dump
ln -s ${SOFTWARE}/htslib/tabix                                          ${SOFTWARE}/bin/tabix
