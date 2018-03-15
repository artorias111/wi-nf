# Use this script to BUILD the necessary SnpEff database!
# MUST USE A BUILD >= WS254
# Specify wormbase BUILD
set -e

BUILD="WS263"
GIT_PARENT_DIR=`git rev-parse --show-toplevel`

# Create directory
mkdir -p ${GIT_PARENT_DIR}/snpeff/${BUILD}

# Update config file
echo "${BUILD}.genome : C. elegans" >> ${GIT_PARENT_DIR}/snpeff/snpEff.config
echo "${BUILD}.MtDNA.codonTable : Invertebrate_Mitochondrial" >> ${GIT_PARENT_DIR}/snpeff/snpEff.config

# Download genome
wget -O ${GIT_PARENT_DIR}/snpeff/${BUILD}/sequences.fa.gz ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.WS253.genomic.fa.gz
# Extract sequence
zcat ${GIT_PARENT_DIR}/snpeff/${BUILD}/sequences.fa.gz > ${GIT_PARENT_DIR}/snpeff/${BUILD}/sequences.fa

# Download and extract protein fasta file
wget -O ${GIT_PARENT_DIR}/snpeff/${BUILD}/protein.fa.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${BUILD}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${BUILD}.protein.fa.gz 
zcat ${GIT_PARENT_DIR}/snpeff/${BUILD}/protein.fa.gz > ${GIT_PARENT_DIR}/snpeff/${BUILD}/protein.fa

# Download gtf
wget -O ${GIT_PARENT_DIR}/snpeff/${BUILD}/genes.gtf.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${BUILD}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${BUILD}.canonical_geneset.gtf.gz
# Build genome
snpEff BUILD -config ${GIT_PARENT_DIR}/snpeff/snpEff.config \
             -dataDir ${GIT_PARENT_DIR}/snpeff \
             -gtf22 -v ${BUILD}