#!/bin/bash

B_ALLELE_VCF=$1
TUMOR_BAM=$2
PATIENT=$3
SEX=$4
DGN_REF=$5

if [[ ! -f ${B_ALLELE_VCF} || ! -f ${TUMOR_BAM} || -z ${PATIENT} || -z ${DGN_REF} || -z ${SEX} ]]; then
  echo "Please provide a *hard-filtered.vcf, tumor bam, patient, sex, & dragen reference"
  exit 1
fi

PREFIX="${PATIENT}_cnv"

echo "/opt/edico/bin/dragen -f \\"
echo "  -r ${DGN_REF} \\"
echo "  --tumor-bam-input ${TUMOR_BAM} \\"
echo "  --intermediate-results-dir /staging/intermediates \\"
echo "  --output-directory $(pwd) \\"
echo "  --output-file-prefix ${PREFIX} \\"
echo "  --enable-map-align false \\"
echo "  --enable-cnv true \\"
echo "  --cnv-normal-b-allele-vcf ${B_ALLELE_VCF} \\"
echo "  --sample-sex ${SEX}"
