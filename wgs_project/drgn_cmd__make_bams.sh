#!/bin/bash

NORMAL_SAMPLE=$1
TUMOR_SAMPLE=$2
PATIENT=$3
FASTQ_LIST=$4
DGN_REF=$5

if [[ -z ${NORMAL_SAMPLE} || -z ${TUMOR_SAMPLE} || -z ${PATIENT} || -z ${FASTQ_LIST} || -z ${DGN_REF} ]]; then
  echo "Please specify patient, normal sample, fastq_list.csv, & dragen reference"
  exit 1
fi

echo "#!/bin/bash"
echo ""
echo "/opt/edico/bin/dragen \\"
echo "  --ref-dir ${DGN_REF} \\"
echo "  --enable-duplicate-marking true \\"
echo "  --enable-map-align-output true \\"
echo "  --enable-variant-caller true \\"
echo "  --output-directory $(pwd) \\"
echo "  --output-file-prefix ${PATIENT} \\"
echo "  --fastq-list-sample-id ${NORMAL_SAMPLE} \\"
echo "  --fastq-list ${FASTQ_LIST} \\"
echo "  --tumor-fastq-list-sample-id ${TUMOR_SAMPLE} \\"
echo "  --tumor-fastq-list ${FASTQ_LIST}"

