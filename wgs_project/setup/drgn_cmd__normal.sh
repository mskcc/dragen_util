#!/bin/bash

NORMAL_SAMPLE=$1
FASTQ_LIST=$2
DGN_REF=$3

if [[ -z ${NORMAL_SAMPLE} || -z ${FASTQ_LIST} || -z ${DGN_REF} ]]; then
  echo "Please specify patient, normal sample, fastq_list.csv, & dragen reference"
  exit 1
fi

OUTPUT_PREFIX=${NORMAL_SAMPLE}_cnv

echo "#!/bin/bash"
echo ""
echo "/opt/edico/bin/dragen \\"
echo "  --ref-dir ${DGN_REF} \\"
echo "  --enable-duplicate-marking true \\"
echo "  --enable-map-align-output true \\"
echo "  --enable-variant-caller true \\"
echo "  --output-directory $(pwd) \\"
echo "  --output-file-prefix ${OUTPUT_PREFIX} \\"
echo "  --fastq-list-sample-id ${NORMAL_SAMPLE} \\"
echo "  --fastq-list ${FASTQ_LIST}"

