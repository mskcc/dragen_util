#!/bin/bash

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" ]]; then
  echo "Please provide path to pairs_with_igoIds.txt, sample_fastq.txt, fastq_list.csv files, & output directory"
  exit 1
fi

IDT_BED=
AGILENT_BED=
REFERENCE=
KNOWN_SITES_1=
KNOWN_SITED_2=
KNOWN_SITES_3=

SAMPLE_PAIRS_FILE="$1" 
DELIVERED_FASTQS="$2"
FASTQ_LIST="$3"
OUTPUT_DIR="$4"

echo "SAMPLE_PAIRS_FILE=${SAMPLE_PAIRS_FILE} DELIVERED_FASTQS=${DELIVERED_FASTQS} FASTQ_LIST=${FASTQ_LIST} OUTPUT_DIR=${OUTPUT_DIR}"

function check_no_split_fastqs {
  SAMPLE_ID=$1

  NUM_FASTQS=$(cat ${FASTQ_LIST} | grep ${SAMPLE_ID} | wc -l)
  EXPECTED=$(cat ${DELIVERED_FASTQS} | grep ${SAMPLE_ID} | wc -l)

  # echo "NUM_FASTQS=${NUM_FASTQS} EXPECTED=${EXPECTED}"
  if [[ "${NUM_FASTQS}" -ne "${EXPECTED}" ]]; then
    echo "${SAMPLE_ID} doesn't have the expected amount"
    exit 1
  fi
}

#######################################
# Determines what the patient ID is from the Tumor & Normal IDs
# Globals:
#   None
# Arguments:
#   TUMOR_ID
#   NORMAL_ID
#######################################
function get_patient_id {
  TUMOR_ID=$1
  NORMAL_ID=$2

  PATIENT_REGEX="s_C_[0-9,A-Z]+"
  
  TUMOR_PATIENT=$(echo $TUMOR_ID | grep -oP ${PATIENT_REGEX})
  NORMAL_PATIENT=$(echo $NORMAL_ID | grep -oP ${PATIENT_REGEX})

  if [[ (${TUMOR_PATIENT} == ${NORMAL_ID}) && (! -z ${TUMOR_PATIENT}) ]]; then
    echo "${TUMOR_PATIENT} != ${NORMAL_ID}"
    exit 1
  fi

  echo $TUMOR_PATIENT
}

function get_bait_file {
  SAMPLE_ID=$1

  target=$(cat ${DELIVERED_FASTQS} | grep ${SAMPLE_ID} | cut -f2 | sort | uniq)
  if [[ "${target}" == "idt" ]]; then
    echo "${IDT_BED}"
  elif [[ "${target}" == "agilent" ]]; then
    echo "${AGILENT_BED}"
  else
    exit 1
  fi
}

function write_bqsr_cmd {
  SAMPLE_ID=$1
  INPUT_BAM=$2
  PARENT_DIR=$3

  BQSR_OUTPUT_DIR=${PARENT_DIR}/bqsr
  mkdir -p ${BQSR_OUTPUT_DIR}
  BQSR_CMD=${BQSR_OUTPUT_DIR}/bqsr___${SAMPLE_ID}.sh
  RECAL_TABLE=${BQSR_OUTPUT_DIR}/${SAMPLE_ID}.recal.table

  SAMPLE_DIR=${PARENT_DIR}/${SAMPLE_ID}
  mkdir -p ${SAMPLE_DIR}
  RECAL_BAM=${SAMPLE_DIR}/${SAMPLE_ID}.bam

  echo "#!/bin/bash" > ${BQSR_CMD}
  echo "" >> ${BQSR_CMD}
  echo "gatk BaseRecalibrator \\" >> ${BQSR_CMD}
  echo "  --java-options '-Xmx28g' \\" >> ${BQSR_CMD}
  echo "  --reference ${REFERENCE} \\" >> ${BQSR_CMD}
  echo "  --known-sites ${KNOWN_SITES_1} \\" >> ${BQSR_CMD}
  echo "  --known-sites ${KNOWN_SITES_2} \\" >> ${BQSR_CMD}
  echo "  --known-sites ${KNOWN_SITES_3} \\" >> ${BQSR_CMD}
  echo "  --verbosity INFO \\" >> ${BQSR_CMD}
  echo "  --input ${INPUT_BAM} \\" >> ${BQSR_CMD}
  echo "  --output ${RECAL_TABLE} && \\" >> ${BQSR_CMD}
  echo "gatk ApplyBQSR \\" >> ${BQSR_CMD}
  echo "  --java-options '-Xmx28g' \\" >> ${BQSR_CMD}
  echo "  --reference ${REFERENCE} \\" >> ${BQSR_CMD}
  echo "  --create-output-bam-index true \\" >> ${BQSR_CMD}
  echo "  --bqsr-recal-file ${RECAL_TABLE} \\" >> ${BQSR_CMD}
  echo "  --input ${INPUT_BAM} \\" >> ${BQSR_CMD}
  echo "  --output ${RECAL_BAM}" >> ${BQSR_CMD}

  chmod 755 ${BQSR_CMD}
}

function write_wrapper {
  PARENT_OUTPUT_DIR=$1
  TUMOR_ID=$2
  NORMAL_ID=$3
  PATIENT_ID=$4
  PARENT_FASTQ_LIST=$5

  # Setup DRAGEN command & directory
  PATIENT_OUTPUT=${PARENT_OUTPUT_DIR}/${PATIENT_ID}
  DRAGEN_OUTPUT_DIR=${PATIENT_OUTPUT}/dragen/${TUMOR_ID}___${NORMAL_ID}
  mkdir -p ${DRAGEN_OUTPUT_DIR}

  PATIENT_FASTQ_LIST=${DRAGEN_OUTPUT_DIR}/fastq_list.csv
  DRAGEN_CMD=${DRAGEN_OUTPUT_DIR}/drgn_cmd__${PATIENT_ID}.sh

  echo "RGPL,RGID,RGSM,RGLB,Lane,Read1File,Read2File" > ${PATIENT_FASTQ_LIST}
  cat ${PARENT_FASTQ_LIST} | grep ",${TUMOR_ID}," >> ${PATIENT_FASTQ_LIST}
  cat ${PARENT_FASTQ_LIST} | grep ",${NORMAL_ID}," >> ${PATIENT_FASTQ_LIST}

  TUMOR_SAMPLE_BAIT=$(get_bait_file ${TUMOR_ID})
  NORMAL_SAMPLE_BAIT=$(get_bait_file ${NORMAL_ID})
  if [ "${TUMOR_SAMPLE_BAIT}" != "${NORMAL_SAMPLE_BAIT}" ]; then
    echo "TUMOR_SAMPLE_BAIT != NORMAL_SAMPLE_BAIT (${TUMOR_SAMPLE_BAIT} != ${NORMAL_SAMPLE_BAIT})"
    exit 1
  fi
  PATIENT_BAIT=${TUMOR_SAMPLE_BAIT}

  if [[ ! -z $(echo $PATIENT_BAIT | grep IDT) ]]; then
    echo "IDT: ${PATIENT_BAIT}. Skipping..."
    return 0
  else
    echo "Agilent: ${PATIENT_BAIT}. Will continue"
  fi

  echo "#!/bin/bash" > ${DRAGEN_CMD}
  echo "" >> ${DRAGEN_CMD}
  echo "/opt/edico/bin/dragen \\" >> ${DRAGEN_CMD}
  echo "  --ref-dir /staging/ref/GRCh37_dna \\" >> ${DRAGEN_CMD}
  echo "  --enable-duplicate-marking true \\" >> ${DRAGEN_CMD}
  echo "  --enable-map-align-output true \\" >> ${DRAGEN_CMD}
  echo "  --enable-variant-caller true \\" >> ${DRAGEN_CMD}
  echo "  --output-directory $DRAGEN_OUTPUT_DIR \\" >> ${DRAGEN_CMD}
  echo "  --output-file-prefix ${PATIENT_ID} \\" >> ${DRAGEN_CMD}
  echo "  --fastq-list-sample-id ${NORMAL_ID} \\" >> ${DRAGEN_CMD}
  echo "  --fastq-list ${PATIENT_FASTQ_LIST} \\" >> ${DRAGEN_CMD}
  echo "  --tumor-fastq-list-sample-id ${TUMOR_ID} \\" >> ${DRAGEN_CMD}
  echo "  --tumor-fastq-list ${PATIENT_FASTQ_LIST} \\" >> ${DRAGEN_CMD}
  echo "  --vc-target-bed ${PATIENT_BAIT} \\" >> ${DRAGEN_CMD}
  echo "  --vc-target-bed-padding 5" >> ${DRAGEN_CMD}
  echo "" >> ${DRAGEN_CMD}

  chmod 755 ${DRAGEN_CMD}

  DGN_NORMAL_BAM=${DRAGEN_OUTPUT_DIR}/${PATIENT_ID}.bam
  DGN_TUMOR_BAM=${DRAGEN_OUTPUT_DIR}/${PATIENT_ID}_tumor.bam

  write_bqsr_cmd ${NORMAL_ID} ${DGN_NORMAL_BAM} ${PATIENT_OUTPUT}
  write_bqsr_cmd ${TUMOR_ID} ${DGN_TUMOR_BAM} ${PATIENT_OUTPUT}
}

while read row; do 
  TUMOR_ID=$(echo $row | cut -d' ' -f1)
  NORMAL_ID=$(echo $row | cut -d' ' -f2)
  TUMOR_IGO_ID=$(echo $row | cut -d' ' -f3)
  NORMAL_IGO_ID=$(echo $row | cut -d' ' -f4)

  PATIENT_ID=$(get_patient_id ${TUMOR_ID} ${NORMAL_ID})

  echo "PATIENT_ID=${PATIENT_ID} TUMOR_ID=${TUMOR_ID} NORMAL_ID=${NORMAL_ID}"

  TUMOR_SAMPLE_BAIT=$(get_bait_file ${TUMOR_ID})
  NORMAL_SAMPLE_BAIT=$(get_bait_file ${NORMAL_ID})
  if [ "${TUMOR_SAMPLE_BAIT}" != "${NORMAL_SAMPLE_BAIT}" ]; then
    echo "TUMOR_SAMPLE_BAIT != NORMAL_SAMPLE_BAIT (${TUMOR_SAMPLE_BAIT} != ${NORMAL_SAMPLE_BAIT})"
    exit 1
  fi
  PATIENT_BAIT=${TUMOR_SAMPLE_BAIT}

  agilent_file=$(basename ${AGILENT_BED})
  is_agilent=$(echo ${PATIENT_BAIT} | grep "${agilent_file}")
  if [[ ! -z ${is_agilent} ]]; then
    echo "Recreating Agilent (${PATIENT_ID}): ${PATIENT_BAIT}"
  else
    echo "[SKIPPING] Not Agilent (${PATIENT_ID}): ${PATIENT_BAIT}"
    continue
  fi
  echo "Running on ${PATIENT_ID}"
  write_wrapper ${OUTPUT_DIR} ${TUMOR_ID} ${NORMAL_ID} ${PATIENT_ID} ${FASTQ_LIST}
done < <(tail -n +2 $SAMPLE_PAIRS_FILE)

