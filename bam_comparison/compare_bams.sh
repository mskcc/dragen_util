ORIGINAL_INPUT=$1
TEST_INPUT=$2
BAM_DIFF_OUT=bam_diff.out

if [[ ! -f ${ORIGINAL_INPUT} || ! -f ${TEST_INPUT} ]]; then
  echo "INVALID INPUTS"
  exit 1
fi

ORIGINAL_INPUT=$(realpath ${ORIGINAL_INPUT})
TEST_INPUT=$(realpath ${TEST_INPUT})

ORIGINAL_NAME="$(basename ${ORIGINAL_INPUT} | cut -d'.' -f1)_original"
TEST_NAME="$(basename ${TEST_INPUT} | cut -d'.' -f1)_test"

ORIGINAL_SORTED=${ORIGINAL_NAME}_sorted.bam
TEST_SORTED=${TEST_NAME}_sorted.bam

echo "Creating ORIGINAL=${ORIGINAL_SORTED} and TEST=${TEST_SORTED}"
SORT1="samtools sort ${ORIGINAL_INPUT} -o ${ORIGINAL_SORTED}"
echo ${SORT1}
eval ${SORT1}

SORT2="samtools sort ${TEST_INPUT} -o ${TEST_SORTED}"
echo ${SORT2}
eval ${SORT2}

bam diff --in1 ${ORIGINAL_SORTED} --in2 ${TEST_SORTED} --mapQual >> ${BAM_DIFF_OUT}

python /igo/home/streidd/working/bioinformatics/bam_comparison/bam_util_to_csv.py ${BAM_DIFF_OUT}
