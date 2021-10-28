#!/bin/bash

INPUT_BAM=$1
REFERENCE=$2
SAMPLE_NM=$3
RECAL_TBL=$(echo ${INPUT_BAM} | cut -d'.' -f1).recal.table

BQSR_SCRIPT=bqsr___${SAMPLE_NM}.sh

echo "/igo/home/igo/resources/gatk-4.1.9.0/gatk BaseRecalibrator \\" >> ${BQSR_SCRIPT}
echo "  --java-options '-Xmx28g' \\" >> ${BQSR_SCRIPT}
echo "  --reference ${REFERENCE} \\" >> ${BQSR_SCRIPT}
echo "  --known-sites /home/igo/resources/bqsr/dbsnp_138_b37/dbsnp_138.b37.vcf \\" >> ${BQSR_SCRIPT}
echo "  --known-sites /home/igo/resources/bqsr/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf \\" >> ${BQSR_SCRIPT}
echo "  --known-sites /home/igo/resources/bqsr/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf \\" >> ${BQSR_SCRIPT}
echo "  --verbosity INFO \\" >> ${BQSR_SCRIPT}
echo "  --input ${INPUT_BAM} \\" >> ${BQSR_SCRIPT}
echo "  --output ${RECAL_TBL}" >> ${BQSR_SCRIPT}

BQSR_BAM_PATH=$(realpath .)/../${SAMPLE_NM}
mkdir -p ${BQSR_BAM_PATH}
BQSR_BAM_FILE_PATH=${BQSR_BAM_PATH}/${SAMPLE_NM}.bam

echo "/igo/home/igo/resources/gatk-4.1.9.0/gatk ApplyBQSR \\" >> ${BQSR_SCRIPT}
echo "  --java-options '-Xmx28g' \\" >> ${BQSR_SCRIPT}
echo "  --reference ${REFERENCE} \\" >> ${BQSR_SCRIPT}
echo "  --create-output-bam-index true \\" >> ${BQSR_SCRIPT}
echo "  --bqsr-recal-file ${RECAL_TBL} \\" >> ${BQSR_SCRIPT}
echo "  --input ${INPUT_BAM} \\" >> ${BQSR_SCRIPT}
echo "  --output ${BQSR_BAM_FILE_PATH}" >> ${BQSR_SCRIPT}

