#!/bin/bash

INPUT_BAM=$1
REFERENCE=$2
SAMPLE_NM=$3
OUTPT_BAM=$4
RECAL_TBL=$(echo ${INPUT_BAM} | cut -d'.' -f1).recal.table

BQSR_SCRIPT=bqsr___${SAMPLE_NM}.sh

echo "/igo/home/igo/resources/gatk-4.1.9.0/gatk BaseRecalibrator \\"
echo "  --java-options '-Xmx28g' \\"
echo "  --reference ${REFERENCE} \\"
echo "  --known-sites /home/igo/resources/bqsr/dbsnp_138_b37/dbsnp_138.b37.vcf \\"
echo "  --known-sites /home/igo/resources/bqsr/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf \\"
echo "  --known-sites /home/igo/resources/bqsr/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf \\"
echo "  --verbosity INFO \\"
echo "  --input ${INPUT_BAM} \\"
echo "  --output ${RECAL_TBL}"

echo "/igo/home/igo/resources/gatk-4.1.9.0/gatk ApplyBQSR \\"
echo "  --java-options '-Xmx28g' \\"
echo "  --reference ${REFERENCE} \\"
echo "  --create-output-bam-index true \\"
echo "  --bqsr-recal-file ${RECAL_TBL} \\"
echo "  --input ${INPUT_BAM} \\"
echo "  --output ${OUTPT_BAM}"

