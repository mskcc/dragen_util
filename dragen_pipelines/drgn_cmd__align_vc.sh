#!/bin/bash

REF=/staging/ref/GRCh37_dna
OUTPUT=$(pwd)
PFX=SAMPLE_NAME
N_RGSM=SAMPLE		# As listed in fastq_list.csv
T_RGSM=TUMOR_SAMPLE	# As listed in fastq_list.cs
FQ_LIST=fastq_list.csv
VC_TARGET_BED=target.bed

/opt/edico/bin/dragen \
  --ref-dir ${REF} \
  --enable-duplicate-marking true \
  --enable-map-align-output true \
  --enable-variant-caller true \
  --output-directory ${OUTPUT} \
  --output-file-prefix ${PFX} \
  --fastq-list-sample-id ${N_RGSM} \
  --fastq-list ${FQ_LIST} \
  --tumor-fastq-list-sample-id ${T_RGSM} \
  --tumor-fastq-list ${FQ_LIST} \
  --vc-target-bed ${VC_TARGET_BED} \
  --vc-target-bed-padding 5

