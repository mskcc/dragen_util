#!/bin/bash

B_ALLELE_VCF=sample.hard-filtered.vcf.gz 
PREFIX=sample
TUMOR_BAM=sample_tumor
SEX=Unknown # (FE)MALE

/opt/edico/bin/dragen -f \
  -r /staging/ref/GRCh37_cnv \
  --tumor-bam-input ${TUMOR_BAM} \
  --intermediate-results-dir /staging/intermediates \
  --output-directory $(pwd) \
  --output-file-prefix ${PREFIX} \
  --enable-map-align false \
  --enable-cnv true \
  --cnv-normal-b-allele-vcf ${B_ALLELE_VCF} \
  --sample-sex MALE
