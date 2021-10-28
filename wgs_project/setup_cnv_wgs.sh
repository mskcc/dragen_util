#!/bin/bash

LOCATION=$(realpath $(dirname "$0"))

pairing_file=$1 	# TUMOR_ID	NORMAL_ID
fastq_list=$2		# RGID,RGSM,RGLB,Lane,Read1File,Read2File
dragen_reference=/staging/ref/GRCh37_cnv

if [[ ! -f ${pairing_file} || -z ${fastq_list} ]]; then
  echo "Please provide a pairing file & fastq_list.csv"
  exit 1
fi

pairing_file=$(realpath ${pairing_file})
fastq_list=$(realpath ${fastq_list})

tail -n+2 ${pairing_file} | while read d; do
  tumor_sample=$(echo ${d} | cut -d' ' -f1)
  normal_sample=$(echo ${d} | cut -d' ' -f2)
  patient=$(echo ${tumor_sample} | cut -d'_' -f1,2,3)

  cnv_wgs_dir=${patient}/dragen/cnv___${tumor_sample}___${normal_sample}
  mkdir -p ${cnv_wgs_dir}
  
  echo "Creating cnv/wgs for patient=${patient} tumor_sample=${tumor_sample} normal_sample=${normal_sample}"
  normal_vcf_script="drgn_create_cnv_wgs___${normal_sample}.sh"
  
  ${LOCATION}/drgn_cmd__normal.sh ${normal_sample} ${fastq_list} ${dragen_reference} >> ${normal_align_vcf_dir}/${normal_vcf_script}
done
