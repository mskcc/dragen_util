#!/bin/bash

LOCATION=$(realpath $(dirname "$0"))

pairing_file=$1 	# TUMOR_ID	NORMAL_ID
fastq_list=$2		# RGID,RGSM,RGLB,Lane,Read1File,Read2File
dragen_reference=/staging/ref/GRCh37_cnv
fasta_reference=/home/igo/resources/bqsr/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta

if [[ ! -f ${pairing_file} || -z ${fastq_list} ]]; then
  echo "Please provide a pairing file & fastq_list.csv"
  printf "PAIRING FILE COLUMNS - \n\tTUMOR_ID	NORMAL_ID\n"
  printf "FASTQ_LIST - \n\tRGID,RGSM,RGLB,Lane,Read1File,Read2File\n"
  exit 1
fi

pairing_file=$(realpath ${pairing_file})
fastq_list=$(realpath ${fastq_list})

tail -n+2 ${pairing_file} | while read d; do
  tumor_sample=$(echo ${d} | cut -d' ' -f1)
  normal_sample=$(echo ${d} | cut -d' ' -f2)
  patient=$(echo ${tumor_sample} | cut -d'_' -f1,2,3)
  align_and_vcf_dir=${patient}/dragen/${tumor_sample}___${normal_sample}

  echo "Creating bam scripts for patient=${patient} normal_sample=${normal_sample} tumor_sample=${tumor_sample}"
  align_and_vcf_script="drgn_create_bams___${normal_sample}.sh"
  cd ${align_and_vcf_dir}
  ${LOCATION}/drgn_cmd__make_bams.sh ${normal_sample} ${tumor_sample} ${patient} ${fastq_list} ${dragen_reference} >> ${align_and_vcf_script}
  cd -

  echo "Creating bqsr script for patient=${patient} normal_sample=${normal_sample} tumor_sample=${tumor_sample}" 
  bqsr_dir=${patient}/bqsr
  mkdir -p ${bqsr_dir}
  
  normal_bqsr_script="drgn_create_bams___${normal_sample}.sh"

  # DRAGEN will create these files
  normal_bam=${align_and_vcf_dir}/${patient}.bam
  tumor_bam=${align_and_vcf_dir}/${patient}_tumor.bam
  cd ${bqsr_dir}
  ${LOCATION}/bqsr_cmd.sh ${normal_bam} ${fasta_reference} ${normal_sample} >> ${normal_bqsr_script}
  ${LOCATION}/bqsr_cmd.sh ${tumor_bam} ${fastq_reference} ${tumor_sample} >> ${tumor_bqsr_script}
  cd -
done
