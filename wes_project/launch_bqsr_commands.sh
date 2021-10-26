#!/bin/bash

PATIENT_DIR=$1

if [[ -z ${PATIENT_DIR} || ! -d ${PATIENT_DIR} ]]; then
  echo "Patient directory doesn't exist: ${PATIENT_DIR}"
  exit 1
fi

echo "Launching commands"
for pDir in $(find ${PATIENT_DIR} -maxdepth 1 -mindepth 1 -type d); do
  sample_dirs=$(find ${pDir} -mindepth 1 -maxdepth 1 -type d -name "s_C*")
  for sd in ${sample_dirs}; do
    bams=$(find ${sd} -type f -name "*.bam")
    SAMPLE=$(basename ${sd})
    if [[ -z ${bams} ]]; then
      echo "BAMs for sample not detected. Will run BQSR for ${SAMPLE}"
      bqsr_dir=${pDir}/bqsr
      sample_script=$(find ${bqsr_dir} -type f -name "*${SAMPLE}.sh")
      if [[ -z ${sample_script} ]]; then
        echo "Couldn't find bqsr script for ${SAMPLE}. Exiting..."
        exit 1
      else
        JOB_NAME=$(basename ${sample_script} | cut -d'.' -f1)
        OUT=${bqsr_dir}/${JOB_NAME}.out
        CMD="bsub -J ${JOB_NAME} -o ${OUT} -n 8 -M 6 $sample_script"
        echo ${CMD}
        eval ${CMD}
      fi
    else
      echo "Sample ${SAMPLE} has bams: ${bams}. Skipping..."
    fi
  done
done

echo "Done."
