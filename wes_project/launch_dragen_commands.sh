#!/bin/bash

PATIENT_DIR=$1

if [[ -z ${PATIENT_DIR} || ! -d ${PATIENT_DIR} ]]; then
  echo "Patient directory doesn't exist: ${PATIENT_DIR}"
  exit 1
fi

echo "Launching commands"
for pDir in $(find ${PATIENT_DIR} -maxdepth 1 -mindepth 1 -type d); do
  dragen_dir=${pDir}/dragen
  for patient_dgn_dir in $(find ${dragen_dir} -mindepth 1 -maxdepth 1 -type d); do
    DGN_SCRIPT=$(find ${patient_dgn_dir} -type f -name "drgn_cmd*.sh")
    if [[  "$(echo $DGN_SCRIPT | tr ' ' '\n' | wc -l)" != "1" ]]; then
      echo "Found more than one dragen script: ${DGN_SCRIPT}"
      exit 1
    fi
    if [[ -z ${DGN_SCRIPT} || ! -f ${DGN_SCRIPT} ]]; then
      echo "Dragen script not found"
      exit 1
    fi
    dragen_size=$(du -sh ${patient_dgn_dir} | cut -f1)
    if [[ ! -z "$(echo ${dragen_size} | grep G)" ]]; then
      echo "Dragen Directory has size ${dragen_size}. May have been run through dragen already: ${patient_dgn_dir}"
      continue
    fi

    cd ${patient_dgn_dir}
    JOB_NAME=$(basename ${DGN_SCRIPT} | cut -d'.' -f1)
    CMD="bsub -J ${JOB_NAME} -o ${patient_dgn_dir}/${JOB_NAME}.out -n48 -q dragen ${DGN_SCRIPT}"
    echo $CMD
    eval $CMD
  done
done

echo "Done."
