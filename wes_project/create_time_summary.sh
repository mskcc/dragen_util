#!/bin/bash

COL_REF_LOAD="Time loading reference"
COL_ALIGN="Time aligning reads"
COL_DUP_MARK="Time duplicate marking"
COL_SORT_MD="Time sorting and marking duplicates"
COL_DRAGStr="Time DRAGStr calibration"
COL_MAP_ALIGN="Time saving map/align output"
COL_RECONFIG="Time partial reconfiguration"
COL_VAR_CALL="Time variant calling"
COL_PARTITION="Time partitioning"
COL_TOTAL="Total runtime"

if [[ -z ${PATIENT_DIR} || ! -d ${PATIENT_DIR} ]]; then
  echo "Patient directory doesn't exist: ${PATIENT_DIR}"
  exit 1
fi

echo "tumor___normal pair,${COL_REF_LOAD},${COL_ALIGN},${COL_DUP_MARK},${COL_SORT_MD},${COL_DRAGStr},${COL_MAP_ALIGN},${COL_RECONFIG},${COL_VAR_CALL},${COL_PARTITION},${COL_TOTAL},tumor_bqsr_time,normal_bqsr_time"
TIME_FILES=$(find ${PATIENT_DIR} -type f -name "*.time_metrics.csv")
for f in ${TIME_FILES}; do
  patient_dir=$(echo ${f} | grep -oP ".*patients/s_C_[0-9A-Z]{6}")
  patient=$(basename $(dirname ${f}))
  ref_load=$(grep "${COL_REF_LOAD}" ${f} | cut -d',' -f5)
  align=$(grep "${COL_ALIGN}" ${f} | cut -d',' -f5)
  dup_mark=$(grep "${COL_DUP_MARK}" ${f} | cut -d',' -f5)
  sort_md=$(grep "${COL_SORT_MD}" ${f} | cut -d',' -f5)
  dragstr=$(grep "${COL_DRAGStr}" ${f} | cut -d',' -f5)
  map_align=$(grep "${COL_MAP_ALIGN}" ${f} | cut -d',' -f5)
  reconfig=$(grep "${COL_RECONFIG}" ${f} | cut -d',' -f5)
  var_call=$(grep "${COL_VAR_CALL}" ${f} | cut -d',' -f5)
  partition=$(grep "${COL_PARTITION}" ${f} | cut -d',' -f5)
  total=$(grep "${COL_TOTAL}" ${f} | cut -d',' -f5)

  sample=$(echo ${patient} | sed 's/___.*//g')
  normal_sample=$(echo ${patient} | sed 's/.*___//g')
  tumor_bqsr=$(find ${patient_dir}/bqsr/ -type f -name "bqsr___*.out" | grep ${sample})
  normal_bqsr=$(find ${patient_dir}/bqsr/ -type f -name "bqsr___*.out" | grep ${normal_sample})
  bqsr_entries=$(cat ${tumor_bqsr} | grep "Run time" | tail -1 | rev | cut -d' ' -f2 | rev)
  bqsr_entries+=",$(cat ${normal_bqsr} | grep "Run time" | tail -1 | rev | cut -d' ' -f2 | rev)"
  echo "${patient},${ref_load},${align},${dup_mark},${sort_md},${dragstr},${map_align},${reconfig},${var_call},${partition},${total},${bqsr_entries}"
done
