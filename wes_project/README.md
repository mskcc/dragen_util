# Samples for DRAGEN

## Dependencies
* LSF Cluster
* Populated `pairs_with_igoIds.txt` `sample_fastq.txt` `fastq_list.csv` 

## Steps
1) Create Patient Setup

Will output a `patients` directory
```
./create_patient_setup.sh pairs_with_igoIds.txt sample_fastq.txt fastq_list.csv $(pwd)
```

2) Launch dragen commands
```
PATIENTS_DIR=./patients
./launch_dragen_commands.sh ${PATIENTS_DIR} 
```

3) Launch BQSR commands (reads from dragen output)
```
PATIENTS_DIR=./patients
./launch_bqsr_commands.sh ${PATIENTS_DIR}
```

4) Extract time metrics
```
PATIENTS_DIR=./patients
./create_time_summary.sh ${PATIENTS_DIR}
```
