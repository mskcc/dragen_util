# DRAGEN CNV/BQSR/WGS 
Creates scripts to run CNV/BQSR/WGS calling w/ DRAGEN

```
./setup_wgs_project.sh pairing_dragen.tsv fastq_list.csv
```
`pairing_dragen.tsv` has headers,
```
TUMOR_ID	NORMAL_ID	SEX
```

`fastq_list.csv` has headers,
```
RGID,RGSM,RGLB,Lane,Read1File,Read2File
```

This creates the following output directories & files,
```
patient_sample_id
├── bqsr
│   ├── drgn_create_bams___patient_sample_id_normal.sh		# 2) Creates BQSR BAMS
│   └── drgn_create_bams___patient_sample_id_tumor.sh
├── cnv
│   └── drgn_call_cnv___patient_sample_id.sh			# 3) Runs Somatic CNV calling on Tumor BQSR BAM
├── dragen
│   └── patient_sample_id_tumor___patient_sample_id_normal
│       └── drgn_create_bams___patient_sample_id.sh		# 1) Creates DRAGEN BAMS
├── patient_sample_id_normal
└── patient_sample_id_tumor

6 directories, 4 files
```
