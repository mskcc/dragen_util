# Liftover Sam 

## Description
Create a liftover SAM file of an input reference file. A liftover SAM file is an alignment of a reference's unlocalized scaffolds to their corresponding primary assembly chromosomes.
* This file can be used for alt-aware alignment in DRAGEN

## Setup
Install dependencies: `faSplit`
```
conda create --name liftover_sam --file requirements.txt
conda activate liftover_sam
```

## Run
```
$ ./find_scaffolds_to_align.sh -r GRCh37.fasta
```
