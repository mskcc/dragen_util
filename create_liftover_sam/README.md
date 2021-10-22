# Liftover Sam 

## Description
Create a liftover SAM file of an input reference file. A liftover SAM file is an alignment of a reference's unlocalized scaffolds to their corresponding primary assembly chromosomes.
* This file can be used for alt-aware alignment in DRAGEN

### Outline 
1. Download hs37d5.fa.gz reference (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/)
2. `faSplit` reference into its scaffolds
3. `bwa idx/mem` unlocalized scaffold to primary assembly (e.g. `/bwa mem -M -t 40 fasta/1.fa fasta/GL000191.1.fa`)
4. concatenate output SAM entries (excluding headers) to `hs37d5_liftover.sam`

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
