#!/bin/bash

md5sum /igo/work/nabors/genomes/GRCh37/GRCh37.fasta # 12a0bed94078e2d9e8c00da793bbc84e  

/opt/edico/bin/dragen --build-hash-table true \
  --ht-reference /igo/work/nabors/genomes/GRCh37/GRCh37.fasta \
  --ht-alt-liftover /staging/ref/GRCh37_alt_aware/hs37d5_liftover.sam \
  --output-directory  /staging/ref/GRCh37_alt_aware >> log___create_ref___hs37d5.out
