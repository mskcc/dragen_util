#!/bin/bash

/opt/edico/bin/dragen --build-hash-table true \
  --enable-cnv true \
  --ht-reference /igo/work/nabors/genomes/GRCh37/GRCh37.fasta \
  --output-directory  /staging/ref/GRCh37_cnv >> log___create_ref___GRCh37_cnv.out
