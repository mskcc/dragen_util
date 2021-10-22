#!/bin/bash

REF=...
LIFTOVER=...
OUTPUT=...

/opt/edico/bin/dragen --build-hash-table true \
  --enable-cnv true \
  --ht-reference ${REF}
  --output-directory  ${OUTPUT}
