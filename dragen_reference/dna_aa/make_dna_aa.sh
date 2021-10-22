#!/bin/bash

REF=...
LIFTOVER=...
OUTPUT=...

/opt/edico/bin/dragen --build-hash-table true \
  --ht-reference ${REF} \
  --ht-alt-liftover ${LIFTOVER} \
  --output-directory  ${OUTPUT}
