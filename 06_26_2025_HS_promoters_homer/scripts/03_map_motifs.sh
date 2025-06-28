#!/bin/bash

annotatePeaks.pl outputs/HS_promoters_mouse.bed mm10 \
  -m outputs/homer/knownResults/known*.motif \
  > outputs/motif_hits.txt