#!/usr/bin/env bash

zcat ../gencode.vM39.annotation.gtf.gz | awk -F '\t' '
$3 == "transcript" &&
$9 ~ /tag "Ensembl_canonical"/ &&
$9 ~ /gene_type "protein_coding"/ &&
$9 ~ /transcript_type "protein_coding"/
' > genes.gtf



awk -F '\t' '
{
    match($9, /gene_id "[^"]+"/)
    print substr($9, RSTART + 9, RLENGTH - 10)
}' genes.gtf |
sort -u |
wc -l