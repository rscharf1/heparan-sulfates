gunzip -c ../GRCm39.primary_assembly.genome.fa.gz \
  > ../GRCm39.primary_assembly.genome.fa

samtools faidx ../GRCm39.primary_assembly.genome.fa

cut -f1,2 ../GRCm39.primary_assembly.genome.fa.fai \
  > ../GRCm39.genome

head GRCm39.genome

bedtools slop \
  -i genes_tss.bed \
  -g ../GRCm39.genome \
  -s \
  -l 1000 \
  -r 200 \
  > genes_promoters.bed

wc -l genes_promoters.bed

awk '{print $3 - $2}' genes_promoters.bed |
sort -n |
uniq -c

bedtools getfasta \
  -fi ../GRCm39.primary_assembly.genome.fa \
  -bed genes_promoters.bed \
  -s \
  -nameOnly \
  > mouse_protein_coding_promoters.fa