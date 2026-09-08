library(data.table)
library(dplyr)
library(ggplot2)

# promoter capture hi-c to identify enhancers interacting with promoters
	# bait is RNA probes that hybridize to promoter regions
	# other... is the region that was bound to the promoter
hic <- fread("inputs/HPC7_Promoter_Capture_Interactions.ibed")

hs_genes <- fread("inputs/HS_relevant_genes.csv")

hic[grep("Ndst", bait_name)]

hs_hic <- hic[grepl(paste(hs_genes$hgnc_symbol, collapse="|"), bait_name, ignore.case = TRUE)]

hs_hic[order(-N_reads)]
hs_hic[order(-score)]

# chip to look at active regions of the genome
chip <- fread("inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph")

chip[V1 == "chr1" & V2 >= 90172816 & V3 <= 90175152][V4 >= 100]

90172816 - 90175152

pdf("outputs/hist.pdf")
hist(chip$V4, breaks = 100, main = "H3K27ac Signal Distribution", xlab = "Signal")
dev.off()