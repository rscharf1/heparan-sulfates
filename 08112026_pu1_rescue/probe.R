library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

dt <- fread("expr_table_clean.txt")

map <- fread("probe_map_clean.txt")

# First join the probe annotation
d <- merge(
  dt,
  map[, .(ID, Symbol)],
  by.x = "ID_REF",
  by.y = "ID"
)

# Mean expression within each condition
d[, PuER := rowMeans(.SD),
  .SDcols = c("GSM1563673", "GSM1563674", "GSM1563675")]

d[, EV := rowMeans(.SD),
  .SDcols = c("GSM1563676", "GSM1563677", "GSM1563678")]

d[, .(ID_REF, Symbol, PuER, EV)] |> head()

pdf("out.pdf")

ggplot(d, aes(x = PuER, y = EV)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "PuER expression",
    y = "Empty vector expression"
  )

dev.off()

to_highlight <- c(
  "PU1"  = "Sfpi1",
  "EXT1" = "Ext1",
  "EXT2" = "Ext2",
  "EXTL2" = "Extl2",
  "EXTL3" = "Extl3",
  "B3GAT2" = "B3gat2",
  "B3GAT3" = "B3gat3",
  "B4GALT1" = "B4galt1",
  "B4GALT3" = "B4galt3",
  "B4GALT7" = "B4galt7",
  "FAM20B" = "Fam20b",
  "GLCE" = "Glce",
  "XYLT2" = "Xylt2"
)

# Make a small table for highlighted genes
hl <- data.table(
  label  = names(to_highlight),
  Symbol = unname(to_highlight)
)

hl <- merge(hl, d, by = "Symbol")

pdf("out2.pdf")

ggplot(d, aes(x = PuER, y = EV)) +
  
  geom_point(
    alpha = 0.25,
    size = 0.8
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  
  # Highlighted genes
  geom_point(
    data = hl,
    color = "red",
    size = 3
  ) +
  
  # Labels
  geom_text_repel(
    data = hl,
    aes(label = label),
    color = "red",
    fontface = "bold"
  ) +
  
  coord_equal() +
  theme_classic() +
  labs(
    x = "PuER expression",
    y = "EV expression"
  )

dev.off()

# Track down codes 
d[grepl("XYLT", Symbol, ignore.case=TRUE)]

map[grepl("pu", Symbol, ignore.case=TRUE)]


map[grepl("Spi1|PU\\.1|Sfpi1", 
          paste(ID, Transcript, Symbol),
          ignore.case = TRUE)]

##########
d[, diff := PuER - EV]

# Sort highest -> lowest
setorder(d, -diff)

# Rank for plotting
d[, rank := .I]

hl <- d[Symbol %in% to_highlight]

# Add the display labels
hl[, label := names(to_highlight)[
  match(Symbol, unname(to_highlight))
]]

pdf("out3.pdf", height = 5)

	ggplot(d, aes(x = rank, y = diff)) +
	  
	  # All probes
	  geom_point(
	    alpha = 0.25,
	    size = 0.8
	  ) +
	  
	  # No-change reference
	  geom_hline(
	    yintercept = 0,
	    linetype = 2
	  ) +
	  
	  # Genes of interest on top
	  geom_point(
	    data = hl,
	    color = "red",
	    size = 3
	  ) +
	  
	  geom_text_repel(
	    data = hl,
	    aes(label = label),
	    color = "red",
	    fontface = "bold"
	  ) +
	  
	  theme_classic() +
	  labs(
	    x = "Probe rank",
	    y = "PuER - EV expression"
	  )

dev.off()

N <- nrow(d)

hl[, percentile := 100 * (N - rank + 1) / N]

hl[1:5, .(
  label,
  Symbol,
  diff = round(diff, 3),
  rank,
  percentile = round(percentile, 2)
)]











