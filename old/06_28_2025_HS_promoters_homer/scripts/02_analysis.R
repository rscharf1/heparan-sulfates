library(data.table)
library(dplyr)

# Look in knownResults.txt, try to pull out these positive controls 

dt <- fread("outputs/chol_genes/knownResults.txt")

dt[dt[[6]] > 0]

dt[[6]] %>% head

dt[order(-dt[[6]])]

dt[order(dt[[4]])] %>% head