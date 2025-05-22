#!/usr/bin/env Rscript

library(dplyr)

# auth <- sprintf("https://%s@github.com/cartographybio/binder_portal.git", token)
auth <- sprintf("https://%s@github.com/rscharf1/heparan-sulfates.git", token)

#List R Scripts
files <- list.files(pattern="\\.R$|\\.py$|\\.sh$|\\.tsv$|\\.csv$|\\.fa$|\\.txt$|\\.md$", recursive=TRUE)

size <- utils:::format.object_size(file.info(files)$size, units = "MB") %>% 
	{stringr::str_split(., pattern=" ", simplify=TRUE)[,1]} %>% 
		as.numeric

files <- files[size < 0.5]

#Add
for(i in seq_along(files)){
	system(sprintf("git add %s", files[i]))
}

#Pull
system(sprintf("git pull %s", auth))

#Time
time <- gsub(" |\\-|\\:", "_", Sys.time())

#Commit
system(sprintf("git commit -m '%s'", time))

#Push
system(sprintf("git push -u %s", auth))