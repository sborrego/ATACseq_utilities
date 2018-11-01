# UNIX
# cat *.txt | sed 's/\://g' > combo_stats.txt

library(dplyr)

wd <- "/Users/stacey/Data/ATACseq/atac_data"

setwd(wd)

combo_stats <- read.delim("combo_stats2.txt", header = FALSE)

combo_stats$cat <- sub("(\\w+).*", "\\1", combo_stats$V1)
combo_stats$data <- sub("\\w+", "\\1", combo_stats$V1)

files <- filter(combo_stats, cat=="Files") %>% select(data)
barcode <- filter(combo_stats, cat=="Barcode") %>% select(data)
reads <- filter(combo_stats, cat=="Reads") %>% select(data)
lib <- filter(combo_stats, cat=="Library") %>% select(data)
cycles1 <- filter(combo_stats, cat=="Cycles1") %>% select(data)
cycles2 <- filter(combo_stats, cat=="Cycles2") %>% select(data)

sample_list <- data.frame(files,
                          barcode,
                          reads,
                          lib,
                          cycles1,
                          cycles2)
colnames(sample_list) <- c("File",
                           "Barcode",
                           "Reads",
                           "Library",
                           "Cycle1",
                           "Cycle2")

