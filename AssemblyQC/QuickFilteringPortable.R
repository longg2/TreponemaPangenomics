#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)

if(length(args) < 4){
	cat("Usage: ./QuickFilteringPortable.R <Quast Transposed TSV Report> <CheckM Lineage Output> <CheckM Marker Lineage> <OutputFile>\n")
	stop("Please provide the required options",call. = F)
}else{
	quastTransposed <- args[1]
	checkMdf <- args[2]
	lineageSearch <- args[3]
	output <- args[4]
}

### Making sure that the packages are installed
packages <- c("dplyr", "tidyr", "ggplot2", "ggpubr")
install.packages(setdiff(c("dplyr", "tidyr", "ggplot2", "ggpubr"),rownames(installed.packages())))

### Loading the packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

### Loading in the data
cat("Loading the data\n")
quastdf <- read.delim(quastTransposed) |> as_tibble()
colnames(quastdf) <- gsub("X..","",colnames(quastdf))

quastdf <- quastdf |>
  rename(Genome.Fraction = Genome.fraction....) |>
  mutate(Assembly = gsub("_contigs", "", Assembly))

checkmdfFiltered <- read.delim(checkMdf, header = T) 
rownamesDF <- rownames(checkmdfFiltered)
colnamesDF <- colnames(checkmdfFiltered)

checkmdfFiltered <- checkmdfFiltered |>
  rename(Assembly = Bin.Id) |>
  filter(grepl(lineageSearch, Marker.lineage,ignore.case = T))

### Merging the dataframes
cat("Merging the dataframes\n")
mergeddf <- quastdf |> inner_join(checkmdfFiltered) |> mutate(Genome.Fraction = as.numeric(Genome.Fraction))

### Making some plots
cat("Making the summary plot\n")
mapvsMarker <- mergeddf |> ggplot(aes(x = Completeness/100, y = Genome.Fraction/100)) +
  theme_classic() +
  geom_point(alpha =0.25) +
  geom_smooth(method = "lm") +
  xlab("Marker Completeness") +
  ylab("Mapped Fraction")

FilterLengths <- mergeddf |> ggplot(aes(x = Total.length, y = Completeness/100)) +
  theme_classic() +
  geom_point(alpha =0.25) +
  geom_hline(yintercept = 90/100, lty = 2, colour = "red") +
  geom_vline(xintercept = 1e6, lty = 2, colour = "red") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  xlab("Assembly Length") +
  ylab("Marker Completeness")

histo <- mergeddf |> ggplot(aes(x = Completeness/100)) +
  theme_classic() +
  geom_histogram(bins = 100, fill = "black") +
  geom_vline(xintercept = 90/100, lty = 2, colour = "red") +
  ylab("Assemblies")+
  xlab("Marker Completeness")

plotOut <- ggarrange(FilterLengths, ggarrange(mapvsMarker, histo, labels = c("B","C")), ncol = 1, labels = c("A", ""))
ggsave(plotOut, file = paste0(output,"_AssemblyFiltering.pdf"), width = 9, height = 6)

# Filtering out the incomplete assemblies
cat("Performing the filtering\n")
mergeddf <- mergeddf |> filter(Genome.Fraction >= 90)

# Now we want to sort out who needs contigs to be filtered
dfLarge <- mergeddf |> filter(Total.length >= 1.2e6, Contamination <= 1)
dfNorm <- mergeddf |> filter(Total.length <= 1.2e6, Total.length >= 1e6, Contamination <= 1)

write.table(dfLarge$Assembly, file = paste0(output,"_Assemblies_Large.list"), col.names = F, row.names = F, quote = F)
write.table(dfNorm$Assembly, file = paste0(output,"_Assemblies.list"), col.names = F, row.names = F, quote = F)
