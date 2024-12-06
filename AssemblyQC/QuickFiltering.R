# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# The script
df <- read.delim("Data/Mapped/transposed_report.tsv") |> as_tibble()
colnames(df) <- gsub("X..","",colnames(df))

df <- df |># select(Assembly,contigs, Total.length, N50, Genome.fraction....) |>
  rename(Genome.Fraction = Genome.fraction....) |>
  mutate(Assembly = gsub("_contigs", "", Assembly))

checkmdfFiltered <- read.delim("Data/Mapped/SummarizedOutput.tab", header = T) 
rownamesDF <- rownames(checkmdfFiltered)
colnamesDF <- colnames(checkmdfFiltered)
#checkmdfFiltered$Strain.heterogeneity <- rownamesDF
#checkmdfFiltered <- checkmdfFiltered[,c(14,1:13)]
#colnames(checkmdfFiltered) <- colnamesDF

checkmdfFiltered <- checkmdfFiltered |>
  rename(Assembly = Bin.Id) |>
  filter(grepl("spiro", Marker.lineage,ignore.case = T))

df <- df |> inner_join(checkmdfFiltered) |> mutate(Genome.Fraction = as.numeric(Genome.Fraction))
hist(df$Genome.Fraction)
hist(df$Total.length)
hist(df$N50)

mapvsMarker <- df |> ggplot(aes(x = Completeness/100, y = Genome.Fraction/100)) +
  theme_classic() +
  geom_point(alpha =0.25) +
  geom_smooth(method = "lm") +
  #geom_hline(yintercept = 90/100, lty = 2, colour = "red") +
  #scale_x_log10() +
  #annotation_logticks(sides = "b") +
  xlab("Marker Completeness") +
  ylab("Mapped Fraction")

FilterLengths <- df |> ggplot(aes(x = Total.length, y = Completeness/100)) +
  theme_classic() +
  geom_point(alpha =0.25) +
  geom_hline(yintercept = 90/100, lty = 2, colour = "red") +
  geom_vline(xintercept = 1e6, lty = 2, colour = "red") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  xlab("Assembly Length") +
  ylab("Marker Completeness")

histo <- df |> ggplot(aes(x = Completeness/100)) +
  theme_classic() +
  geom_histogram(bins = 100, fill = "black") +
  geom_vline(xintercept = 90/100, lty = 2, colour = "red") +
  ylab("Assemblies")+
  xlab("Marker Completeness")

df |> filter(Completeness >= 90) |> ggplot(aes(x = Completeness/100,y = Contamination/100)) +
  theme_classic() +
  geom_point(alpha =0.25) 

ggarrange(FilterLengths, ggarrange(mapvsMarker, histo, labels = c("B","C")), ncol = 1, labels = c("A", ""))
ggsave("FilteringSummaryTrimmed100.png", width = 9, height = 6)

# Filtering out the incomplete assemblies
df <- df |> filter(Genome.Fraction >= 90)

hist(df$Total.length)
hist(df$Genome.Fraction)
hist(df$N50)

df |> ggplot(aes(x = Contamination, y = Total.length)) +
  theme_classic() +
  geom_point(alpha = 0.25) +
  geom_vline(xintercept = 1, lty = 2, colour = "red") +
  geom_hline(yintercept = 1.2e6, lty = 2, colour = "red") +
  ylab("Assembly Length")+
  xlab("Marker Contaminataion")

# Now we want to sort out who needs contigs to be filtered
dfLarge <- df |> filter(Total.length >= 1.2e6)
dfNorm <- df |> filter(Total.length <= 1.2e6, Total.length >= 1e6, Contamination <= 1)
dfNorm |> filter(Contamination <= 1)

write.table(dfLarge$Assembly, file = "LargerAssemblies.list", col.names = F, row.names = F, quote = F)
write.table(dfNorm$Assembly, file = "NormalAssemblies.list", col.names = F, row.names = F, quote = F)
