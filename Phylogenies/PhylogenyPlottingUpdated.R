library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)
library(ggpubr)
library(ggstance)
library(ggnewscale)
library(ggrepel)
library(phytools)
library(ape)
library(gtools)
library(lubridate)
library(latex2exp)
library(pheatmap)
library(xtable)
library(grid)
library(colorspace)
library(emmeans)
library(scales)

# For the colours
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')
#coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3")
#SS14Colours <- c(lighten(coloursClusters[1], 0.15 * 2:1), coloursClusters[1], darken(coloursClusters[1], 0.15 * 1:2))
#names(SS14Colours) <- c(1,3,4,13,20)
#NicholsColours <- c(lighten(coloursClusters[2], 0.1 * 5:1), coloursClusters[2], darken(coloursClusters[2], 0.1 * 1:5))
#names(NicholsColours) <- c(2,5:9,15,19,21,22,24)
#TPEColours <- c(lighten(coloursClusters[3], 0.15 * 3:1), coloursClusters[3], darken(coloursClusters[3], 0.15 * 1:3))
#names(TPEColours) <- c(10,11,12,14,17,23,25)
#TENColours <- c(lighten(coloursClusters[4], 0.25), darken(coloursClusters[4], 0.25))
#names(TENColours) <- c(16,18)
#
#coloursRank1 <- c(SS14Colours, NicholsColours, TPEColours, TENColours)

# Functions
metadataParsing <- function(file){
  metadata <- read.delim(file, header = T)[-1,] 
  colnames(metadata) <- c("Genome", "name", "value")
  
  metadata <- metadata |>
    filter(name %in% c("biotic_relationship", "collection_date", "description",
                       "geo_loc_name", "sub_species", "strain", "host", 
                       "host_disease", "host_health_state", "isolate",
                       "isolation_source", "source_type", "pathogenicity")) |> 
    pivot_wider(id_cols = Genome,values_from = value,names_from = name, values_fill = NA,
                values_fn = function(x){paste(unlist(x),collapse = "|",sep = "|")})
 return(metadata)
}

fixingDates <- function(date){
  date <- gsub("[A-z]","", date) # Making sure it's only numbers
  lengthString <- nchar(date)
  lengthString <- ifelse(is.na(lengthString), 0, lengthString) # In case it's NA
  if(lengthString == 10){ # if it's the full amount!
    date <- gsub("[^[:alnum:] ]","-", date) # Making sure it's POSIX
    return(date)
  }else if(lengthString == 7){ # If year-month
    date <- gsub("[^[:alnum:] ]","-", date) # Making sure it's POSIX
    date <- paste0(date, "-01")
    return(date) 
  }else if(lengthString == 4){ # If year-
    date <- paste0(date, "-01-01")
    return(date)
  }else{ # If it's none of the above
    return(NA)
  }
}

# Taking care of the metadata
#ann_colors <- list(ST = c("Outgroup" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Geridu" = colour[5],"Ancient" = colour[1], "NF" = colour[20], "NIPH" = colour[22]))
# Loading in the metadata
metadata <- metadataParsing("Data/T160D20240604OtherMetaData.tab")
diseaseToSub <- c("syp" = "pallidum", "con" = "pallidum", "yaw" = "pertenue",
                  "end" = "endemicum", "bej" = "endemicum", "unk" = "Unknown")

metadata[nrow(metadata) + 1,] <- metadata[2,]
metadata[nrow(metadata),1] <- "Reference"
# Sorting out the subspecies
metadata <- metadata |> mutate(host_disease = tolower(host_disease),host_disease = replace(host_disease,
                                                                                           is.na(host_disease) | host_disease == "not collected", "unknown" )) |>
  mutate(host_disease = replace(host_disease, grepl("^syph", host_disease), "syphilis")) |>
  mutate(sub_species = tolower(sub_species), host_disease = tolower(host_disease))
metadata$sub_species[is.na(metadata$sub_species)] <- diseaseToSub[substr(metadata$host_disease[is.na(metadata$sub_species)],1,3)] |> unname()

# Some other small things that need doing
metadata <- metadata |> mutate(geo_loc_name = gsub(":.*","", geo_loc_name),
                               Genome = gsub("\\.\\d$","", Genome),
                               geo_loc_name = replace(geo_loc_name, grepl("not", geo_loc_name), NA)) 
metadata$collection_date <- sapply(metadata$collection_date, fixingDates)
metadata <- metadata |> mutate(collection_date = decimal_date(as.Date(collection_date)))

# Now to load in the Poppunk data
#Poppunkclusters <- read.csv("Data/Assigned_clusters.csv") |>
#  rename(Genome = Taxon) |> mutate(Genome = sapply(Genome, function(x){
#	  tmp <- strsplit(x, split = "_") |> unlist() |> unname()
#	  tmp <- paste(tmp[1:(length(tmp)-1)], sep = "_", collapse = "_" )
#	  return(tmp)
#  }))

PoppunkLineages <- read.csv("Data/AllMasked_lineages.csv") |>
#PoppunkLineages <- read.csv("Data/TomBealeMapped_lineages.csv") |>
  dplyr:::rename(Genome = id) |> mutate(Genome = sapply(Genome, function(x){
	  tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	  tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
  	  tmp <- gsub("_NA","",tmp)
	  return(tmp)
  })) |> distinct()

# Finally, we'll need the incomplete MLST results so that we can include 
# comparative notes
STData <- read.delim("Data/MLSTResults.txt") |>
	bind_rows(read.delim("Data/BealeMLST.txt")) |> 
  mutate(ST = gsub("[^[:alnum:] ]","", ST)) |>
  mutate(Sample = gsub("\\.\\d.*","", Sample)) |>
  rename(Genome = Sample)

NFdata <- STData |> filter(ST == "NF") |> apply(MARGIN = 1, function(x){
					   if(any(x[3:4] %in% "-")){
						   x[2] = "Failed"
					   }else{
						   x[2] = "Unknown"
					   }
					   return(x)
			       }) |> t() |> as.data.frame()
STData <- STData |> filter(ST != "NF") |> select(Genome,ST, contains("TP")) |> bind_rows(NFdata |> select(Genome, ST, contains("TP")) )
  #mutate(ST = replace(ST, ST == "NF" & !any(c(TP0136 == "-",TP0548 == "-", TP0705 == "-")), "Novel"))

# Need to include SS14/Nichols information so we can distinguish between the two!
subInfo <- read.delim("Data/TpallidumProfile.tab")[,c(1,5)] |> mutate(ST = as.character(ST))

bealeMetadata <- read.csv("../MetadataExtra/Beale2021.csv") |> 
  dplyr::rename(Genome = SRR.ENA_Accession, geo_loc_name = Geo_Country,
         collection_date = Sample_Year) |>
  mutate(collection_date = as.numeric(collection_date))

metadata <- metadata |> filter(!grepl("Reference", Genome)) |>
  bind_rows(bealeMetadata) |>
  filter(!grepl("GCA_037762225|ERR4853567|ERR7123579|ERR7123580|ERR7123581|ERR7123582|ERR7123585|ERR7123588|ERR7123590|ERR7123591", Genome)) |>
  unite(collection_date, contains("collection"),na.rm = T) |>
  unite(geo_loc_name, contains("geo"),na.rm = T) |>
  mutate(geo_loc_name = replace(geo_loc_name, nchar(geo_loc_name) == 0, "Other")) |>
  mutate(collection_date = as.numeric(collection_date))

fullMeta <- metadata |> full_join(PoppunkLineages) |> left_join(STData) |> left_join(subInfo) |>
  mutate(clonal_complex = replace(clonal_complex, ST == "Failed", "Failed")) |>
  mutate(clonal_complex = replace(clonal_complex, ST == "Unknown", "Unknown")) |>
  mutate(clonal_complex = replace(clonal_complex, nchar(clonal_complex) == 0, "NotTPA")) |>
  rename(TPASub = clonal_complex) |> filter(!is.na(Genome),!is.na(Rank_1_Lineage)) |> distinct() |>
  as.data.frame() #|>

rownames(fullMeta) <- fullMeta$Genome

# Want to do some filtering so that the geography is less confusing
smallCountries <- fullMeta |> count(geo_loc_name, name = "Count") |> filter(Count < 10) |>
  pull(geo_loc_name)
fullMeta <- fullMeta |>
  mutate(geo_loc_name = replace(geo_loc_name,
                                geo_loc_name %in% smallCountries | is.na(geo_loc_name),
                                "Other"))
# IDing the Canadian Samples and simplifying the host information
fullMeta <- fullMeta |>
  mutate(Canada = ifelse(grepl("_draft$",Genome), "Canadian", "Global")) |>
  mutate(Canada = replace(Canada, is.na(Canada), "Global")) |>
  mutate(host = replace(host, grepl("^ERR|SRR", Genome), "Homo sapiens")) |>
  mutate(host = ifelse(host == "Homo sapiens", "Human", ifelse(host == "NA"| is.na(host), "Unknown", "Primate"))) |>  #|>
  mutate(host = replace(host, is.na(host), "Unknown")) #|>
rownames(fullMeta) <- gsub("_NA","",rownames(fullMeta))
 
#fullMeta[nrow(fullMeta),4:8] <- 1

# Loading in the info from Lieberman et al. 2021
tmp <- fullMeta |> group_by(Rank_1_Lineage) |> summarize(min(collection_date, na.rm = T), max(collection_date, na.rm = T))

# Dealing with the colours
coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3")
conversionList = c("SS14", "Nichols", "TPE", "TEN")
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols","TPE","TEN")

tmp <- PoppunkLineages |> filter(Rank_5_Lineage == 1) |> pull(Rank_1_Lineage) |> unique() |> sort()
SS14Colours <- c(lighten(coloursClusters[1], 0.15 * 2:1), darken(coloursClusters[1], 0.15 * 1:2))
names(SS14Colours) <- tmp
customLegendSS14 <- legendGrob(names(SS14Colours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = SS14Colours, pch = 1, cex = 2/3), nrow = 3)

tmp <- PoppunkLineages |> filter(Rank_5_Lineage == 2) |> pull(Rank_1_Lineage) |> unique() |> sort()
NicholsColours <- c(lighten(coloursClusters[2], 0.1 * 4:1), coloursClusters[2], darken(coloursClusters[2], 0.1 * 1:4))
names(NicholsColours) <- tmp
customLegendNichols <- legendGrob(names(NicholsColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = NicholsColours, pch = 1, cex = 2/3), nrow = 3)

tmp <- PoppunkLineages |> filter(Rank_5_Lineage == 3) |> pull(Rank_1_Lineage) |> unique() |> sort()
TPEColours <- c(lighten(coloursClusters[3], 0.15 * 3:1), coloursClusters[3], darken(coloursClusters[3], 0.15 * 1:3))
names(TPEColours) <- tmp
customLegendTPE <- legendGrob(names(TPEColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = TPEColours, pch = 1, cex = 2/3), nrow = 3)

tmp <- PoppunkLineages |> filter(Rank_5_Lineage == 4) |> pull(Rank_1_Lineage) |> unique() |> sort()
TENColours <- c(lighten(coloursClusters[4], 0.25), darken(coloursClusters[4], 0.25))
names(TENColours) <- tmp
customLegendTEN <- legendGrob(names(TENColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = TENColours, pch = 1, cex = 2/3), nrow = 3)
coloursRank1 <- c(SS14Colours, NicholsColours, TPEColours, TENColours)

# Making sure the numbers add up everywhere
genomesUsed <- read.table("../AMR/GenomeList.list", col.names = F) |> unlist() |> unname() # This is making sure I'm only including the genomes of interest
fullMeta <- fullMeta |> filter(Genome %in% genomesUsed) 

#full_join(Poppunkclusters) |> 
# Making for simple clustering
#sortedClusters <- fullMeta$Cluster |> unique() |> mixedsort()

# Going to simplify the clusters so that only major groups pop-up
#singletons <- fullMeta |> count(Cluster,name = "Count") |>
#  filter(Count == 1) |> pull(Cluster)
#
#fullMeta <- fullMeta |>
#  mutate(Cluster = replace(Cluster, Cluster %in% singletons, "Singletons"))
#### Comparing the Poppunk Lineages and the MLST Results ####
# Now we need to show that the poppunk clusters are consistent!
summarizedClusters <- fullMeta |> filter(!grepl("Reference|037762225", Genome)) |> group_by(Rank_5_Lineage,Rank_1_Lineage) |> count(TPASub, name = "Genomes") |>
  #mutate(Genomes = Genomes/sum(Genomes)) |>
  pivot_wider(id_cols = c(Rank_5_Lineage, Rank_1_Lineage), names_from = TPASub, values_from = Genomes, values_fill = 0) |>
  as.data.frame()

summarizedClusters$Total <- rowSums(summarizedClusters[,-c(1:2)])
tmpList <- colSums(summarizedClusters[,-c(1,2)]) |> unlist() |> as.list()
tmpList$Rank_1_Lineage <- NA
tmpList$Rank_5_Lineage <- "Total"
summarizedClusters[nrow(summarizedClusters) + 1,] <- tmpList[colnames(summarizedClusters)]
rm(tmpList)

summarizedClusters[,c("Rank_5_Lineage","Rank_1_Lineage", "SS14-like", "Nichols-like", "NotTPA", "Unknown", "Failed", "Total")] |> xtable(auto = T, row.names = F) |> print(file = "PopPUNKvsGrillovaRank1.tex")

# Making that table
tmp <- fullMeta |> filter(!grepl("Reference|037762225|draft", Genome)) |>
	filter(ST != "Failed") 

unknownSTs <- tmp |> filter(ST == "Unknown") |>  select(ST, contains("TP0")) |> distinct() |> arrange(TP0136, TP0548, TP0705)
unknownSTs$ST <- paste0("t",1:nrow(unknownSTs))

tmp |> left_join(unknownSTs, by = c("TP0136","TP0548", "TP0705")) |> 
	mutate(ST.x = replace(ST.x, ST.x == "Unknown", NA)) |> 
	unite(ST, ST.x, ST.y, na.rm = T) |>
	group_by(Rank_5_Lineage, Rank_1_Lineage) |> distinct(ST) |> 
	summarize(ST = paste(mixedsort(ST), collapse = ",")) |>
	arrange(Rank_5_Lineage, Rank_1_Lineage) |> xtable(auto = T, row.names = F) |> print(file = "LineagetoST.tex")

#### The Latex Tables ####
# Understanding the lineages by composition # NOTE NAs are from either the reference genome (included already) or the ancient genome (DIFFERENT PAPER)

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(geo_loc_name, name = "Count") |>
  pivot_wider(id_cols = geo_loc_name, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |>
  xtable(auto = T) |> print(file = "LineageGeography.tex")

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(sub_species, name = "Count") |>
  pivot_wider(id_cols = sub_species, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |>
  xtable(auto = T) |> print(file = "LineageSubspecies.tex")

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(host, name = "Count") |> 
  pivot_wider(id_cols = host, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |>
  xtable(auto = T) |> print(file = "LineageHost.tex")

#### Now to quickly do a SNP distances test ####
snpMat <- read.delim(file = "Data/snpDistances.tab", header = F) |>
	mutate(V1 = sapply(V1, function(x){
	tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	return(tmp)
				}), V1 = gsub("\\.\\d|_NA","",V1)) |>
	mutate(V2 = sapply(V2, function(x){
	tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	return(tmp)
				}), V2 = gsub("\\.\\d|_NA","",V2)) |>
	mutate(V1 = gsub("Reference", "GCA_000017825", V1), V2 = gsub("Reference", "GCA_000017825", V2))

heatmapMatrix <- snpMat |> pivot_wider(id_cols = V1, names_from = V2, values_from = V3) |> as.data.frame()

rownames(heatmapMatrix) <- heatmapMatrix[,1]
heatmapMatrix <- heatmapMatrix[,-1]

heatmapTest <- fullMeta |> filter(Genome %in% rownames(heatmapMatrix)) |> select(Rank_5_Lineage, Rank_1_Lineage) |> 
	mutate(SubLineage = as.factor(Rank_1_Lineage),Lineage = as.factor(Rank_5_Lineage), .keep = "none")

test <- list(Lineage = coloursClusters[1:4], SubLineage = c(SS14Colours, NicholsColours, TPEColours, TENColours))
names(test[[1]]) <- 1:4

pheatmap(heatmapMatrix, annotation_row = heatmapTest, clustering_method = "ward.D2", annotation_legend = F, border_color = NA,
	 annotation_col = heatmapTest, show_rownames = F, show_colnames = F, annotation_colors = test, filename = "Figures/SNPHeatmap.pdf",
	 annotation_names_col = F,
	 width = 9, height = 9)
dev.off()

# Now to do the tests *within* the groups
snpDF <- snpMat |> mutate(Start = heatmapTest[V1,"Lineage"], End = heatmapTest[V2,"Lineage"]) |> 
	select(Start, End, V3) |> as_tibble() |> mutate(Start = conversionList[as.numeric(Start)], Start = factor(Start, levels = c("SS14", "Nichols", "TPE", "TEN"))) |> 
	mutate(End = conversionList[as.numeric(End)], End = factor(End, levels = c("SS14", "Nichols", "TPE", "TEN"))) |> rename(SNP = V3)

summarizedSNPdf <- snpDF |> group_by(Start, End) |> summarize(Median = median(SNP), Mean = mean(SNP), SE = sd(SNP)/sqrt(length(SNP)))

summarizedSNPdf |> ggplot(aes(x = End, y = Mean, colour = Start, group = Start)) +
	theme_classic()+
	scale_colour_manual(values = coloursClusters, "Clade") +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) +
	geom_line(show.legend = F) +
	geom_errorbar(aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), width  = 0.1,show.legend = F) +
	geom_point() +
	geom_text_repel(aes(label = round(Mean, 2)), show.legend = F) +
	xlab("Clade") +
	theme(legend.position = "bottom")
ggsave("Figures/InteractionPlot.pdf", width = 6, height = 4)


snpDF |> lm(formula = SNP ~ Start * End) |> summary()

snpDF |> filter(Start == End) |> select(Start, SNP) |>
	lm(formula = SNP ~ Start) |> emmeans(specs = "Start") |> pairs() 

#pairsData <- snpDF |> filter(Start == End) |> select(Start, SNP) |>
#	lm(formula = SNP ~ Start) |> emmeans(specs = "Start") |> pairs() |> as.data.frame() |>
#	mutate(contrast = gsub("Start","", contrast))

#pairsData |> ggplot(aes(y = contrast, x = estimate)) +
#		theme_classic() +
#		geom_vline(xintercept = 0, lty = 2) +
#  		geom_tile(aes(width = 1.96*SE*2, height = 0.5), fill = coloursClusters[3], alpha = 0.5) +
#		geom_point() +
#		ylab("Comparisons") +
#		xlab("Mean within clade SNPs")
#ggsave("Figures/SNPemmeans.pdf", width = 6, height = 4)
#### Reading in the full phylogeny ####
tree <- read.tree("Data/AllMaskedNoTrimming.nwk")
Tiplabels <- tree$tip.label

# Need to strip the dates off the ends of the tips for this tree
Tiplabels <- sapply(Tiplabels, function(x){
	tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	return(tmp)
	})

Tiplabels[which(grepl("Ref",Tiplabels))] <- "Reference"
Tiplabels <- gsub("\\.\\d|_NA","",Tiplabels)
Tiplabels[which(Tiplabels == "Reference")] <- "GCA_000017825"
tree$tip.label <- Tiplabels
tree <- midpoint(tree)

# Let's make the actual tree now
#poppunkTree <- ggtree(tree, right = T) %<+% fullMeta +
#  theme_tree() +
#  geom_rootedge(1e-3) +
#  geom_treescale(linesize = 1, offset = 2) +
#  geom_tippoint(aes(colour = Cluster)) +
#  scale_colour_manual(values = colour, name = "Poppunk Cluster") +
#  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
#  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))
#poppunkTree
#
#ggsave("TrimmedMLPhylo.png", width = 6, height = 4)
#
# Now to plot the lineages
poppunkLin5 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_5_Lineage))) +
  scale_colour_manual(values = coloursClusters, name = "Rank 5") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin4 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_4_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 4") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin3 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_3_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 3") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin2 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_2_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 2") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))
#ggsave(poppunkLin2, file = "Lineage2.png", width = 9, height = 6)

ggarrange(poppunkLin5, poppunkLin4, poppunkLin3, poppunkLin2, labels = "auto")

ggsave("Figures/PhylowithLineagesFull.pdf", width = 12, height = 9)

#### Reading in the Trimmed Phylogeny ####
# Let's do the same thing, but, now with only trimmed!
tree <- read.tree("Data/TrimmedBealeNoDraft.treefile")
Tiplabels <- tree$tip.label

# Need to strip the dates off the ends of the tips for this tree
Tiplabels <- sapply(Tiplabels, function(x){
	tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	return(tmp)
	})

Tiplabels[which(grepl("Ref",Tiplabels))] <- "Reference"
Tiplabels <- gsub("\\.\\d|_NA","",Tiplabels)
Tiplabels[which(Tiplabels == "Reference")] <- "GCA_000017825"
tree$tip.label <- Tiplabels
tree <- midpoint(tree)

customLegendShape <- legendGrob(c("Human", "Non-Human Primate", "Unknown"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), pch = c(16,17,15), gp = gpar(col = "black", cex = 2/3), nrow = 4)

# Now to plot the lineages
poppunkLin5 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_5_Lineage))) +
  scale_colour_manual(values = coloursClusters, name = "Rank 5") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin4 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_4_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 4") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin3 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_3_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 3") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin2 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_2_Lineage))) +
  scale_colour_manual(values = colour, name = "Rank 2") +
  guides(colour = guide_legend(nrow = 3, title.vjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(0.3,0.3), legend.background = element_rect(colour = "black"))

poppunkLin1 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  geom_tippoint(aes(colour = as.factor(Rank_1_Lineage))) +
  guides(custom = guide_custom(customLegendSS14,title = "SS14 - Lineage 1", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols - Lineage 2", order = 2),
	 custom = guide_custom(customLegendTPE,title = "TPE - Lineage 3", order = 3),custom = guide_custom(customLegendTEN,title = "TEN - Lineage 4", order = 4),
	 colour = guide_none()) +
  scale_colour_manual(values = coloursRank1)  +
  theme(legend.position = "bottom")

tree$tip.label[!(tree$tip.label %in% rownames(fullMeta))]

ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-4) +
  geom_treescale(linesize = 1, offset = 2) +
  #geom_nodelab(geom = "label") +
  #geom_text(aes(label = ifelse(node > 200, node, NA))) +
  geom_tippoint(aes(colour = as.factor(Rank_1_Lineage), shape = host)) +
  scale_colour_manual(values = coloursRank1, name = "Lineage") +
  guides(custom = guide_custom(customLegendSS14,title = "SS14 Sub-lineages", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols Sub-lineages", order = 2),
	 custom = guide_custom(customLegendTPE,title = "TPE Sub-lineages", order = 3),custom = guide_custom(customLegendTEN,title = "TEN Sub-lineages", order = 4),
	 custom = guide_custom(customLegendShape, title = "Host", order = 5), shape = guide_none(), colour = guide_none()) +
  #guides(custom = guide_custom(customLegendShape, title = "Host", order = 5),colour = guide_legend(nrow = 3), shape = guide_none()) +
  theme(legend.position = "bottom", legend.title.position = "top") +
  geom_cladelab(node = 161, label = "Nichols", align = T, aes(color = colour[2])) +
  geom_cladelab(node = 159, label = "SS14", align = T, aes(color = colour[1])) +
  geom_cladelab(node = 219, label = "TPE", align = T, aes(color = colour[8])) +
  geom_cladelab(node = 309, label = "TEN", align = T, aes(color = colour[12]))
#ggsave(file = "Figures/Lineage1TrimmedBoostrap.pdf", width = 12, height = 9)

ggsave(file = "Figures/Lineage1Trimmed.pdf", width = 9, height = 6)

ggarrange(poppunkLin5, poppunkLin4, poppunkLin3, poppunkLin2, labels = "auto")

ggsave("Figures/PhylowithLineagesTrimmed.pdf", width = 12, height = 9)

as_data_frame(tree) |> arrange(branch.length) |> filter(branch.length > 1e-6, branch.length < 1e-3) |> filter(nchar(label) > 3) |> pull(label) |> write.table("FilteredTips.list", quote = F, row.names = F)
