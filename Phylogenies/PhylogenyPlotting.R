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
library(svglite)

# For the colours
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')

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

TpalLintoSub <- function(x){
	if(x == 1){
		return("SS14")
	}else if(x == 5){
		return("TEN")
	}else if(x == 4){
		return("TPE")
#	}else if(x == c(3)){
#		return("Nichols - Madagascar")
	}else if(x %in% c(2,3)){
		return("Nichols")
	}else{
		stop("I don't know this one!")
	}
}

customLegendLineages <- function(lineages, startColour, stepsize, numRows = 3){
	# If no step size
	if(missing(stepsize)){
		stepsize <- 1/length(lineages)
	}

	if(length(lineages) == 1){
		linColours <- startColour
	}else if(length(lineages) %% 2){ # If an odd number
		middle <- median(lineages)
		indSlice <- 1:(middle - 1)
		linColours <- c(lighten(startColour, stepsize * rev(indSlice)), startColour, darken(startColour, stepsize * indSlice))

	}else{ # If an even number
		indSlice <- 1:(length(lineages)/2)
		linColours <- c(lighten(startColour, stepsize * rev(indSlice)), darken(startColour, stepsize * indSlice))
	}

	# Now to make the legend
	names(linColours) <- lineages
	customLegend <- legendGrob(names(linColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = linColours, pch = 1, cex = 2/3), nrow = numRows)

	return(list("Colours" = linColours, "Legend" = customLegend))

}

### Taking care of the metadata ####
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

poppunkLineage <- read.csv("Data/RevisedAllMask_lineages.csv") |>
  dplyr:::rename(Genome = id) |> mutate(Genome = sapply(Genome, function(x){
	  tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	  tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
  	  tmp <- gsub("_NA","",tmp)
	  return(tmp)
  })) |> distinct()

# Finally, we'll need the incomplete MLST results so that we can include 
# comparative notes
STData <- read.delim("Data/updatedMLST.tab") |>
  mutate(ST = gsub("[^[:alnum:] ]","", ST)) |>
  mutate(Sample = gsub("\\.\\d.*","", Sample)) |>
  rename(Genome = Sample)

NFdata <- STData |> filter(ST == "NF") |> apply(MARGIN = 1, function(x){
					   if(any(x[3:5] %in% "-")){
						   x[2] = "Failed"
					   }else{
						   x[2] = "Unknown"
					   }
					   return(x)
			       }) |> t() |> as.data.frame()
STData <- STData |> filter(ST != "NF") |> select(Genome,ST, contains("TP")) |> bind_rows(NFdata |> select(Genome, ST, contains("TP")) )

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

fullMeta <- metadata |> full_join(poppunkLineage) |> left_join(STData) |> left_join(subInfo) |>
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
 
# Loading in the info from Lieberman et al. 2021
tmp <- fullMeta |> group_by(Rank_1_Lineage) |> summarize(min(collection_date, na.rm = T), max(collection_date, na.rm = T))

# Dealing with the colours now
conversionList = c("SS14", "Nichols", "Nichols - Madagascar", "TPE", "TEN")
coloursClusters <- c("#F0454B", "#2B8737", "#8A600D", "#0369AC", "#92278F")
#coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3") # Without Mad
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- conversionList

# Making a quick function to setup the colours much easier
tmp <- poppunkLineage |> filter(Rank_5_Lineage == 1) |> pull(Rank_1_Lineage) |> unique() |> sort()
customLegendSS14 <- customLegendLineages(tmp, coloursClusters[1])

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 2) |> pull(Rank_1_Lineage) |> unique() |> sort()
customLegendNichols <- customLegendLineages(tmp, coloursClusters[2])

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 3) |> pull(Rank_1_Lineage) |> unique() |> sort()
customLegendNicholsMad <- customLegendLineages(tmp, coloursClusters[3])

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 4) |> pull(Rank_1_Lineage) |> unique() |> sort()
customLegendTPE <- customLegendLineages(tmp, coloursClusters[4])

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 5) |> pull(Rank_1_Lineage) |> unique() |> sort()
customLegendTEN <- customLegendLineages(tmp, coloursClusters[5])


customLegendShape <- legendGrob(c("Human", "Non-Human Primate", "Unknown"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
			   pch = 15:17, gp = gpar(col = "black", pch = c(15:17), cex = 2/3), nrow = 3)

coloursRank1 <- c(customLegendSS14$Colours,customLegendNichols$Colours,customLegendNicholsMad$Colours,customLegendTPE$Colours,customLegendTEN$Colours)

savedLegendCustom <- guides(custom = guide_custom(customLegendSS14$Legend,title = "SS14", order = 1),custom = guide_custom(customLegendNichols$Legend,title = "Nichols", order = 2),
       custom = guide_custom(customLegendNicholsMad$Legend,title = "Nichols - Madagascar", order = 3),
       custom = guide_custom(customLegendTPE$Legend,title = "TPE", order = 4), custom = guide_custom(customLegendTEN$Legend,title = "TEN", order = 5),
       custom = guide_custom(customLegendShape, title= "Host", order = 6),
       colour = guide_none(), shape = guide_none()) 

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

summarizedClusters[,c("Rank_5_Lineage","Rank_1_Lineage", "SS14-like", "Nichols-like", "NotTPA", "Unknown", "Failed", "Total")] |> write.csv("PopPUNKvsGrillovaRank1.csv", row.names = F, quote = F)

	#xtable(auto = T, row.names = F) |> print(file = "PopPUNKvsGrillovaRank1.tex")

# Making a MLST Profile table
knownSTs <- fullMeta |> filter(!(ST %in% c("Failed", "Unknown"))) |>
	mutate(across(contains("TP0"), ~ gsub(pattern = "\\*|\\?", replacement = "", x = .x)))|>
	select(Genome,ST, Rank_5_Lineage, Rank_1_Lineage, contains("TP0")) |>
	unite(Profile, contains("TP"), sep = ".") |> arrange(Profile)

unknownSTs <- fullMeta |> filter(!grepl("Reference|037762225|draft", Genome)) |> filter(ST %in% c("Failed", "Unknown")) |>
	mutate(across(contains("TP0"), ~ gsub(pattern = ".*\\*|.*\\?", replacement = "Y", x = .x))) |>
	mutate(across(contains("TP0"), ~ gsub(pattern = "-", replacement = "X", x = .x))) |>
	select(Genome,ST, Rank_5_Lineage, Rank_1_Lineage, contains("TP0")) |>
	unite(Profile, contains("TP"), sep = ".") |> arrange(Profile)
	
knownSTs |> bind_rows(unknownSTs) |> arrange(Rank_5_Lineage, Rank_1_Lineage, Profile) |> write.csv(file = "LineagetoProfile.csv")

tmp <- unknownSTs |> select(-Genome) |> group_by(Rank_5_Lineage, Rank_1_Lineage) |>
	arrange(Profile) |>
	reframe(NovelProfile = paste(unique(Profile)[!grepl("X", unique(Profile))], collapse = ";"),
	PartialProfile = paste(unique(Profile)[grepl("X", unique(Profile))], collapse = ";"))

knownSTs |> select(-Genome) |> group_by(Rank_5_Lineage, Rank_1_Lineage) |>	
	summarize(Profile = paste(mixedsort(unique(Profile)), collapse = ";")) |> full_join(tmp) |>
	write.csv("LineagetoProfileCompiled.csv")

#### The CSV Tables ####
# Understanding the lineages by composition # NOTE NAs are from either the reference genome (included already) or the ancient genome (DIFFERENT PAPER)

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(geo_loc_name, name = "Count") |>
  pivot_wider(id_cols = geo_loc_name, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |> write.csv(file = "LineageGeography.csv", row.names = F)

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(sub_species, name = "Count") |>
  pivot_wider(id_cols = sub_species, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |>write.csv(file = "LineageSubspecies.csv", row.names = F)

fullMeta |> group_by(Rank_5_Lineage, Rank_1_Lineage) |> filter(!is.na(Rank_5_Lineage)) |> count(host, name = "Count") |> 
  pivot_wider(id_cols = host, names_from = Rank_1_Lineage, values_from = Count, values_fill = 0) |>write.csv(file = "LineageHost.csv", row.names = F)

#### Now to quickly do a SNP distances test ####
snpMat <- read.delim(file = "Data/snpMatrixRevisionNoOdd.tab.gz", header = F) |>
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
	filter(V1 != "Reference", V2 != "Reference")
	#mutate(V1 = gsub("Reference", "GCA_000017825", V1), V2 = gsub("Reference", "GCA_000017825", V2))

heatmapMatrix <- snpMat |> pivot_wider(id_cols = V1, names_from = V2, values_from = V3) |> as.data.frame()

rownames(heatmapMatrix) <- heatmapMatrix[,1]
heatmapMatrix <- heatmapMatrix[,-1]

heatmapTest <- fullMeta |> filter(Genome %in% rownames(heatmapMatrix)) |> select(Rank_5_Lineage, Rank_1_Lineage) |> 
	mutate(Rank_5_Lineage = conversionList[as.numeric(Rank_5_Lineage)]) |>
	mutate(SubLineage = as.factor(Rank_1_Lineage),Lineage = as.factor(Rank_5_Lineage), .keep = "none") |>
	select(-SubLineage)

test <- list(Lineage = coloursClusters, SubLineage = c(customLegendSS14$Colours, customLegendNichols$Colours,
						       customLegendNicholsMad$Colours, customLegendTPE$Colours, customLegendTEN$Colours))
names(test[[1]]) <- conversionList

pheatmap(heatmapMatrix, annotation_row = heatmapTest, clustering_method = "ward.D2", annotation_legend = T, border_color = NA,
	 annotation_col = heatmapTest, show_rownames = F, show_colnames = F, annotation_colors = test, filename = "Figures/SNPHeatmap.pdf",
	 annotation_names_col = F,
	 width = 6, height = 6)
dev.off()

# Now to do the tests *within* the groups
snpDF <- snpMat |> filter(V1 != V2) |> mutate(Start = heatmapTest[V1,"Lineage"], End = heatmapTest[V2,"Lineage"]) |> 
	select(Start, End, V3) |> mutate(Start = factor(as.character(Start), conversionList),
					 End = factor(as.character(End),conversionList)) |>
	as_tibble() |> rename(SNP = V3) 

summarizedSNPdf <- snpDF |> group_by(Start, End) |> summarize(Median = median(SNP), Mean = mean(SNP), SE = sd(SNP)/sqrt(length(SNP)))

summarizedSNPdf |> ggplot(aes(x = End, y = Mean, colour = Start, group = Start)) +
	theme_classic()+
	scale_colour_manual(values = coloursClusters, "") +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) +
	geom_line(show.legend = F) +
	geom_errorbar(aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), width  = 0.1,show.legend = F) +
	geom_point() +
	geom_text_repel(aes(label = round(Mean, 2)), show.legend = F) +
	xlab("Phylogenetic Clade") +
	ylab("Mean SNPs") +
	theme(legend.position = "bottom")
ggsave("Figures/InteractionPlot.pdf", width = 6, height = 4)

snpDF |> lm(formula = SNP ~ Start * End) |> summary()

snpDF |> filter(Start == End) |> 
	lm(formula = SNP ~ Start) |> emmeans(specs = "Start") |> pairs() |> confint()

snpDF |> filter(Start != End) |> unite(Start, End, col = "Comparison", sep = ":") |>
	filter(Comparison %in% c("Nichols:SS14", "Nichols:TPE", "Nichols:TEN", "Nichols:Nichols - Madagascar", "SS14:Nichols - Madagascar",
				 "SS14:TPE", "SS14:TEN", "Nichols - Madagascar:TPE","Nichols - Madagascar:TEN", "TPE:TEN")) |>
#snpDF |> filter(Start != End) |> 
	lm(formula = SNP ~ Comparison) |> emmeans(specs = "Comparison") |> as.data.frame()

#### Looking at the ANI Results ####
# Loading in the results
aniTable <- read.delim("Data/ANIResults.tab.gz", header = F, col.names = c("Query", "Reference", "ANI", "Matches", "Total")) |> as_tibble() |> 
	mutate(Query = basename(Query), Reference = basename(Reference), Query = gsub("\\.\\d_.*|.fasta","", Query), Reference = gsub("\\.\\d_.*|.fasta","", Reference)) |>
	select(-Matches, -Total) |> filter(!grepl("GCF", Query), !grepl("GCF", Reference))

# Preparing the heatmap
aniWide <- aniTable |> pivot_wider(Query, names_from = Reference, values_from = ANI)
aniMat <- aniWide[,-1] |> as.matrix()
rownames(aniMat) <- aniWide$Query

heatmapTest <- fullMeta |> filter(Genome %in% rownames(aniMat)) |> select(Rank_5_Lineage, Rank_1_Lineage) |> 
	mutate(Rank_5_Lineage = conversionList[as.numeric(Rank_5_Lineage)]) |>
	mutate(SubLineage = as.factor(Rank_1_Lineage),Lineage = as.factor(Rank_5_Lineage), .keep = "none") |>
	select(-SubLineage)
distanceMatrix <- as.dist(1-aniMat/100, upper = T, diag = F)

heatCol <- list(Lineage = coloursClusters, SubLineage = c(customLegendSS14$Colours, customLegendNichols$Colours,
						       customLegendNicholsMad$Colours, customLegendTPE$Colours, customLegendTEN$Colours))
names(heatCol[[1]]) <- conversionList

pheatmap(distanceMatrix, annotation_row = heatmapTest, clustering_method = "ward.D2", annotation_legend = T, border_color = NA,
	 annotation_col = heatmapTest, show_rownames = F, show_colnames = F, annotation_colors = heatCol,
	 annotation_names_col = F,
	 width = 6, height = 6)

# Let's run some quick stats like we did in the SNP analysis

aniDF <- aniTable |> mutate(Query = heatmapTest[Query, "Lineage"], Reference = heatmapTest[Reference, "Lineage"]) |>
	mutate(Query = factor(as.character(Query), conversionList),
					 Reference = factor(as.character(Reference),conversionList))

aniDF |> lm(formula = ANI ~ Query * Reference) |> summary()

aniDF |> filter(Query != Reference) |> unite(Query, Reference, col = "Comparison", sep = ":") |>
	filter(Comparison %in% c("Nichols:SS14", "Nichols:TPE", "Nichols:TEN", "Nichols:Nichols - Madagascar", "SS14:Nichols - Madagascar",
				 "SS14:TPE", "SS14:TEN", "Nichols - Madagascar:TPE","Nichols - Madagascar:TEN", "TPE:TEN")) |>
#snpDF |> filter(Query != Reference) |> 
	lm(formula = ANI ~ Comparison) |> emmeans(specs = "Comparison")  |> summary()
#### Reading in the full phylogeny ####
tree <- read.tree("Data/ConstantSitesWhole.treefile")
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

poppunkLin1 <- ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-6) +
  geom_treescale(linesize = 1, offset = 2) +
  #geom_nodepoint(pch = 23, colour = ifelse(bootstrapValues >= 75, "white", NA), fill = ifelse(bootstrapValues >= 90, "black", ifelse(bootstrapValues >= 75, "grey", NA)), size = 2) +
  #geom_text(aes(label = ifelse(node > 200, node, NA))) +
  geom_tippoint(aes(colour = as.factor(Rank_1_Lineage), shape = host)) +
  scale_colour_manual(values = coloursRank1, name = "Lineage") +
  savedLegendCustom +
  theme(legend.position = "bottom", legend.title.position = "top")

ggsave(poppunkLin1, file = "Figures/PhylowithLineage1Full.pdf", width = 9, height = 6)

#### Reading in the Trimmed Phylogeny ####
# Let's do the same thing, but, now with only trimmed!
tree <- read.tree("Data/ConstantSites.treefile")
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
  savedLegendCustom +
  scale_colour_manual(values = coloursRank1)  +
  theme(legend.position = "bottom")

tree$tip.label[!(tree$tip.label %in% rownames(fullMeta))]
bootstrapValues <- as.numeric(tree$node.label)
ggtree(tree, right = T) %<+% fullMeta +
  theme_tree() +
  geom_rootedge(1e-6) +
  geom_treescale(linesize = 1, offset = 2) +
 # geom_nodepoint(pch = 23, colour = ifelse(bootstrapValues >= 75, "white", NA), fill = ifelse(bootstrapValues >= 90, "black", ifelse(bootstrapValues >= 75, "grey", NA)), size = 2) +
  #geom_text(aes(label = ifelse(node > 300, node, NA))) +
  geom_tippoint(aes(colour = as.factor(Rank_1_Lineage), shape = host)) +
  scale_colour_manual(values = coloursRank1, name = "Lineage") +
  savedLegendCustom +
  #guides(custom = guide_custom(customLegendShape, title = "Host", order = 5),colour = guide_legend(nrow = 3), shape = guide_none()) +
  theme(legend.position = "bottom", legend.title.position = "top") +
  geom_cladelab(node = 326, label = "TPA", align = T, textcolor = "black") +
  geom_cladelab(node = 327, label = "Nichols", align = F, textcolor = coloursClusters[2]) +
  geom_cladelab(node = 325, label = "SS14", align = F, textcolor  = coloursClusters[1]) +
  geom_cladelab(node = 394, label = "Nichols - Madagascar", align = F, textcolor  = coloursClusters[3]) +
  geom_cladelab(node = 402, label = "TPE", align = T, textcolor  = coloursClusters[4]) +
  geom_cladelab(node = 423, label = "TEN", align = T, textcolor = coloursClusters[5])
#ggsave(file = "Figures/Lineage1TrimmedBoostrap.pdf", width = 12, height = 9)

ggsave(file = "Figures/Lineage1TrimmedBoot.pdf", width = 12, height = 9)

ggarrange(poppunkLin5, poppunkLin4, poppunkLin3, poppunkLin2, labels = "auto")

ggsave("Figures/PhylowithLineagesTrimmed.pdf", width = 12, height = 9)

as_data_frame(tree) |> arrange(branch.length)|> filter(nchar(label) > 3) |> pull(branch.length) |> log10() |> hist()

as_data_frame(tree) |> arrange(branch.length) |> filter(branch.length < 6e-4, nchar(label) > 3) |> pull(label) |> write.table("FilteredTips.list", quote = F, row.names = F)

as_data_frame(tree) |> arrange(-branch.length) |> filter(branch.length >= 6e-4, nchar(label) > 3) |> write.table("RemovedTips.list", quote = F, row.names = F)

# Getting the metadata sorted in order of the tree for the supplemental data
tree2 <- ladderize(tree, F)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
orderedTips <- tree2$edge[is_tip, 2]

fullMeta[tree$tip.label[orderedTips],] |> write.table(file = "OrderedMetadata.tab", sep = "\t", row.names = F, quote = F)

#### Preparing the table showing when each strain was used ####
whenUsed <- fullMeta |> select(Genome, Rank_5_Lineage) |> mutate(Lineage = conversionList[Rank_5_Lineage]) |> 
	select(-Rank_5_Lineage) |>
	mutate(Pangenome = ifelse(grepl("^G", Genome), T, F), PopPUNK = T, CoreSNPs = T, TrimmedPhylo = Genome %in% Tiplabels) |>
	as_tibble() 
whenUsed |> write.csv("SuppUsed.csv", row.names = F, quote = F)

whenUsed |> select(-Genome) |> pivot_longer(-Lineage, names_to = "Analysis", values_to = "PA") |>
	group_by(Analysis, Lineage) |> summarize(Genomes = sum(PA)) |>
	pivot_wider(id_cols = Lineage, names_from = Analysis, values_from = Genomes) |>
	write.csv("SuppUsedSummarized.csv", row.names = F, quote = F)
