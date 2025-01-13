library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(scales)
library(ggExtra)
library(ggnewscale)
library(latex2exp)
library(parallel)
library(pbapply)
library(purrr)
library(gtools)
library(xtable)
library(ggforce)
library(lubridate)
library(eulerr)
library(xtable)
library(grid)
library(colorspace)
library(emmeans)
###
#### Functions ####
BlastParsing <- function(blastFile, ncores){
	tmp <- read.delim(blastFile, header = F,
			  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) |> as_tibble()
	
	tmpList <- split(tmp, as.factor((tmp$Query)))
	
	op <- pboptions(type = "timer")
	FilteredAmbiguousGenes<- pblapply(cl = ncores,tmpList, function(x){ 
		       	filtResults <- x |> filter(Evalue == min(Evalue))
	
			if(nrow(filtResults) > 1){
				filtResults <- filtResults |> filter(Bitscore == max(Bitscore))
				if(nrow(filtResults) > 1){
					filtResults <- filtResults |> filter(PIdent == max(PIdent))
					}
					if(nrow(filtResults) > 1){
						filtResults <- filtResults |> filter(Mismatches == min(Mismatches))
						}
				}
			return(filtResults[1,])
	 }) |> bind_rows() |> filter(PIdent >= 90)
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

panHeap <- function(panDF, iter = NULL){
  # Doing some set up here
  roaryMatrix <- panDF |> select(-Gene) |> as.matrix.data.frame()
  N = 1:ncol(roaryMatrix)

  coefList <- NULL
  
  if(length(iter) == 0){
    iter = ncol(roaryMatrix)
    if(iter <= 1000){
      cat(paste0("Doing all possible permutations (",iter,")\n"))
    }else{
      iter <- 1000
      cat(paste0("Doing ",iter," permutations\n"))
    }
  }
  
  # The actual work
  for(i in 1:iter){
    newGenes <- NULL
    forMat <- roaryMatrix
    forMat <- forMat[,sample(1:ncol(forMat), ncol(forMat))]
    
    newGenes <- NULL
    for(j in 1:ncol(forMat)){
      foundGenes <- forMat[,j]
      newGenes <- c(newGenes,sum(foundGenes)) # How many new genes do we have?
      
      if(j == 1){
        totalGenes = newGenes
      }else{
        totalGenes = c(totalGenes, sum(newGenes))
      } 
      indtosetzero <- which(foundGenes == 1)
      forMat[indtosetzero,] <- 0
    }
    
    # Now to figure out how to calculate Heaps value
    coefficients <- coef(nls(totalGenes ~ k*N^y, start = list(k = 1, y = 0)))
    coefList[[i]] <- coefficients
  }  
  return(coefList)
  
}

TpalLintoSub <- function(x){
	if(x == 1){
		return("SS14")
	}else if(x == 4){
		return("TEN")
	}else if(x == 3){
		return("TPE")
	}else if(x %in% c(2)){
		return("Nichols")
	}else{
		stop("I don't know this one!")
	}
}

calculateEllipses <- function(dat){
	# Comes from https://stackoverflow.com/questions/30825387/getting-the-parameters-of-a-data-ellipse-produced-by-the-car-package-in-r
	ell.info <- cov.wt(dat) # Getting the covariance matrix
	eigen.info <- eigen(ell.info$cov) # Need to calculate the eigen vectors
	lengths <- sqrt(eigen.info$values * 2 * qf(0.95,2, nrow(dat) - 1)) # How long are the eigenvectors for a 95%

	if(is.nan(lengths[2])){
		warning("Not enough variablity for an ellipse. Preparing a line instead.\n b will be equal to 1e-6")
		lengths[2] <- 1e-6
	}
	
	# Angle calculation for the major axis
	angle <- atan2(eigen.info$vectors[1,2],eigen.info$vectors[2,2])

	angle = -angle
	df <- data.frame("Centre_x" = ell.info$center[1], "Centre_y" = ell.info$center[2], "a" = lengths[1], "b" = lengths[2], "angle" = angle)
	return(df)
}

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')
smallColours <-  c("#3BB273", "#4D9DE0","#7768AE","#E1BC29","#E15554", "#D1A7A0")
theme_set(theme_classic()) # saving myself some commands

##### Reading in the metadata #####
# Loading in the metadata
metadata <- read.delim("Data/T160D20240604OtherMetaData.tab", header = T)[-1,] 
colnames(metadata) <- c("Genome", "name", "value")

metadata <- metadata |>
  filter(name %in% c("biotic_relationship", "collection_date", "description",
                     "geo_loc_name", "sub_species", "strain", "host", 
                     "host_disease", "host_health_state", "isolate",
                     "isolation_source", "source_type", "pathogenicity")) |> 
  pivot_wider(id_cols = Genome,values_from = value,names_from = name, values_fill = NA,
              values_fn = function(x){paste(unlist(x),collapse = "|",sep = "|")})

diseaseToSub <- c("syp" = "pallidum", "con" = "pallidum", "yaw" = "pertenue",
                  "end" = "endemicum", "bej" = "endemicum", "unk" = "Unknown")

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

# Want to do some filtering so that the geography is less confusing
smallCountries <- metadata |> count(geo_loc_name, name = "Count") |> filter(Count < 10) |>
  pull(geo_loc_name)
metadata <- metadata |>
  mutate(geo_loc_name = replace(geo_loc_name,
                                geo_loc_name %in% smallCountries | is.na(geo_loc_name),
                                "Other"))
# AMR info
amrInfo <- read.delim("Data/AMRHits.tab") |> 
#amrInfo <- read.delim("Data/NewAssembled/AMRHits.tab") |> 
  select(Genome, Drug.Class) |>
  pivot_wider(id_cols = Genome, names_from = Drug.Class, 
              values_from = Drug.Class,
              values_fn = function(x){as.logical(length(x))},
              values_fill = F) 

# Including the Poppunk info
poppunkLineage <- read.csv("Data/AllMasked_lineages.csv") |>
  rename(Genome = id) |> mutate(Genome = sapply(Genome, function(x){
	  tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	  tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	  tmp <- gsub("_NA$","", tmp)
	  return(tmp)
  })) |> distinct()

# Loading in the genomes lengths
genomeStats <- read.delim("Data/PanSequenceStats.tab", header = T) |> as_tibble() |>
  rename(Genome = File) |> mutate(Genome = sapply(Genome, function(x){
	  tmp <- strsplit(x, split = "_") |> unlist() |> unname()
	  tmp <- paste(tmp[1:2], sep = "_", collapse = "_" )
	  tmp <- gsub("_NA$|\\.\\d$","", tmp)
	  return(tmp)
  })) |> distinct()

metadata <- metadata |> full_join(poppunkLineage) |> left_join(amrInfo) |> left_join(genomeStats) # Merging it all together

# Dealing with the colours now
coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3")
conversionList = c("SS14", "Nichols", "TPE", "TEN")
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols","TPE","TEN")

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 1) |> pull(Rank_1_Lineage) |> unique() |> sort()
SS14Colours <- c(lighten(coloursClusters[1], 0.15 * 2:1), darken(coloursClusters[1], 0.15 * 1:2))
names(SS14Colours) <- tmp
customLegendSS14 <- legendGrob(names(SS14Colours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = SS14Colours, pch = 1, cex = 2/3), nrow = 3)

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 2) |> pull(Rank_1_Lineage) |> unique() |> sort()
NicholsColours <- c(lighten(coloursClusters[2], 0.1 * 4:1), coloursClusters[2], darken(coloursClusters[2], 0.1 * 1:4))
names(NicholsColours) <- tmp
customLegendNichols <- legendGrob(names(NicholsColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = NicholsColours, pch = 1, cex = 2/3), nrow = 3)

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 3) |> pull(Rank_1_Lineage) |> unique() |> sort()
TPEColours <- c(lighten(coloursClusters[3], 0.15 * 3:1), coloursClusters[3], darken(coloursClusters[3], 0.15 * 1:3))
names(TPEColours) <- tmp
customLegendTPE <- legendGrob(names(TPEColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = TPEColours, pch = 1, cex = 2/3), nrow = 3)

tmp <- poppunkLineage |> filter(Rank_5_Lineage == 4) |> pull(Rank_1_Lineage) |> unique() |> sort()
TENColours <- c(lighten(coloursClusters[4], 0.25), darken(coloursClusters[4], 0.25))
names(TENColours) <- tmp
customLegendTEN <- legendGrob(names(TENColours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = TENColours, pch = 1, cex = 2/3), nrow = 3)
coloursRank1 <- c(SS14Colours, NicholsColours, TPEColours, TENColours)

customLegendShape <- legendGrob(c("Human", "Non-Human Primate", "Unknown", "Lineage"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"),
				pch = c(16,17,15,18), gp = gpar(col = "black", cex = 2/3), nrow = 3)

coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols","TPE","TEN")
############# Determining if it's an open or closed pan-genome #####
roaryOutput <- as_tibble(read.delim("Data/gene_presence_absence.Rtab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("^X|ConfAncient|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\.", "-", colnames(roaryOutput)[2:ncol(roaryOutput)])

# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.95 * ncol(roaryOutput)) # Have incomplete assemblies, want to play it safe

# Now to figure out how to make the models!
# Likely want to do multiple permutations to show that it doesn't matter who's first

heapResults <- panHeap(roaryOutput) |> bind_rows()

heapSummarized <- heapResults |> 
  pivot_longer(everything(), values_to = "Value", names_to = "variable") |>
  group_by(variable) |>
  summarize(meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
                         meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
			 loValue = meanValue - 1.96 * seValue, hiValue = meanValue + 1.96 * seValue, "") |> as.data.frame()


heapSummarizedPlot <- heapSummarized |> mutate(variable = factor(variable, labels = c("K", TeX('$\\gamma$'))))

# Now to take a look at the individual sub-species
subSep <- pblapply(list.files("Data/IndivPan/", pattern = "*tab", full.names = T), function(x){
	sampleName <- gsub("\\.rtab","",basename(x))
	panGenome <- as_tibble(read.delim(x))
	heapResults <- panHeap(panGenome) |> bind_rows()

	heapSummarized <- heapResults |> 
	  pivot_longer(everything(), values_to = "Value", names_to = "variable") |>
	  group_by(variable) |>
  	  summarize(meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
                         meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
			 loValue = meanValue - 1.96 * seValue, hiValue = meanValue + 1.96 * seValue, Sample = sampleName, N = ncol(panGenome)) |> as.data.frame()
	return(heapSummarized)
			 }) |> bind_rows()

sepLin <- subSep |> bind_rows(heapSummarized |> mutate(Sample = "All", N = 544)) |> filter(Sample != "tpa") |> group_by(Sample, variable) |> 
	mutate(Sample = ifelse(Sample == "tpe", "TPE", ifelse(Sample == "ten", "TEN", ifelse(Sample == "ss14", "SS14", ifelse(Sample == "All", "All","Nichols"))))) |> 
	select(Sample, N, variable, meanValue, loValue, hiValue) |> group_by(Sample, variable) |>
	pivot_longer(-c(Sample, N, variable)) |> unite(New, variable, name) |> pivot_wider(id_cols = c(Sample, N), names_from = New, values_from = value) |>
	reframe(N = 1:N, y = k_meanValue*(N)^y_meanValue, ylo = k_loValue*(N)^y_loValue, yhi = k_hiValue*(N)^y_hiValue) |>
       	ggplot(aes(x = N, y = y, ymin = ylo, ymax = yhi, colour = Sample, fill = Sample)) +
	theme_classic() +
	geom_ribbon(alpha = 0.5, colour = NA) +
	scale_fill_manual(values = c(coloursClustersConversion, "All" = "#00B2E3")) +
	scale_colour_manual(values = c(coloursClustersConversion, "All" = "#00B2E3")) +
	geom_line() +
	xlab("Genomes") +
	ylab("Genes") +
	theme(legend.position = "bottom")

ggsave("Figures/HeapRegression.pdf",sepLin, width = 6, height= 4)

############# Core Gene Presence #####
accessGenes <- roaryOutput$Gene[!coreGenes]
coreGenes <- roaryOutput$Gene[coreGenes]

write.table(accessGenes, file = "AccessoryGenes.list", row.names = F, col.names = F, quote = F)
# Incorporating the ancient into the thresholds
roaryOutputwithAncient <- roaryOutput

pangenomeCounts <- roaryOutputwithAncient |> select(-Gene) |>
  summarize_all(sum) |>
  pivot_longer(everything(), names_to = "Genome", values_to = "Count") |>
  mutate(Genome = gsub("-.*","",Genome))
metadata <- metadata |> left_join(pangenomeCounts) |> filter(Genome != "GCA_037762225") |>
	mutate(SubspeciesSam = conversionList[Rank_5_Lineage]) |>
  	mutate(SubspeciesSam = factor(SubspeciesSam, levels = c("SS14","Nichols", "TPE", "TEN")))

geneBox <- metadata |> #filter(!grepl("ERR", Sample)) |> 
  ggplot(aes(y = Count, x = SubspeciesSam, colour = SubspeciesSam)) +
  geom_jitter(height = 0, width = 0.4) +
  geom_boxplot(outliers = F, fill = NA, colour = "black") + 
  scale_colour_manual(values = coloursClustersConversion, "Lineage") + 
  xlab("PopPUNK Lineage") +
  theme(legend.position = "bottom", legend.title.position = "top")  +
  ylab("Genes")

ggsave("Figures/BoxplotCounts.pdf",geneBox,  width = 9, height = 6)

# Now determining if genome length plays an important role

geneCountModel <- metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> lm(formula =  Count ~ SubspeciesSam + sum_len)

metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> lm(formula =  Count ~ SubspeciesSam + sum_len) |> emmeans(specs = "SubspeciesSam") #|> pairs()|> as.data.frame()

metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> ggplot(aes(x = sum_len, y = Count, colour = as.factor(Rank_1_Lineage), group = as.factor(Rank_5_Lineage))) +
	geom_point() +
	scale_colour_manual("Lineage", values = coloursRank1) +
	guides(custom = guide_custom(customLegendSS14,title = "SS14", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols", order = 2),
	       custom = guide_custom(customLegendTPE,title = "TPE", order = 3), custom = guide_custom(customLegendTEN,title = "TEN", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none()) +
	new_scale_colour() +
	geom_smooth(method = "lm", aes(colour = as.factor(Rank_5_Lineage)), show.legend = F) +
	scale_colour_manual(values= c("red", "green", "blue", "purple"))  +
	theme(legend.position = "bottom", legend.title.position = "top")

metadata |> group_by(Rank_5_Lineage) |> filter(!is.na(Count))|> summarize(mean(Count))
metadata |> filter(!is.na(Count))|> summarize(median(Count))

# Loading in the gene information
GenesID <- read.delim("Data/BestBlastHits.tab", header = F)[,c(1,2)]
colnames(GenesID) <- c("Gene", "Accession")
GenesID <- read.delim("Data/ProteinOutput.list") |> rename(Accession = ProteinID) |> 
	right_join(GenesID |> mutate(Accession = gsub("\\.\\d$","",Accession)), relationship = "many-to-many")

# Quickly looking to see if TPE has a larger accessory content
subSep <- pblapply(list.files("Data/IndivPan/", pattern = "*tab", full.names = T), function(x){
	sampleName <- gsub("\\.rtab","",basename(x))
	panGenome <- as_tibble(read.delim(x))
	
	coreGenes <- rowSums(panGenome[,-1]) >= floor(0.95 * ncol(panGenome)) # Have incomplete assemblies, want to play it safe
	
	core <- sum(coreGenes)
	access <- sum(!coreGenes)
	return(data.frame(Sample = sampleName, Core = core, Access = access))	
			 }) |> bind_rows()

##### Doing the Core Genome PCoA ####
coreDf <- roaryOutputwithAncient |> filter(Gene %in% coreGenes)

# Preparing the data
coreData <- as.data.frame(coreDf)
rownames(coreData) <- coreData[,1]
coreData <- coreData[,-1]
genomes <- colnames(coreData)
coreData <- apply(coreData, MARGIN = 1,FUN = as.numeric)
rownames(coreData) <- genomes

coreDataSub <- coreData
maxCount <- coreDataSub |> colSums() |> max()
ind <- coreDataSub |> colSums() == maxCount
coreDataSub <- coreDataSub[,!ind]

# If I want to use a MDS plot
geneDist <- dist(coreDataSub, method = "binary")

fit <- cmdscale(geneDist,eig = T, k = 4, add =T)
coord <- fit$points |> as_tibble()
coord <- coord[,1:4]
#coord <- coord[,c(fit$eig > 1)[1:100]]

coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100

coord <- coord |> #mutate(shape = replace(shape, grepl("draft", Genome), "Canadian")) |>
  mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata) |>
  select(Genome, host, sub_species, Rank_1_Lineage, Rank_5_Lineage, V1,V2, host) |>
  mutate(sub_species = factor(sub_species, levels = sort(unique(sub_species)))) |>
  mutate(Rank_1_Lineage = factor(Rank_1_Lineage, levels = sort(unique(Rank_1_Lineage)))) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage, levels = sort(unique(Rank_5_Lineage)))) |>
  mutate(host = ifelse(host == "Homo sapiens", "Human", ifelse(host == "NA"| is.na(host), "Unknown", "Non-Human Primate"))) |>  #|>
  mutate(host = replace(host, is.na(host), "Unknown")) #|>
  #left_join(amrInfo)

pCore <- coord |> #mutate(Assembled = ifelse(grepl("ERR", Genome), T, F)) |>
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	#geom_jitter(width = 0.15, height = 0.15,alpha = 0.75) +
  #geom_density_2d() +
	geom_point(alpha = 0.5) +
	scale_colour_manual("Lineage", values = colour) +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5),
	       shape = guide_legend(title = "Source",nrow = 2, title.position = "top", title.hjust = 0.5)) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	guides(custom = guide_custom(customLegendSS14,title = "SS14 - Lineage", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols - Lineage", order = 2),
	       custom = guide_custom(customLegendTPE,title = "TPE - Lineage", order = 3), custom = guide_custom(customLegendTEN,title = "TEN - Lineage", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none()) +
	theme(legend.position = "bottom", legend.title.position = "top") 

ggsave(pCore, file = "Figures/CorePCoA.pdf", width = 9, height = 6)

# Saving the core gene counts
subDF <- coreData |> as.data.frame()
subDF$Genome <- rownames(subDF)
geneCountsCore <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = sapply(Rank_5_Lineage,TpalLintoSub)) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes")

geneCountsCore |> group_by(Category) |> pivot_wider(Gene,names_from = Category, values_from = Genomes)  |>
	left_join(GenesID) |> distinct() |> write.csv("CoreCounts.csv", row.names =F, quote = F)

#### Accessory Genome Analysis ####
accessDf <- roaryOutputwithAncient |> filter(!(Gene %in% coreGenes))

# Preparing the data
accessData <- accessDf |> as.data.frame()
rownames(accessData) <- accessData[,1]
accessData <- accessData[,-1]
genomes <- colnames(accessData)
accessData <- apply(accessData, MARGIN = 1,FUN = as.numeric)
rownames(accessData) <- genomes
  
geneDist <- dist(accessData, method = "binary")

fit <- cmdscale(geneDist,eig = T, k = 4, add =T)
coord <- fit$points |> as_tibble()
coord <- coord[,1:4]
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100

coord <- coord |>
  mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata) |>
  select(Genome, host, sub_species, Rank_1_Lineage, Rank_5_Lineage, V1,V2, host) |>
  mutate(sub_species = factor(sub_species, levels = sort(unique(sub_species)))) |>
  mutate(Rank_1_Lineage = factor(Rank_1_Lineage, levels = sort(unique(Rank_1_Lineage)))) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage, levels = sort(unique(Rank_5_Lineage)))) |>
  mutate(host = ifelse(host == "Homo sapiens", "Human", ifelse(host == "NA"| is.na(host), "Unknown", "Non-Human Primate"))) |>  
  mutate(host = replace(host, is.na(host), "Unknown"))

# Plotting the results
ellipseData <- lapply(split.data.frame(coord, coord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,6:7])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows()
p1 <- coord |> #mutate(Assembled = ifelse(grepl("ERR", Genome), T, F)) |>
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_1_Lineage, group = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	scale_colour_manual("Lineage", values = coloursRank1) +
	guides(custom = guide_custom(customLegendSS14,title = "SS14 Sub-lineages", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols Sub-lineages", order = 2),
	       custom = guide_custom(customLegendTPE,title = "TPE Sub-lineages", order = 3), custom = guide_custom(customLegendTEN,title = "TEN Sub-lineages", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none()) +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	geom_point(data= ellipseData, shape = 23, size = 3, inherit.aes = F, aes(x = Centre_x, y = Centre_y, fill = Rank_5_Lineage), show.legend = F) +
	scale_fill_manual(values= c("red", "green", "blue", "purple")) +
	scale_colour_manual(values= c("red", "green", "blue", "purple")) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)"))  +
	theme(legend.position = "bottom", legend.title.position = "top")

p1PCOAAccess <- p1

# Finding the Nichols genomes which are hanging out with the SS14 clade
coord |> filter(V1 > 0, Rank_5_Lineage == 2) |> select(Genome) |> left_join(metadata) |> write.table(sep = "\t", file = "OddNichols.tab", quote = F, row.names = F)

NicholsWithSS14 <- coord |> filter(V1 > 0) |> filter(Rank_5_Lineage == 2) |> pull(Genome)
subDF <- accessData |> as.data.frame()
subDF$Genome <- rownames(subDF)
NicholsWithSS14Access <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = sapply(Rank_5_Lineage,TpalLintoSub)) |> mutate(Category = replace(Category, Genome %in% NicholsWithSS14, "Nichols-SS14")) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes") |> filter(grepl("Nichols|SS14", Category)) |> 
	pivot_wider(Gene,names_from = Category, values_from = Genomes) |> left_join(GenesID) |> distinct() 

subDF <- coreData |> as.data.frame()
subDF$Genome <- rownames(subDF)
NicholsWithSS14Core <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = sapply(Rank_5_Lineage,TpalLintoSub)) |> mutate(Category = replace(Category, Genome %in% NicholsWithSS14, "Nichols-SS14")) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes") |> filter(grepl("Nichols|SS14", Category)) |> 
	pivot_wider(Gene,names_from = Category, values_from = Genomes) |> left_join(GenesID) |> distinct() 

# Testing what happens if we remove gene_310
geneDist <- dist(accessData[,-which(colnames(accessData) == "group_310")], method = "binary")

fit <- cmdscale(geneDist,eig = T, k = 4, add =T)
coordMod <- fit$points |> as_tibble()
coordMod <- coordMod[,1:4]
coordMod$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100

coordMod <- coordMod |> #mutate(shape = replace(shape, grepl("draft", Genome), "Canadian")) |>
  mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata) |>
  select(Genome, host, sub_species, Rank_1_Lineage, Rank_5_Lineage, V1,V2, host) |>
  mutate(sub_species = factor(sub_species, levels = sort(unique(sub_species)))) |>
  mutate(Rank_1_Lineage = factor(Rank_1_Lineage, levels = sort(unique(Rank_1_Lineage)))) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage, levels = sort(unique(Rank_5_Lineage)))) |>
  mutate(host = ifelse(host == "Homo sapiens", "Human", ifelse(host == "NA"| is.na(host), "Unknown", "Non-Human Primate"))) |>  #|>
  mutate(host = replace(host, is.na(host), "Unknown")) #|>
  #left_join(amrInfo)

# Plotting the results

#baryCentre <- coordMod |> group_by(Rank_5_Lineage) |> summarize(V1 = mean(V1), V2 = mean(V2))

ellipseData <- lapply(split.data.frame(coordMod, coord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,6:7])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows()
p1 <- coordMod |> #mutate(Assembled = ifelse(grepl("ERR", Genome), T, F)) |>
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_1_Lineage, group = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	scale_colour_manual("Lineage", values = coloursRank1) +
	guides(custom = guide_custom(customLegendSS14,title = "Lineage 1", order = 1),custom = guide_custom(customLegendNichols,title = "Lineage 2", order = 2),
	       custom = guide_custom(customLegendTPE,title = "Lineage 3", order = 3), custom = guide_custom(customLegendTEN,title = "Lineage 4", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none()) +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	geom_point(data= ellipseData, shape = 23, size = 3, inherit.aes = F, aes(x = Centre_x, y = Centre_y, fill = Rank_5_Lineage), show.legend = F) +
	scale_colour_manual(values= c("red", "green", "blue", "purple")) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)"))  +
	theme(legend.position = "bottom", legend.title.position = "top")
ggsave("Figures/AccessoryPCOAno17.pdf", p1, width = 9, height = 6)


ggarrange(p1PCOAAccess, p1, common.legend = T, legend = "bottom", labels = "auto", align = "hv")

ggsave("Figures/AllPCOA.pdf", width = 12, height = 9)

# Finding out which genes are unique to SS14 vs Nichols vs TEN vs TPE
subDF <- accessData |> as.data.frame()
subDF$Genome <- rownames(subDF)
geneCountsAccess <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = sapply(Rank_5_Lineage,TpalLintoSub)) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes")

geneCountsAccess |> group_by(Category) |> pivot_wider(Gene,names_from = Category, values_from = Genomes)  |>
	left_join(GenesID) |> distinct() |> write.csv("AccessoryCounts.csv", row.names =F, quote = F)
