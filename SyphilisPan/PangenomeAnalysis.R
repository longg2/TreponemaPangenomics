library(FactoMineR)
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
dedupMatrixWithCounts <- function(df){
   
  dfDedup <- df |> unique.matrix() # Easy, just getting the matrix
  
  duplicatedCounts <- df |> as.data.frame() |> 
    group_by(across(everything())) |>
    summarise(copies = n(), .groups = "drop") |>
    mutate(GroupNumber = 1:length(copies)) #|>  # Getting the actual counts
  
  # Now I need the ability to associate the duplicates to their clusters
  duplicatedSections <- df |> as.data.frame() |> 
    mutate(Genome = rownames(df)) |>
    group_by(across(-Genome)) |>
    group_split(.keep = F) |> lapply(function(x){unname(unlist(x))}) 
  
  # This makes the proper clusters
  groupingVariable <- NULL
  
  for(i in 1:length(duplicatedSections)){
    tmp <- sapply(duplicatedSections[[i]], function(x){
      vect <- rep(i,length(x))
      #names(vect) <- x
      return(vect)
    })
    groupingVariable <- c(groupingVariable, tmp)
  }
    #filter(copies > 1) |> arrange(-copies)
  
  # Getting the counts associated
  dfDedupCounted <- dfDedup |> as.data.frame() |> # appending the counts!
    full_join(duplicatedCounts) |> data.frame()
  
  rownames(dfDedupCounted) <- rownames(dfDedup)
  
  return(list(DedupMatrix = dfDedup, DedupMatrixCounts = dfDedupCounted, Groups = groupingVariable))
  
}

clusteringData <- function(coord, clusters.n = NA){
  # Need to know where to properly look
  indicesWithCoord <- which(sapply(coord, class) == "numeric") |> unname()
  
  if(is.na(clusters.n)){
    optimalCluster <- Optimal_Clusters_Medoids(coord[,indicesWithCoord],
                                           max_clusters = 10, 
                                           criterion = "silhouette", distance_metric = "euclidean")
  }else{
    clusters <- Clara_Medoids(coord[,indicesWithCoord], 
                              clusters.n, distance_metric = "euclidean", samples = 10, sample_size = 0.25)
    #coord$clusters <- predict(clusters, newdata = coord[,indicesWithCoord])
    return(clusters)
  }
}

expandingDuplicates <- function(Groups, df){
  tmp <- lapply(1:max(Groups), function(x){
    vecGen <- names(Groups[Groups == x]) # Getting the groups
    dataInd <- which(df$Genome %in% vecGen) # Assumes its got a Genome variable
    
    if(length(dataInd) == 0){
      return(NULL)
    }else if(length(vecGen) == 1){
      return(df[dataInd,])
    }
    
    # Now to actually print out the variable. First, we're going to get the column
    tmp <- df[rep(dataInd, length(vecGen)), -which(colnames(df) == "Genome")]
    tmp$Genome <- vecGen
    return(tmp)
  }) |> bind_rows()
  
  return(tmp) 
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

  # Setting up the number of permutations
  if(length(iter) == 0){
    iter = 1000
    if(iter < 1000){
      cat(paste0("Doing a minimum of 1000 permutations\n"))
    }else{
      iter >= 1000
      cat(paste0("Doing ",iter," permutations\n"))
    }
  }
 
  # If there are less permutations possible than iterations, we need to change it up!
  if(ncol(roaryMatrix) <= 9 & iter < factorial(ncol(roaryMatrix))){
	cat(paste("Calculating all of the possibilities"))
	perms <- permutations(ncol(roaryMatrix), ncol(roaryMatrix)) 
	perms <- split(perms, seq(1:nrow(perms)))
  }else{
	cat(paste("Randomizing the order of the genomes",iter,"times\n"))
	perms <- pblapply(1:iter, cl = 12, function(x) sample(1:ncol(roaryMatrix), ncol(roaryMatrix)))
  }

  # Now to calculate Heaps Law
  cat(paste("Performing the Heaps' law calculation \n"))
  coefDF <- pblapply(perms, cl = 12, function(perm){
		   #Making the matrix
		   permedMat <- roaryMatrix[,perm]
		   newGenes <- NULL

		   #iterating through the matrix
    		   for(j in 1:ncol(permedMat)){
    		   	foundGenes <- permedMat[,j]
    		   	newGenes <- c(newGenes,sum(foundGenes)) # How many new genes do we have?
    		   	
    		   	if(j == 1){
    		   	  totalGenes = newGenes
    		   	}else{
    		   	  totalGenes = c(totalGenes, sum(newGenes))
    		   	} 
    		   	indtosetzero <- which(foundGenes == 1)
    		   	permedMat[indtosetzero,] <- 0
    		   }

		   # Calculating the coefficients
    		   coefficients <- coef(nls(totalGenes ~ k*N^y, start = list(k = 1, y = 0)))
		   
		   return(coefficients)
  })  |> bind_rows()

  return(coefDF)
  
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

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')
smallColours <-  c("#3BB273", "#4D9DE0","#7768AE","#E1BC29","#E15554", "#D1A7A0")
conversionList = c("SS14", "Nichols", "Nichols - Madagascar", "TPE", "TEN")
conversionListNoSplit = c("SS14", "Nichols", "Nichols", "TPE", "TEN")
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

# Including the Poppunk info
poppunkLineage <- read.csv("Data/RevisedAllMask_lineages.csv") |>
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

metadata <- metadata |> full_join(poppunkLineage) |> left_join(genomeStats) |>
	mutate(SubspeciesSam = conversionListNoSplit[Rank_5_Lineage]) |>
  	mutate(SubspeciesSam = factor(SubspeciesSam, levels = c("SS14","Nichols", "TPE", "TEN")))

# Dealing with the colours now
coloursClusters <- c("#F0454B", "#2B8737", "#8A600D", "#0369AC", "#92278F")
#coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3") # Without Mad
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols", "Nichols - Madagascar","TPE","TEN")

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

coloursRank1 <- c(customLegendSS14$Colours,customLegendNicholsMad$Colours,customLegendNichols$Colours,customLegendTPE$Colours,customLegendTEN$Colours)

customLegendShape <- legendGrob(c("Human", "Non-Human Primate", "Unknown", "Lineage"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"),
				pch = c(16,17,15,18), gp = gpar(col = "black", cex = 2/3), nrow = 3)

savedLegendCustom <- guides(custom = guide_custom(customLegendSS14$Legend,title = "SS14", order = 1),custom = guide_custom(customLegendNichols$Legend,title = "Nichols", order = 2),
       custom = guide_custom(customLegendNicholsMad$Legend,title = "Nichols - Madagascar", order = 3),
       custom = guide_custom(customLegendTPE$Legend,title = "TPE", order = 4), custom = guide_custom(customLegendTEN$Legend,title = "TEN", order = 5),
       custom = guide_custom(customLegendShape, title = "Type", order = 6), colour = guide_none(), shape = guide_none()) 

# Making it easy to get the colours
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols", "Nichols - Madagascar","TPE","TEN")
############# Core Gene Presence #####
roaryOutput <- read.delim("Data/gene_presence_absence.Rtab") |> as_tibble()
#roaryOutput <- as_tibble(read.delim("Data/NewAssembled/gene_presence_absence.Rtab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("^X|ConfAncient|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\..*", "", colnames(roaryOutput)[2:ncol(roaryOutput)])

# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.95 * ncol(roaryOutput)) # Have incomplete assemblies, want to play it safe
accessGenes <- roaryOutput$Gene[!coreGenes]
coreGenes <- roaryOutput$Gene[coreGenes]

write.table(accessGenes, file = "AccessoryGenes.list", row.names = F, col.names = F, quote = F)

pangenomeCounts <- roaryOutput |> select(-Gene) |>
  summarize_all(sum) |>
  pivot_longer(everything(), names_to = "Genome", values_to = "Count") |>
  mutate(Genome = gsub("-.*","",Genome))

metadata <- metadata |> left_join(pangenomeCounts) 

# Now determining if genome length plays an important role

geneCountModel <- metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage]) |>
	lm(formula =  Count ~ Rank_5_Lineage + sum_len)
summary(geneCountModel)

metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> ggplot(aes(x = sum_len, y = Count, colour = as.factor(Rank_1_Lineage), group = as.factor(Rank_5_Lineage))) +
	geom_point() +
	scale_colour_manual("Lineage", values = coloursRank1) +
	savedLegendCustom +
	new_scale_colour() +
	geom_smooth(method = "lm", aes(colour = as.factor(Rank_5_Lineage)), show.legend = F) +
	scale_colour_manual(values= coloursClusters)  +
	theme(legend.position = "bottom", legend.title.position = "top")

metadata |> filter(Count > 940, sum_len <1.145e6, sum_len > 1.135e6) |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage]) |>
	lm(formula =  Count ~ Rank_5_Lineage + sum_len) |> emmeans(specs = "Rank_5_Lineage") |> pairs()

metadata |> group_by(Rank_5_Lineage) |> filter(!is.na(Count))|> summarize(median(Count))
metadata |> filter(!is.na(Count))|> summarize(median(Count))

# Loading in the gene information
GenesID <- read.delim("Data/RefSeqHitsTop.tab", header = F)[,c(1,2)]
colnames(GenesID) <- c("Gene", "Accession")
GenesID <- read.delim("Data/TranslatedTopHits.tab") |> rename(Accession = ProteinID) |> 
	right_join(GenesID |> mutate(Accession = gsub("\\.\\d$","",Accession)), relationship = "many-to-many")

# Quickly looking to see if TPE has a larger accessory content
subSep <- pblapply(conversionList, function(x){

		# Finding which columns are from our subtype
		colInd <- which(colnames(roaryOutput) %in% (metadata |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage]) |> 
							    filter(Rank_5_Lineage == x) |> pull(Genome)))
		
		# Making sure we're only working with those columns
		subRoary <- roaryOutput[,c(1,colInd)]

		# Removing genes we don't have
		subRoary <- subRoary[subRoary[,-1] |> rowSums() != 0,]

		# Doing the analysis!
		coreGenes <- rowSums(subRoary[,-1]) >= floor(0.95 * ncol(subRoary)) # Have incomplete assemblies, want to play it safe
		accessoryGenes <- !coreGenes # Have incomplete assemblies, want to play it safe
		
		out <- data.frame("Lineage" = x, "Core" = sum(coreGenes), "Accessory" = sum(accessoryGenes))
		return(out)
	}) |> bind_rows()

##### Doing the Core Genome PCoA ####
coreDf <- roaryOutput |> filter(Gene %in% coreGenes)

# Preparing the data
coreData <- as.data.frame(coreDf)
rownames(coreData) <- coreData[,1]
coreData <- coreData[,-1]
genomes <- colnames(coreData)
coreData <- apply(coreData, MARGIN = 1,FUN = as.numeric)
rownames(coreData) <- genomes

# Saving the core gene counts
subDF <- coreData |> as.data.frame()
subDF$Genome <- rownames(subDF)
geneCountsCore <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = conversionList[Rank_5_Lineage]) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes")

groupGeneCountsCore <- geneCountsCore |> group_by(Category) |> pivot_wider(Gene,names_from = Category, values_from = Genomes)  |>
	left_join(GenesID) |> distinct() 

groupGeneCountsCore |> write.csv("CoreCounts.csv", row.names =F, quote = F)
geneFiltering <- groupGeneCountsCore |> filter(grepl("tpr", Gene) | grepl("major outer sheath.*", Name)) |> pull(Gene)

coreDataSub <- coreData
maxCount <- coreDataSub |> colSums() |> max()
ind <- coreDataSub |> colSums() == maxCount
coreDataSub <- coreDataSub[,!ind]

# If I want to use a MDS plot
geneDist <- dist(coreDataSub, method = "binary")

fit <- cmdscale(geneDist,eig = T, k = 4, add =T)
coord <- fit$points |> as_tibble()
coord <- coord[,1:4]

coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100

# Making the dataframe
coord <- coord |>
  mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata) |>
  select(Genome, host, sub_species, Rank_1_Lineage, Rank_5_Lineage, V1,V2, host) |>
  mutate(sub_species = factor(sub_species, levels = sort(unique(sub_species)))) |>
  mutate(Rank_1_Lineage = factor(Rank_1_Lineage, levels = sort(unique(Rank_1_Lineage)))) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage, levels = sort(unique(Rank_5_Lineage)))) |>
  mutate(host = ifelse(host == "Homo sapiens", "Human", ifelse(host == "NA"| is.na(host), "Unknown", "Non-Human Primate"))) |>  
  mutate(host = replace(host, is.na(host), "Unknown"))

# Making the ellipse
ellipseData <- lapply(split.data.frame(coord, coord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,6:7])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows() |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage])

# Plotting the results
p1 <- coord |> 
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_1_Lineage, group = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	scale_colour_manual("Lineage", values = coloursRank1) +
	savedLegendCustom +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	geom_point(data= ellipseData, shape = 23, size = 3, inherit.aes = F, aes(x = Centre_x, y = Centre_y, fill = Rank_5_Lineage), show.legend = F) +
	scale_fill_manual(values= coloursClustersConversion) +
	scale_colour_manual(values= coloursClustersConversion) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)"))  +
	theme(legend.position = "bottom", legend.title.position = "top")

ggsave(p1, file = "Figures/CorePCoA.pdf", width = 9, height = 6)

#### Accessory Genome Analysis ####
accessDf <- roaryOutput |> filter(!(Gene %in% coreGenes))

# Preparing the data
accessData <- accessDf|> as.data.frame()
rownames(accessData) <- accessData[,1]
accessData <- accessData[,-1]
genomes <- colnames(accessData)
accessData <- apply(accessData, MARGIN = 1,FUN = as.numeric)
rownames(accessData) <- genomes

# Finding out which genes are unique to SS14 vs Nichols vs TEN vs TPE
subDF <- accessData |> as.data.frame()
subDF$Genome <- rownames(subDF)
geneCountsAccess <- subDF |> mutate(Genome = gsub("-.*","", Genome)) |> left_join(metadata |> select(Genome, Rank_5_Lineage)) |>
	mutate(Category = conversionList[Rank_5_Lineage]) |>
  	select(-c(Genome, Rank_5_Lineage)) |> group_by(Category) |> summarize(across(everything(), sum)) |>
  	pivot_longer(-Category, names_to = "Gene", values_to = "Genomes")

geneCountsGrouped <- geneCountsAccess |> group_by(Category) |> pivot_wider(Gene,names_from = Category, values_from = Genomes)  |>
	left_join(GenesID) |> distinct() 

geneCountsGrouped |> write.csv("AccessoryCounts.csv", row.names =F, quote = F)

# Making the plot!
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
tmpCoord <- coord |> filter(!(Rank_5_Lineage == 2 & V1 > 0))
ellipseData <- lapply(split.data.frame(coord, coord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,6:7])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows() |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage])

p1 <- coord |> 
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_1_Lineage, group = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75, show.legend = F) +
	scale_colour_manual("Lineage", values = coloursRank1) +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
#	geom_ellipse(lty = 4, data = ellipseDataNoOdd, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	geom_point(data= ellipseData, shape = 23, size = 3, inherit.aes = F, aes(x = Centre_x, y = Centre_y, fill = Rank_5_Lineage), show.legend = F) +
	savedLegendCustom +
	scale_fill_manual(values= coloursClustersConversion) +
	scale_colour_manual(values= coloursClustersConversion) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)"))  +
	theme(legend.position = "bottom", legend.title.position = "top")

p1PCOAAccess <- p1

# Removing the escaped tpr genes
geneDist <- dist(accessData[,-which(colnames(accessData) %in% "group_313")], method = "binary")

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
ellipseData <- lapply(split.data.frame(coordMod, coord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,6:7])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows() |> mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage])
p1 <- coordMod |> #mutate(Assembled = ifelse(grepl("ERR", Genome), T, F)) |>
  ggplot(aes(x = V1, shape = host, y = V2, label = Genome, colour = Rank_1_Lineage, group = Rank_5_Lineage)) +#, shape = Assembled)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	scale_colour_manual("Lineage", values = coloursRank1) +
	savedLegendCustom +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	geom_point(data= ellipseData, shape = 23, size = 3, inherit.aes = F, aes(x = Centre_x, y = Centre_y, fill = Rank_5_Lineage), show.legend = F) +
	scale_fill_manual(values= coloursClustersConversion) +
	scale_colour_manual(values= coloursClustersConversion) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)"))  +
	theme(legend.position = "bottom", legend.title.position = "top")
#ggsave("Figures/AccessoryPCOAno17.pdf", p1, width = 9, height = 6)

ggarrange(p1PCOAAccess, p1, common.legend = T, legend = "bottom", labels = "auto", align = "hv")

ggsave("Figures/AllPCOA.pdf", width = 12, height = 9)

############# Determining if it's an open or closed pan-genome #####
# Let's try removing the tpr/repeat genes
filtHeap <- GenesID |> filter(grepl("tpr|group_313", Gene) | grepl("major outer sheath", Name))

# Now to figure out how to make the models!
# Likely want to do multiple permutations to show that it doesn't matter who's first

heapFilteredResults <- panHeap(roaryOutput[!(roaryOutput$Gene %in% c(filtHeap,"group_313")),], 1e5)

heapResults <- panHeap(roaryOutput, 1e5)

heapSummarized <- heapResults |> 
  pivot_longer(everything(), values_to = "Value", names_to = "variable") |>
  group_by(variable) |>
  summarize(meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
                         meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
			 loValue = meanValue - 1.96 * seValue, hiValue = meanValue + 1.96 * seValue) |> as.data.frame()

heapSummarizedPlot <- heapSummarized |> mutate(variable = factor(variable, labels = c("K", TeX('$\\gamma$'))))

heapResults |> pivot_longer(everything(), values_to = "Value", names_to = "variable") |>
  mutate(variable = factor(variable, labels = c("K", TeX('$\\gamma$')))) |>
  ggplot(aes(x = Value, fill = variable)) +
  geom_density() +
  geom_rect(data = heapSummarizedPlot, inherit.aes = F,aes(ymin = -Inf, ymax = Inf,
                                       xmin = meanValue - 1.96*sdValue, xmax = meanValue + 1.96*sdValue), alpha = 0.5) +
  geom_vline(data = heapSummarizedPlot, aes(xintercept = meanValue), lty = 2) +
  facet_wrap("variable", scales = "free", labeller = label_parsed)+
  scale_fill_manual(values = smallColours, "Parameter") +
  theme(legend.position = "none")

ggsave("Figures/HeapParameters.pdf", width = 6, height = 4)

# Now to take a look at the individual sub-species (Updated: more efficient)
colnames(roaryOutput) <- roaryOutput |> colnames() |> gsub(pattern = "-.*", replacement = "")

subSep <- pblapply(c("SS14", "Nichols", "TPE", "TEN"), function(x){

		# Finding which columns are from our subtype
		colInd <- which(colnames(roaryOutput) %in% (metadata |> filter(SubspeciesSam == x) |> pull(Genome)))
		
		# Making sure we're only working with those columns
		subRoary <- roaryOutput[,c(1,colInd)]

		# Removing genes we don't have
		subRoary <- subRoary[subRoary[,-1] |> rowSums() != 0,]

		# Doing the analysis!
		heapResults <- panHeap(subRoary, 1e5)

		heapSummarized <- heapResults |> 
		  pivot_longer(everything(), values_to = "Value", names_to = "variable") |>
		  group_by(variable) |>
  		  summarize(meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
        	                 meanValue = mean(Value), sdValue = sd(Value), seValue = sdValue/sqrt(length(Value)),
				 loValue = meanValue - 1.96 * seValue, hiValue = meanValue + 1.96 * seValue, Sample = x, N = length(colInd)) |> as.data.frame()
		return(heapSummarized)
	}) |> bind_rows()

sepLin <- subSep |> bind_rows(heapSummarized |> mutate(Sample = "All", N = 544)) |> filter(Sample != "tpa") |> group_by(Sample, variable) |> 
	select(Sample, N, variable, meanValue, loValue, hiValue) |> group_by(Sample, variable) |>
	pivot_longer(-c(Sample, N, variable)) |> unite(New, variable, name) |> pivot_wider(id_cols = c(Sample, N), names_from = New, values_from = value) |>
	reframe(N = 1:N, y = k_meanValue*(N)^y_meanValue, ylo = k_loValue*(N)^y_loValue, yhi = k_hiValue*(N)^y_hiValue) |>
       	ggplot(aes(x = N, y = y, ymin = ylo, ymax = yhi, colour = Sample, fill = Sample)) +
	theme_classic() +
	geom_ribbon(alpha = 0.5, colour = NA) +
	scale_fill_manual(values = c(coloursClustersConversion, "All" = "black")) +
	scale_colour_manual(values = c(coloursClustersConversion, "All" = "black")) +
	geom_line() +
	#geom_hline(yintercept = sum(coreGenes), lty = 2, colour = "red") +
	scale_y_continuous(breaks = extended_breaks(n = 10)) +
	scale_x_continuous(breaks = extended_breaks(n = 10)) +
	xlab("Genomes") +
	ylab("Genes") +
	theme(legend.position = "bottom")

ggsave("Figures/HeapRegression.pdf",sepLin, width = 6, height= 4)

# Testing the significance of the tpr genes in the pan-genome characteristics

subSep <- pblapply(c("SS14", "Nichols", "TPE", "TEN"), function(x){

		# Finding which columns are from our subtype
		colInd <- which(colnames(roaryOutput) %in% (metadata |> filter(SubspeciesSam == x) |> pull(Genome)))
		
		# Making sure we're only working with those columns
		subRoary <- roaryOutput[,c(1,colInd)]

		# Removing genes we don't have
		subRoary <- subRoary[subRoary[,-1] |> rowSums() != 0,]

		# Doing the analysis!
		heapResults <- panHeap(subRoary, 1e5) |> mutate(Lineage = x, Filt = F)

		return(heapResults)
	}) |> bind_rows()

subSepFilt <- pblapply(c("SS14", "Nichols", "TPE", "TEN"), function(x){

		# Finding which columns are from our subtype
		colInd <- which(colnames(roaryOutput) %in% (metadata |> filter(SubspeciesSam == x) |> pull(Genome)))
		
		# Making sure we're only working with those columns
		subRoary <- roaryOutput[!(roaryOutput$Gene %in% c(filtHeap,"group_313")),c(1,colInd)]

		# Removing genes we don't have
		subRoary <- subRoary[subRoary[,-1] |> rowSums() != 0,]

		# Doing the analysis!
		heapResults <- panHeap(subRoary,1e5) |> mutate(Lineage = x, Filt = T)

		return(heapResults)
	}) |> bind_rows()

bind_rows(subSep, subSepFilt) |> mutate(Lineage = factor(Lineage, levels = c("SS14", "Nichols", "TPE", "TEN"))) |> lm(formula = y ~ Lineage * Filt) |> summary()
