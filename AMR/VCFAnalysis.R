library(FactoMineR)
library(factoextra)
library(colorspace)
library(lubridate)
library(pbapply)
library(ggrepel)
library(scales)
library(car)
library(ggrepel)
library(ggnewscale)
library(ggforce)
library(grid)
#library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(xtable)
vcfReading <- function(file){
  sampleName <- gsub("\\..*","",basename(file))
  vcfFile <- read.delim(file, header = T, comment.char = "#") |> mutate(CHROM = sampleName, ALT = as.character(ALT), REF = as.character(REF)) |>
 	mutate(ALT = replace(ALT, grepl("TRUE", ALT), "T"), ALT = replace(ALT, grepl("FALSE", ALT), "F"),
		REF = replace(REF, grepl("TRUE", REF), "T"), REF = replace(REF, grepl("FALSE", REF), "F")) 

  colnames(vcfFile) <- gsub(";.*","",colnames(vcfFile))

  altList <- strsplit(split = ",", vcfFile$ALT)

  # Getting the non-sample columns sorted out
  nonSample <- grepl("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT", colnames(vcfFile))

  vcfFile[,!nonSample] <- sapply(vcfFile[,!nonSample], function(x){
		 x <- as.numeric(gsub(":.*","",x))
		 z <- c()
		 for(i in 1:length(x)){
			#z[i] <- ifelse(x[i] == 0, "REF", altList[[i]][x[i]])
			z[i] <- ifelse(x[i] == 0, vcfFile$REF[i], altList[[i]][x[i]])
		 }
		 return(z)
					  })

  firstSample <- max(which(nonSample)) + 1
  vcfFile <- vcfFile[,c(1:2,firstSample:ncol(vcfFile))]
  
  return(vcfFile)
}
AAalignmenttoVCF <- function(file){
	# Read in the dataframe
	tmp <- readAAStringSet(file, format = "fasta") |> as.matrix() |> as.data.frame()
	rownames(tmp) <- gsub(":\\d$", "", rownames(tmp))
	gene <- gsub(".*/|AA.fasta", "", file)
	colnames(tmp) <- paste(sep = ":", gene, 1:ncol(tmp))

	# Filter out the monomorphic sites
	variants <- sapply(tmp, function(x) table(x) |> unlist() |> length() > 1)
	tmp <- t(tmp[,variants]) |> as.data.frame()
	tmp$CHROM <- gene
	tmp$POS <- gsub(".*:","", rownames(tmp)) |> as.numeric()
	
	return(tmp)
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
	#df <- data.frame("Centre" = ell.info$center, "xmax" = xmax, "xmin" = xmin, "ymax" = ymax, "ymin" = ymin)
	return(df)
}
TetWorkup <- function(vcf){
	colnames(vcf) <- gsub("^X|_R_", "", colnames(vcf))
	vcfDF <- vcf 
	
	vcfDF <- vcfDF |> unite(TMP,CHROM,POS, sep = ":")# |> pull(TMP)
	rownames(vcfDF) <- vcfDF$TMP
	vcfDF <- vcfDF |> select(-TMP)
	
	# Since this is a messy alignment, need to be more careful with how I'm removing certain genes
	noMissing <- sapply(vcfDF, function(x){
		      y <- table(x) |> unlist()# |> names()
		      propMissing <- y["*"]/sum(y)
		      ifelse(propMissing >= 0.25, F, T)
						  })
	
	vcfDF <- vcfDF[,noMissing]
	vcfDF$Position <- rownames(vcfDF)
	return(vcfDF)

}
#### Metadata Info ####
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
            '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
            '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
            '#000075', '#808080', '#f0f0ff', '#000000')

countryColours <- c("Canada" = darken("#e6194b", 0.25),
                    "Peru" = "#e6194b", 
                    "USA" = lighten("#e6194b", 0.25),
                    "China" = darken("#f58231",0.25),
                    "Japan" =  "#f58231", 
                    "Sri Lanka" = lighten("#f58231", 0.25),
                    "Madagascar" = darken("#3cb44b", 0.25), 
                    "Tanzania" = "#3cb44b", 
                    "Zimbabwe" = lighten("#3cb44b", 0.25),
                    "Hungary" = "#000075", 
                    "Ireland" = darken("#4363d8",0.25), 
                    "Italy" = "#4363d8",
                    "Portugal" = lighten("#4363d8", 0.25), 
                    "Russia" = darken("#008080", 0.25), 
                    "Spain" = "#008080",
                    "UK" = lighten("#008080", 0.25), 
                    "Other" = "#808080")
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
# Loading in the metadata
metadata <- read.delim("OtherData/T160D20240604OtherMetaData.tab", header = T)[-1,] 
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

# Beale Metadata
beaeleMetadata <- read.csv("../MetadataExtra/Beale2021.csv") |> 
  dplyr:::rename(Genome = SRR.ENA_Accession, geo_loc_name = Geo_Country,
         collection_date = Sample_Year) |>
  mutate(collection_date = as.numeric(collection_date))

# Reading in the Lineage data
poppunkLineage <- read.csv("AllMasked_lineages.csv") |> mutate(id = gsub("_query","",id)) |>
#poppunkLineage <- read.csv("TomBealeMapped_lineages.csv") |> mutate(id = gsub("_query","",id)) |>
	as_tibble() |> distinct()
colnames(poppunkLineage)[1] <- "Genome"

metadata <- metadata |> bind_rows(beaeleMetadata) #|> left_join(poppunkLineage) |> filter(!is.na(Rank_1_Lineage))

# Want to do some filtering so that the geography is less confusing
smallCountries <- metadata |> count(geo_loc_name, name = "Count") |> filter(Count < 10) |>
  pull(geo_loc_name)
metadata <- metadata |>
  mutate(geo_loc_name = replace(geo_loc_name,
                                geo_loc_name %in% smallCountries | is.na(geo_loc_name),
                                "Other"))

countries <- metadata$geo_loc_name |> unique() |> sort()
ind <- which(countries == "Other")
countries <- c(countries[c(1:(ind - 1), (ind + 1):length(countries), ind)])
rm(ind)

metadata$geo_loc_name <- factor(metadata$geo_loc_name,
                                          levels = countries)

# Dealing with the colours
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

customLegendShape <- legendGrob(c("Genome", "Lineage"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), pch = c(16,23), gp = gpar(col = "black", cex = 2/3), nrow = 2)

customLegendGov <- legendGrob(c("SS14","Nichols", "TPE", "TEN"), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                           pch = 20, gp = gpar(col = coloursClusters, pch = 1, cex = 2/3), nrow = 3)
##### Looking at the 16S rRNA Results ####
files <- list.files("16STetracycline", full.names = T, pattern = "vcf")
tmp <- pblapply(files, cl = 10, vcfReading)  # Reading in the data
names(tmp) <- gsub("16STetracycline/|.vcf.gz","",files)
#tmp <- tmp[sapply(tmp, nrow) > 0] # Filtering out the empty file

tmp <- pblapply(tmp,cl = 10, function(x){
	x <- x |> pivot_longer(cols = -c(CHROM, POS)) |>
	mutate(name = gsub("^X|\\..*|_R_","",name)) |>
	filter(!is.na(value)) |> distinct() #|>

	#pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") #|>
	return(x)
					  })

names(tmp) <- gsub(".vcf.gz","",basename(files))
vcfDF <- tmp |> bind_rows() |>
	 pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") |>
	 mutate(CHROM = gsub("NT", "",CHROM)) #|> filter(if_any(everything(), ~ !grepl("AB",.x)))

# Finding the strains to remove due to missing genes
filteringSamples <- vcfDF |> group_by(CHROM) |> select(-POS) |> mutate(across(everything(),~ .x == "AB")) |>
	group_by(CHROM) |> summarize(across(everything(), ~ sum(.x)), Length = length(CHROM)) |> pivot_longer(-c(CHROM, Length)) |>
	filter(value/Length == 1) |> pull(name) |> unique()
vcfDF <- vcfDF[,!(colnames(vcfDF) %in% filteringSamples)]

# finding those which are <5% indels
keeping <- vcfDF |> select(-CHROM, -POS) |> pivot_longer(everything()) |> group_by(name) |> summarize(value = sum(grepl("AB", value))/length(name)) |>
	filter(value == 0.00) |> pull(name)
# Getting things ready for the MCA
coord <- vcfDF |> unite(TMP, CHROM,POS, sep = ":") |> pull(TMP)
vcfDF <- vcfDF |> as.data.frame() |> select(-c(CHROM, POS))
rownames(vcfDF) <- coord
vcfDF <- vcfDF[,keeping]

vcfDF <- vcfDF[,colnames(vcfDF) %in% poppunkLineage$Genome]

# Did this in a weird way. Aligned each example from CARD separately against our 16S genes. Will have to look separately

# C. acnes. Only has G1032C.
grepSearch <- paste("Cacnes", 1032, sep = ":", collapse = "|")
grep(grepSearch, rownames(vcfDF)) # Comes up empty!

# Ecoli rrnB. Has A964G, G1053A, C1504T, A1055G. 
grepSearch <- paste("rrnB", c(964,1053,1504,1055), sep = ":", collapse = "|")
grep(grepSearch, rownames(vcfDF)) # Comes up empty!

# H. pylori. Has A965G, G966T, A967C. Also has 926 - 928.
grepSearch <- paste("Hpylori", c(926:928, 965:967), sep = ":", collapse = "|")
rowInd <- grep(grepSearch, rownames(vcfDF)) # Comes up empty!
vcfDF[rowInd[1],] |> t() |> as.data.frame() |> table()
vcfDF[rowInd[2],] |> t() |> as.data.frame() |> table()
vcfDF[rowInd[3],] |> t() |> as.data.frame() |> table()
vcfDF[rowInd,] |> t() |> as.data.frame() |> filter(`Hpylori:926` != "A" | `Hpylori:928` != "T" | `Hpylori:965` != "T") |> 
	write.csv("Hpylori.csv", quote = F)

# Looking at the 965 - 967 hotspot in Tpal from Tantalos et al. 2024 It's actually at 967 - 969 here
grepSearch <- paste("16STpal", c(967:969), sep = ":", collapse = "|")
rowInd <- grep(grepSearch, rownames(vcfDF)) # Comes up empty!
vcfDF[rowInd[1],] |> t() |> as.data.frame() |> table()
vcfDF[rowInd[2],] |> t() |> as.data.frame() |> table()
vcfDF[rowInd,] |> t() |> as.data.frame() |> table()
vcfDF[rowInd,] |> t() |> as.data.frame() |> filter(`16STpal:967` != "T" | `16STpal:969` != "A") |> write.csv("Hotspot.csv", quote = F)

##### Penicillin NT Results! ####
files <- list.files("Penicillin/NTMapped", full.names = T, pattern = "vcf")
tmp <- pblapply(files, cl = 10, vcfReading)  # Reading in the data
names(tmp) <- gsub("Penicillin/|.vcf","",files)
#tmp <- tmp[sapply(tmp, nrow) > 0] # Filtering out the empty file

tmp <- pblapply(tmp,cl = 10, function(x){
	x <- x |> pivot_longer(cols = -c(CHROM, POS)) |>
	mutate(name = gsub("^X|\\..*|_R_","",name)) |>
	filter(!is.na(value)) |> distinct() #|>

	#pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") #|>
	return(x)
					  })

names(tmp) <- gsub(".vcf.gz","",basename(files))
vcfDF <- tmp |> bind_rows() |>
	 pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") |>
	 mutate(CHROM = gsub("NT", "",CHROM)) #|> filter(if_any(everything(), ~ !grepl("AB",.x)))

# Finding the strains to remove due to missing genes
filteringSamples <- vcfDF |> group_by(CHROM) |> select(-POS) |> mutate(across(everything(),~ .x == "AB")) |>
	group_by(CHROM) |> summarize(across(everything(), ~ sum(.x)), Length = length(CHROM)) |> pivot_longer(-c(CHROM, Length)) |>
	filter(value/Length == 1) |> pull(name) |> unique()
vcfDF <- vcfDF[,!(colnames(vcfDF) %in% filteringSamples)]

# finding those which are <5% indels
keeping <- vcfDF |> select(-CHROM, -POS) |> pivot_longer(everything()) |> group_by(name) |> summarize(value = sum(grepl("AB", value))/length(name)) |>
	filter(value == 0.00) |> pull(name)
# Getting things ready for the MCA
coord <- vcfDF |> unite(TMP, CHROM,POS, sep = ":") |> pull(TMP)
vcfDF <- vcfDF |> as.data.frame() |> select(-c(CHROM, POS))
rownames(vcfDF) <- coord
vcfDF <- vcfDF[,keeping]
#vcfDF <- vcfDF |> select(-contains("draft")) # Removing draft genomes

test <- apply(vcfDF, MARGIN = 1, FUN = function(x){
		       counts <- table(x) |> unlist() |> sort(decreasing = T)
		       
		       minGen <- 2
		       if(length(counts) == 2 & any(counts < minGen)){
				return(FALSE)
		       }else if(length(counts) > 2 & sum(counts[2:length(counts)]) < minGen){
				return(FALSE)
		       }else if(length(counts) == 1){
				return(FALSE)
		       }else{
		       		return(T)
		       }
					  }) 

# Making sure we're not including any stragglers 
vcfDFMCA <- vcfDF[test,colnames(vcfDF) %in% poppunkLineage$Genome]

mcaResults <- MCA(t(vcfDFMCA), graph = F)

mcaCoord <- mcaResults$ind$coord |> as.data.frame()
mcaCoord$Genome <- rownames(mcaCoord)
mcaCoord <- mcaCoord |> mutate(Genome = gsub(":\\d","", Genome)) |> left_join(metadata) |> left_join(poppunkLineage) |> 
	mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage], Rank_5_Lineage = factor(Rank_5_Lineage, levels = conversionList, ordered = T)) |>
	distinct()

# Want to see about adding some arrows to the plot!
mcaCoordVar <- mcaResults$var$coord |> as.data.frame()
mcaCoordVar$SNP <- rownames(mcaCoordVar)

rowsImp <- rownames(mcaResults$var$contrib |> as.data.frame() |> arrange(-`Dim 2`) |> filter(`Dim 1` >= 1/nrow(mcaCoordVar) * 100 | `Dim 2` >= 1/nrow(mcaCoordVar) * 100))

mcaCoordVar <- mcaCoordVar[rowsImp,]

baryCentre <- mcaCoord |> group_by(Rank_5_Lineage) |> 
	summarize(`Dim 1` = mean(`Dim 1`), `Dim 2` = mean(`Dim 2`))

# Let's get the confidence ellipses calculated manually since stat_ellipse is dying on me....
ellipseData <- lapply(split.data.frame(mcaCoord, mcaCoord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,1:2])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows()

ntPBPLin1 <- mcaCoord |> ggplot(aes(x = `Dim 1`, y = `Dim 2`, colour = as.factor(Rank_1_Lineage), group = as.factor(Rank_5_Lineage))) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	guides(custom = guide_custom(customLegendSS14,title = "SS14 Sub-lineages", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols Sub-lineages", order = 2),
	       custom = guide_custom(customLegendTPE,title = "TPE Sub-lineages", order = 3), custom = guide_custom(customLegendTEN,title = "TEN Sub-lineages", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none(), fill = guide_none()) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	scale_colour_manual(values = coloursRank1) +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	scale_colour_manual(values= c("red", "green", "blue", "purple")) +
	geom_point(baryCentre, inherit.aes = F, mapping = aes(fill = as.factor(Rank_5_Lineage), x= `Dim 1`, y =`Dim 2`), shape = 23, colour = "black", size = 2) +
	scale_fill_manual(values = c("red", "green", "blue", "purple")) +
	theme_classic() +
	xlab(paste0("Dim 1 (",round(mcaResults$eig[1,2], digits = 2),"%)")) +
	ylab(paste0("Dim 2 (",round(mcaResults$eig[2,2], digits = 2),"%)")) +
	theme(legend.position = "bottom", legend.title.position = "top") 

ntResults <- fviz_mca_var(mcaResults, col.var = "contrib", ggtheme = theme_classic(), gradient.cols = c("blue", "green", "red")) + theme(plot.title = element_blank())

# Now to see the distribution of these SNPs across PopPUNK Lineages
vcfDFMCA$SNPLocation <- rownames(vcfDFMCA)
snpDist <- vcfDFMCA |> pivot_longer(-SNPLocation, names_to = "Genome", values_to = "SNP") |> left_join(poppunkLineage) |> 
	select(SNPLocation, Rank_5_Lineage, Genome, SNP) |>
	group_by(SNPLocation, Rank_5_Lineage) |> count(SNP, name = "Count") |>
	pivot_wider(c(SNPLocation, SNP), names_from = Rank_5_Lineage, values_from = Count, values_fill = 0)

#uniquePos <- mcaCoordVar$SNP |> gsub(pattern = "_.*",replacement = "") |> unique()

snpDist |> xtable() |> autoformat() |> print.xtable(booktabs = T, include.rownames = F, file = "MCASNPs.tex")

snpDist <- vcfDFMCA |> pivot_longer(-SNPLocation, names_to = "Genome", values_to = "SNP") |> left_join(poppunkLineage) |> 
	select(SNPLocation, Rank_1_Lineage, Genome, SNP) |>
	group_by(SNPLocation, Rank_1_Lineage) |> count(SNP, name = "Count") |>
	pivot_wider(c(SNPLocation, SNP), names_from = Rank_1_Lineage, values_from = Count, values_fill = 0)

#uniquePos <- mcaCoordVar$SNP |> gsub(pattern = "_.*",replacement = "") |> unique()
snpDist |> write.table("MCASNPsLin1.tab", row.names = F, quote = F, sep = "\t") #|> xtable() |> autoformat() |> print.xtable(booktabs = T, include.rownames = F, file = "MCASNPsLin1.tex")
vcfDFMCA

##### Dealing with the AAs now for Penicillin ####
files <- list.files("Penicillin/AA", full.names = T, pattern = "vcf")
tmp <- pblapply(files, cl = 10, vcfReading)  # Reading in the data
names(tmp) <- gsub("Penicillin/|.vcf","",files)

tmp <- pblapply(tmp,cl = 10, function(x){
	x <- x |> pivot_longer(cols = -c(CHROM, POS)) |>
	mutate(name = gsub("^X|\\..*|_R_","",name)) |>
	filter(!is.na(value)) |> distinct() #|>
	#pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") #|>
	return(x)
					  })

names(tmp) <- gsub("Penicillin/|.vcf","",files)
vcfDF <- tmp |> bind_rows() |>
	 pivot_wider(c(CHROM, POS), values_fn = function(x) paste(x, collapse = ","), values_fill = "AB") |>
	 mutate(CHROM = gsub("AA", "",CHROM)) #|> filter(if_any(everything(), ~ !grepl("AB",.x)))

coord <- vcfDF |> unite(TMP, CHROM,POS, sep = ":") |> pull(TMP)
vcfDF <- vcfDF |> as.data.frame() |> select(-c(CHROM, POS))
rownames(vcfDF) <- coord

vcfDF <- vcfDF |> select(-contains("draft")) # Removing draft genomes
# Quickly testing what will happen if I remove the genomes with missing genes
noMissing <- pbsapply(vcfDF, cl = 10, function(x){
	      y <- table(x) |> unlist() |> names()
		ifelse(sum(c("X","AB","*") %in% y), F, T)
					  })
vcfDF <- vcfDF[,noMissing]

test <- pbapply(vcfDF, cl = 10, MARGIN = 1, FUN = function(x){
		       counts <- table(x) |> unlist() |> sort(decreasing = T)
		       
		       minGen <- 2
		       if(length(counts) == 2 & any(counts < minGen)){
				return(FALSE)
		       }else if(length(counts) > 2 & sum(counts[2:length(counts)]) < minGen){
				return(FALSE)
		       }else if(length(counts) == 1){
				return(FALSE)
		       }else{
		       		return(T)
		       }
					  }) 
vcfDFMCA <- vcfDF[test,]

# Making the MCA
mcaResults <- MCA(t(vcfDFMCA), graph = F)

mcaCoord <- mcaResults$ind$coord |> as.data.frame()
mcaCoord$Genome <- rownames(mcaCoord)
mcaCoord <- mcaCoord |> mutate(Genome = gsub(":\\d","", Genome)) |> left_join(metadata) |> left_join(poppunkLineage) |>
	mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage], Rank_5_Lineage = factor(Rank_5_Lineage, levels = conversionList, ordered = T)) |>
	distinct()

# Want to see about adding some arrows to the plot!
mcaCoordVar <- mcaResults$var$coord |> as.data.frame()
mcaCoordVar$SNP <- rownames(mcaCoordVar)

rowsImp <- rownames(mcaResults$var$contrib |> as.data.frame() |> arrange(-`Dim 2`) |> filter(`Dim 1` >= 1/nrow(mcaCoordVar)  * 100| `Dim 2` >= 1/nrow(mcaCoordVar) * 100))
mcaCoordVar <- mcaCoordVar[rowsImp,]

baryCentre <- mcaCoord |> group_by(Rank_5_Lineage) |> summarize(`Dim 1` = mean(`Dim 1`), `Dim 2` = mean(`Dim 2`),)

ellipseData <- lapply(split.data.frame(mcaCoord, mcaCoord$Rank_5_Lineage), function(x){
			     tmp <- calculateEllipses(x[,1:2])
			     tmp$Rank_5_Lineage <- x$Rank_5_Lineage[1]
			     return(tmp)
	       }) |> bind_rows()
aaPBPLin1 <- mcaCoord |> ggplot(aes(x = `Dim 1`, y = `Dim 2`, colour = as.factor(Rank_1_Lineage), group = as.factor(Rank_5_Lineage))) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	#geom_point(data = baryCentre, shape = 23, size = 3, show.legend = F) +
	geom_jitter(width = 0.0015, height = 0.0015,alpha = 0.75) +
	guides(custom = guide_custom(customLegendSS14,title = "SS14 Sub-lineages", order = 1),custom = guide_custom(customLegendNichols,title = "Nichols Sub-lineages", order = 2),
	       custom = guide_custom(customLegendTPE,title = "TPE Sub-lineages", order = 3), custom = guide_custom(customLegendTEN,title = "TEN Sub-lineages", order = 4),
	       custom = guide_custom(customLegendShape, title = "Type", order = 5), colour = guide_none(), shape = guide_none(), fill = guide_none()) +
#	geom_point(data = mcaCoordVar, aes(x = `Dim 1`, y = `Dim 2`), inherit.aes = F, shape = "square") +
#	geom_text_repel(data = mcaCoordVar, aes(x = `Dim 1`, y = `Dim 2`, label = SNP), inherit.aes = F) +
	scale_x_continuous(breaks = pretty_breaks(n = 5)) +
	scale_y_continuous(breaks = pretty_breaks(n = 5)) +
	scale_colour_manual(values = coloursRank1) +
	new_scale_colour() +
	geom_ellipse(lty = 2, data = ellipseData, inherit.aes = F, aes(x0 = Centre_x, y0 = Centre_y, a = a, b = b, angle = angle, colour = Rank_5_Lineage), show.legend = F) +
	scale_colour_manual(values= c("red", "green", "blue", "purple")) +
	geom_point(baryCentre, inherit.aes = F, mapping = aes(fill = as.factor(Rank_5_Lineage), x= `Dim 1`, y =`Dim 2`), shape = 23, colour = "black", size = 2) +
	scale_fill_manual(values = c("red", "green", "blue", "purple")) +
	theme_classic() +
	xlab(paste0("Dim 1 (",round(mcaResults$eig[1,2], digits = 2),"%)")) +
	ylab(paste0("Dim 2 (",round(mcaResults$eig[2,2], digits = 2),"%)")) +
	theme(legend.position = "bottom", legend.title.position = "top") 

ggsave("Figures/Beale/PBPAA.pdf", aaPBPLin1, width = 9, height = 6)

# Plotting paper
ggarrange(ntPBPLin1, aaPBPLin1, nrow = 2, common.legend = T, legend = "bottom", labels = "auto", align = "hv")
ggsave("Figures/Beale/PBPAll.pdf", width = 9, height = 6)

# Plotting Variables
aaResults <- fviz_mca_var(mcaResults, col.var = "contrib", ggtheme = theme_classic(), gradient.cols = c("blue", "green", "red")) + theme(plot.title = element_blank())
ggarrange(ntResults, aaResults, ncol = 1, align = "hv", labels = "auto")
ggsave("Figures/Beale/MCAVarPositions.pdf", width = 12, height = 9)

vcfDFMCA$SNPLocation <- rownames(vcfDFMCA)
snpDist <- vcfDFMCA |> pivot_longer(-SNPLocation, names_to = "Genome", values_to = "SNP") |> left_join(poppunkLineage) |> 
	select(SNPLocation, Rank_5_Lineage, Genome, SNP) |>
	group_by(SNPLocation, Rank_5_Lineage) |> count(SNP, name = "Count") |>
	pivot_wider(c(SNPLocation, SNP), names_from = Rank_5_Lineage, values_from = Count, values_fill = 0)

#uniquePos <- mcaCoordVar$SNP |> gsub(pattern = "_.*",replacement = "") |> unique()

snpDist |> xtable() |> autoformat() |> print.xtable(booktabs = T, include.rownames = F, file = "MCAAAs.tex")

snpDist <- vcfDFMCA |> pivot_longer(-SNPLocation, names_to = "Genome", values_to = "SNP") |> left_join(poppunkLineage) |> 
	select(SNPLocation, Rank_1_Lineage, Genome, SNP) |>
	group_by(SNPLocation, Rank_1_Lineage) |> count(SNP, name = "Count") |>
	pivot_wider(c(SNPLocation, SNP), names_from = Rank_1_Lineage, values_from = Count, values_fill = 0)

#uniquePos <- mcaCoordVar$SNP |> gsub(pattern = "_.*",replacement = "") |> unique()
snpDist |> write.table("MCAAAsLin1.tab", row.names = F, quote = F, sep = "\t") #|> xtable() |> autoformat() |> print.xtable(booktabs = T, include.rownames = F, file = "MCASNPsLin1.tex")

