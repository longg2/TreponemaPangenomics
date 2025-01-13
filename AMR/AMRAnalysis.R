# Required R Packages
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(pbapply)
library(lubridate)
library(emmeans)
library(colorspace)
library(grid)

# Functions to make life easier
RGItxtParser <- function(file){
  filename <- gsub("\\..*","",basename(file))
  tmp <- read.delim(file, header = T)[,-c(1,2)]
  tmp$Genome <- filename
  return(tmp)
}
BlastParser <- function(file){
  if(file.info(file)$size <= 95) return(NA)
  filename <- gsub("\\..*","",basename(file))
  tmp <- read.delim(file, header = T)
  tmp$Genome <- filename
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
coloursClusters <- c("#F0454B", "#2B8737", "#0369AC", "#92278F", "#00B2E3")
#### Loading in the metadata ####
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

# Want to do some filtering so that the geography is less confusing
smallCountries <- metadata |> count(geo_loc_name, name = "Count") |> filter(Count < 10) |>
  pull(geo_loc_name)
metadata <- metadata |>
  mutate(geo_loc_name = replace(geo_loc_name,
                                geo_loc_name %in% smallCountries | is.na(geo_loc_name),
                                "Other"))

# Reading in Beale 2021 metadata
beaeleMetadata <- read.csv("../MetadataExtra/Beale2021.csv") |> 
  rename(Genome = SRR.ENA_Accession, geo_loc_name = Geo_Country,
         collection_date = Sample_Year) |>
  mutate(collection_date = as.numeric(collection_date))

# Reading in the Lineage data
poppunkLineage <- read.csv("AllMasked_lineages.csv") |> mutate(id = gsub("_query","",id)) |>
#poppunkLineage <- read.csv("TomBealeMapped_lineages.csv") |> mutate(id = gsub("_query","",id)) |>
  rename(Genome = id) |> as_tibble() |> distinct()
#### The actual Script ####
RGIFiles <- list.files(c("AMRData"), full.names = T) # Getting the list of files

# Parsing the RGI results
RGIResults <- pblapply(RGIFiles, cl = 10, RGItxtParser) |> bind_rows() |> as_tibble() |>
  rename(LengthSubject = Percentage.Length.of.Reference.Sequence)

### Temp ###
RGIResults <- RGIResults |> 
	filter(!grepl("GCA_037762225|ERR4853567|ERR7123579|ERR7123580|ERR7123581|ERR7123582|ERR7123585|ERR7123588|ERR7123590|ERR7123591", Genome)) # These are making sure we remove either an engineered strain or extremely odd assemblies
### End Temp ###

RGIStrict <- RGIResults |> filter(Cut_Off == "Strict")
RGILoose <- RGIResults |> filter(Cut_Off == "Loose") |>
  mutate(LengthHit = abs(Hit_Start - Hit_End)) 

### Pulling out the strict!
RGIStrict |> count(Drug.Class) # Seeing which drug classes are in each
RGIStrict |> count(ARO)

### Loose Hits!
# Let's see if we can filter the Loose. Currently: Hit >= 1000bp and >= 50% ID
LongerLooseHits <- RGILoose |> 
  select(-c(Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence)) |>
  filter(LengthSubject >= 90) |>
  filter(Best_Identities >= 50, Best_Hit_Bitscore >= 100)

# What are the drug classes of the best hits?
LongerLooseHits |> count(Drug.Class)
AMRInfoOnly <- LongerLooseHits |> select(Drug.Class, Resistance.Mechanism, AMR.Gene.Family,
                          Antibiotic, Best_Hit_ARO, ARO) |> distinct()

RGILoose |> filter(LengthSubject >= 90) |> ggplot(aes(y = Best_Hit_Bitscore, x = Best_Identities)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 100, colour = "red") +
  geom_vline(xintercept = 50, colour = "red") +
  theme_classic() 
  #ggtitle("Tetracyline & Doxycycline Only")

# Now to look at the Tetracylines specifically
tmp <- RGILoose |> filter(grepl("tetra|doxy", Antibiotic))

tmp |> ggplot(aes(x = Best_Hit_Bitscore, y = Best_Identities)) +
  geom_point(alpha = 0.1) +
  geom_vline(xintercept = 175, colour = "red") +
  theme_classic() +
  ggtitle("Tetracyline & Doxycycline Only")
ggsave(file = "Figures/TetraDoxyFiltering.png", width = 6, height = 4)

tmp |> ggplot(aes(y = Best_Hit_Bitscore, x = LengthHit)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 175, colour = "red") +
  theme_classic() +
  ggtitle("Tetracyline & Doxycycline Only")

bestHitTet <- tmp |> filter(LengthSubject >= 90, Best_Hit_Bitscore >= 200)

### Plotting!
# Let's make a nice plot showing resistance amounts across classes

DrugClasses <- RGIStrict |> bind_rows(LongerLooseHits) |> #bind_rows(bestHitTet) |>
  select(Drug.Class,Cut_Off, Genome) |> group_by(Drug.Class, Cut_Off,Genome) |>
  distinct() |> group_by(Drug.Class,Cut_Off) |> count(Drug.Class, name = "Genomes") |> 
  arrange(-Genomes) |> mutate(Drug.Class = gsub(" antibiotic", "", Drug.Class))# |> 
  #mutate(Drug.Class = factor(Drug.Class, levels = Drug.Class))

DrugClasses$Drug.Class <- factor(DrugClasses$Drug.Class, levels = DrugClasses$Drug.Class |> unique())

ClassPlot <- DrugClasses |> ggplot(aes(x = Drug.Class, y = Genomes, fill = Cut_Off)) +
  theme_classic() +
  geom_col() +
  scale_fill_manual(values = c("Strict" = colour[2], "Loose" = colour[1])) +
  geom_hline(yintercept= length(unique(RGIResults$Genome)), lty = 2) +
  #geom_hline(yintercept= length(RGIFiles), lty = 2) +
  guides(fill = guide_legend(title = "Cut Off")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.25,0.5), legend.background = element_rect(colour = "black")) +
  xlab("Antibiotic Class")

ggsave(file = "Figures/Beale/DrugClasses.pdf", ClassPlot, width = 6, height = 4)

DrugClasses |> group_by(Drug.Class) |> select(-Cut_Off) |>
	summarize(Genomes = sum(Genomes), Prop = Genomes/813, PropMissing = 1 - Prop) |> 
	pivot_longer(-Drug.Class, names_to = "Status") |> filter(Status != "Genomes") |>
	ggplot(aes(x = Drug.Class, y = value, fill = Status)) +
	theme_void() +
	geom_col() +
	scale_fill_manual(values = c("#00b2de", "grey90")) +
	#geom_hline(yintercept = 1) +
	coord_polar(theta = "y") +
	facet_wrap("Drug.Class", scales = "free_x", ncol = 2) +
	guides(fill = guide_none())

ggsave(file = "Figures/Beale/DrugClassesPie.pdf", width = 2, height = 4)


# Should probably also show the drugs themselves and the methods!
indivDrug <- RGIStrict |> bind_rows(LongerLooseHits) |># bind_rows(bestHitTet) |>
  select(Antibiotic,Cut_Off, Genome) |>
  mutate(Antibiotic = strsplit(Antibiotic, "; ")) |> unnest(Antibiotic) |>
  group_by(Antibiotic, Cut_Off,Genome) |>
  distinct() |> group_by(Antibiotic,Cut_Off) |> count(Antibiotic, name = "Genomes") |> 
  arrange(-Genomes)# |> 

indivDrug$Antibiotic <- factor(indivDrug$Antibiotic,
                                          levels = unique(indivDrug$Antibiotic))

ResistancePlot <- indivDrug |> ggplot(aes(x = Antibiotic, y = Genomes,
                                           fill = Cut_Off)) +
  theme_classic() +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Strict" = colour[2], "Loose" = colour[1])) +
  geom_hline(yintercept= length(unique(RGIResults$Genome)), lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Antibiotic")

DrugMethod <- RGIStrict |> bind_rows(LongerLooseHits) |># bind_rows(bestHitTet) |>
  select(Resistance.Mechanism,Cut_Off, Genome) |>
  mutate(Resistance.Mechanism = strsplit(Resistance.Mechanism, "; ")) |> unnest(Resistance.Mechanism) |>
  group_by(Resistance.Mechanism, Cut_Off,Genome) |>
  distinct() |> group_by(Resistance.Mechanism,Cut_Off) |> count(Resistance.Mechanism, name = "Genomes") |> 
  arrange(-Genomes)# |> 

DrugMethod$Resistance.Mechanism <- factor(DrugMethod$Resistance.Mechanism,
                                          levels = unique(DrugMethod$Resistance.Mechanism))

ResistancePlot <- DrugMethod |> ggplot(aes(x = Resistance.Mechanism, y = Genomes,
                                           fill = Cut_Off)) +
  theme_classic() +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Strict" = colour[2], "Loose" = colour[1])) +
  geom_hline(yintercept= length(unique(RGIResults$Genome)), lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Resistance Mechanism")

# Showing the Best identities and bitscores
IdentBestHitPlot <- RGIStrict |> bind_rows(LongerLooseHits) |> #bind_rows(bestHitTet) |>
  ggplot(aes(x = Best_Hit_Bitscore, y = Best_Identities, colour = Cut_Off, shape = Drug.Class)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  scale_colour_manual(values = c("Strict" = colour[2], "Loose" = colour[1])) +
  guides(shape = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))

ggarrange(IdentBestHitPlot, ggarrange(ClassPlot, ResistancePlot, labels = c("B", "C"), ncol = 1,
          common.legend = F, legend = "none", align = "hv"),
          labels = c("A",""), legend = "bottom", common.legend = T)

ggsave("Figures/Beale/AMRSummarized.png", width = 9, height = 6)

# Saving the results
RGIStrict |> bind_rows(LongerLooseHits) |> #bind_rows(bestHitTet) |>
  write.table(file = "BealeAMRHits.tab", sep = "\t", row.names = F, quote = F)

amrInfo <- RGIStrict |> bind_rows(LongerLooseHits) |>
  select(Genome, Drug.Class) |>
  pivot_wider(id_cols = Genome, names_from = Drug.Class, 
              values_from = Drug.Class,
              values_fn = function(x){as.logical(length(x))},
              values_fill = F)
colnames(amrInfo) <- gsub(" antibiotic","", colnames(amrInfo))

amrInfo |> left_join(poppunkLineage) |> group_by(Rank_5_Lineage) |> select(where(is.logical)) |> pivot_longer(-Rank_5_Lineage) |>
	group_by(Rank_5_Lineage, name) |> count(value, name = "Genomes") |> pivot_wider(id_cols = c(Rank_5_Lineage, name), names_from = value, values_from = Genomes) |>
	ungroup() |> arrange(name, Rank_5_Lineage)

#### Looking at the macrolide resistance results ####
# Want to pull out the distribution of 2100 vs 2101 mutations and see if it matters
RGIStrict |> filter(Drug.Class == "macrolide antibiotic") |> select(Genome, SNPs_in_Best_Hit_ARO) |> rename(SNP = SNPs_in_Best_Hit_ARO) |>
	distinct() |>
	count(SNP)

# Some quick conversions making the results easier to read
conversionList = c("SS14", "Nichols", "TPE", "TEN")
coloursClustersConversion <- coloursClusters
names(coloursClustersConversion) <- c("SS14","Nichols","TPE","TEN")

# Merging the amrinfo and metadata together
macrolideModelData <- amrInfo |> #select(Genome, macrolide) |> 
  left_join(metadata |> select(Genome, collection_date, geo_loc_name, host)) |>
  left_join(beaeleMetadata, by = c("Genome"),
            relationship = "one-to-one") |>
  unite(collection_date, contains("collection"),na.rm = T) |>
  unite(geo_loc_name, contains("geo"),na.rm = T) |>
  mutate(geo_loc_name = replace(geo_loc_name, nchar(geo_loc_name) == 0, "Other")) |>
  mutate(collection_date = as.numeric(collection_date))|>
  mutate(across(where(is.logical), as.integer)) |>
  left_join(poppunkLineage |> select(Genome, contains("Rank"))) |>
	mutate(Rank_5_Lineage = conversionList[Rank_5_Lineage], Rank_5_Lineage = factor(Rank_5_Lineage, levels = conversionList, ordered = T))

smallCountries <- macrolideModelData |> count(geo_loc_name, name = "Count") |> filter(Count < 10) |>
  pull(geo_loc_name)
macrolideModelData <- macrolideModelData |>
  mutate(geo_loc_name = replace(geo_loc_name,
                                geo_loc_name %in% smallCountries | is.na(geo_loc_name),
                                "Other"))

countries <- macrolideModelData$geo_loc_name |> unique() |> sort()
ind <- which(countries == "Other")
countries <- c(countries[c(1:(ind - 1), (ind + 1):length(countries), ind)])
rm(ind)

macrolideModelData$geo_loc_name <- factor(macrolideModelData$geo_loc_name,
                                          levels = countries)
# Getting country dispersal
macrolideModelData |> group_by(Rank_5_Lineage) |> count(geo_loc_name, name = "Genomes") |>
	ggplot(aes(x = Rank_5_Lineage,y = Genomes, fill = geo_loc_name)) +
		geom_col() +
  		scale_fill_manual(values = countryColours, name = "Country") +
		theme_void() +
		coord_polar(theta = "y") +
		facet_wrap("Rank_5_Lineage", scales = "free")
macrolideModelData |> group_by(Rank_5_Lineage) |> count(geo_loc_name, name = "Genomes") |>
	mutate(Prop = Genomes/sum(Genomes)) |> filter(Rank_5_Lineage == "Nichols")


# Looking at the proportion/country
geoData <- macrolideModelData |> group_by(geo_loc_name) |>
  count(macrolide, name = "Genomes") |> mutate(Genomes = Genomes/sum(Genomes))
geoData$macrolide <- sapply(geoData$macrolide, function(x){ifelse(x,"Resistant",
                                                                  "Susceptible")})
topLabels <- macrolideModelData |> count(geo_loc_name, name = "Count")

geoPlot <- geoData |> ggplot(aes(x = geo_loc_name, y = Genomes, fill = macrolide)) +
  geom_col(position = "stack") +
  geom_text(inherit.aes = F, data = topLabels, aes(x = geo_loc_name, y = 1.05, label = Count), size = 8/.pt) +
  theme_classic() +
  scale_fill_manual(values = c("Resistant" = colour[1],
                               "Susceptible" = colour[4]),
                    name = "Azithromycin Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position= c(0.25,0.55),
        legend.background = element_rect(colour = "black")) +
  ylab("p(Genomes)") +
  xlab("Country")

ggsave("Figures/Beale/GeographyPlot.pdf", geoPlot, width = 6, height = 4)
  
# Performing the regression!
globglm <- macrolideModelData |> filter(collection_date >= 2000) |>
  glm(formula = macrolide ~ collection_date , family = "binomial") #|> 

BIC(globglm)
  #summary()

indivdata <- macrolideModelData |> filter(collection_date >= 2000) |>
  filter(!(geo_loc_name %in% c("Tanzania", "Sri Lanka", "Madagascar"))) 

## Testing if Lineage plays an important role here
indivglm <- indivdata |> mutate(geo_loc_name = relevel(geo_loc_name, "Other")) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage)) |>
  glm(formula = macrolide ~ collection_date + geo_loc_name + Rank_5_Lineage, family = "binomial") 

#indivglm <- indivdata |> mutate(geo_loc_name = relevel(geo_loc_name, "Other")) |>
#  glm(formula = macrolide ~ collection_date + geo_loc_name , family = "binomial") 
BIC(indivglm)

emmResult <- emmeans(indivglm, "geo_loc_name")

logOdds <- emmResult |> as_tibble() |> mutate(LCL = emmean - 1.96*SE, UCL = emmean + 1.96 * SE) |>
  mutate(geo_loc_name = factor(as.character(geo_loc_name), levels = names(countryColours))) |>
  ggplot(aes(x = emmean, y = geo_loc_name)) +
  geom_tile(aes(width = 1.96*SE*2, height = 0.5, fill = geo_loc_name), alpha = 0.5) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic() +
  scale_fill_manual(values = countryColours, name = "Country") +
  theme(legend.position = "none") +
  xlab("Log-Odds of Azithromycin Resistance") +
  ylab("Country")

emmResult <- emmeans(indivglm, "Rank_5_Lineage")

logOddsRank5 <- emmResult |> as_tibble() |> mutate(LCL = emmean - 1.96*SE, UCL = emmean + 1.96 * SE) |>
  mutate(Rank_5_Lineage = factor(Rank_5_Lineage)) |>
  ggplot(aes(x = emmean, y = Rank_5_Lineage)) +
  geom_tile(aes(width = 1.96*SE*2, height = 0.5, fill = Rank_5_Lineage), alpha = 0.5) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic() +
  scale_fill_manual(values = coloursClustersConversion, name = "Lineage") +
  theme(legend.position = "none") +
  xlab("Log-Odds of Azithromycin Resistance") +
  ylab("Lineage")

ggsave("Figures/Beale/EMMeansMacrolide.png", logOdds, width = 6, height = 4)

uniqueIndivDates <- indivdata |> select(collection_date, geo_loc_name, Rank_5_Lineage) |> mutate(Rank_5_Lineage = as.factor(Rank_5_Lineage)) |> distinct()
uniqueIndivDates <- cbind(uniqueIndivDates, predict(indivglm, type = "response", uniqueIndivDates, se.fit = T) |> as.data.frame())
  
# Doing this to make sure that the model looks nice in the regression line
indivRegPlot <- uniqueIndivDates |> group_by(Rank_5_Lineage, geo_loc_name) |>
  reframe(x = as.data.frame(spline(collection_date, y = fit))$x,
          y = as.data.frame(spline(collection_date, y = fit))$y,
          SE = as.data.frame(spline(collection_date, y = se.fit))$y)
indivRegPlot <- indivRegPlot |> 
  mutate(y = replace(y, y < 0, 0), y = replace(y, y > 1, 1)) |>
  mutate(loSE = y-1.96*SE, loSE = replace(loSE, loSE < 0, 0)) |>
  mutate(hiSE = y+1.96*SE, hiSE = replace(hiSE, hiSE > 1, 1))

logPlot <- macrolideModelData |> filter(collection_date >= 2000) |>
  mutate(Rank_5_Lineage = as.factor(Rank_5_Lineage)) |>
  mutate(geo_loc_name = factor(as.character(geo_loc_name), levels = names(countryColours))) |>
  ggplot(aes(x = collection_date, y = macrolide, colour = geo_loc_name)) +
  geom_jitter(height = 0, width = 1, alpha = 0.5 ) +
  geom_line(data = indivRegPlot, inherit.aes = F, 
              aes(x = x, y = y, colour = geo_loc_name),
              lwd = 1) +
  geom_ribbon(data = indivRegPlot, aes(y = y, x = x, 
                                       ymin = loSE, ymax = hiSE, fill = geo_loc_name),
              alpha = 0.1, colour = NA) +
  coord_cartesian(ylim = c(0,1)) +
  facet_wrap("Rank_5_Lineage") +
  theme_classic() +
  ylab("p(Azithromycin Resistance)") +
  xlab("Sampling Date") +
  scale_color_manual(values = countryColours, name = "Country") +
  scale_fill_manual(values = countryColours, name = "Country") +
  guides(colour = guide_legend(nrow = 4, title.hjust = 0.5),
         fill = guide_none()) +
  theme(legend.position= c(0.7,0.25),
        legend.background = element_rect(colour = "black")) 

# Plotting
ggarrange(logPlot,ggarrange(logOdds, logOddsRank5, labels = c("b","c"), align = "hv"),
          ncol = 1, labels = c("a", ""))

ggsave("Figures/Beale/UpdatedAMRPlot.pdf", width = 12, height = 9)
##### Looking at overall Drug Prevalence ####
rgiPlot <- RGIResults |> filter(Best_Hit_Bitscore >= 30, Best_Identities >= 20) |>
  filter(Cut_Off == "Strict" | LengthSubject >= 90) |>
  mutate(Drug.Class = replace(Drug.Class, grepl(";", Drug.Class), "multidrug")) |>
  mutate(Drug.Class = gsub(" antibiotic", "", Drug.Class)) |>
  mutate(Model_type = gsub(" model", "", Model_type)) |>
  ggplot(aes(y = Best_Hit_Bitscore, x = Best_Identities, colour = Drug.Class, shape = Model_type)) +
  geom_point() +
  guides(shape = guide_legend(ncol = 1, title = "Model"),
         colour = guide_legend(nrow = 4, title = "Drug Class")) +
  theme_classic() +
  scale_colour_manual(values = colour) +
  scale_x_continuous(limits = c(20,100)) +
  scale_y_log10(limits = c(30, 1e4)) +
  annotation_logticks(sides = "l") +
  ylab("Bitscore") +
  xlab("Percent Identity") +
  theme(legend.position = "bottom")

ggsave("Figures/Beale/RGIIndepthAll.png", rgiPlot, width = 12, height = 9)

RGIResults |> filter(Best_Hit_Bitscore >= 30, Best_Identities >= 20) |>
  filter(Cut_Off == "Strict" | LengthSubject >= 90) |> filter(Drug.Class == "tetracycline antibiotic") |>
  ggplot(aes(y = Best_Hit_Bitscore, x = Best_Identities, colour = Best_Hit_ARO, alpha = LengthSubject)) +
  geom_point() +
  theme_classic() +
  scale_colour_manual(values = colour) +
  guides(alpha = guide_legend(nrow = 3, title = "Relative Query Length (%)"), colour = guide_legend(nrow = 3, title = "Gene")) +
  ylab("Bitscore") +
  xlab("Percent Identity")+
  theme(legend.position = "bottom")

ggsave("Figures/RGITetracycline.png", width = 9, height = 6)
