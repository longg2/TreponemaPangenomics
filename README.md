# Genomic approaches improve Treponema pallidum lineage classifications and predict emerging antimicrobial resistance insights

This repository contains all R scripts and metadata required to perform the pangenomic, phylogenetic, and genome analyses performed in the study. Several packages will be required, however, no other programs will be needed.

To get this point, you'll need to use the scripts found in [LongBioinformatics](https://github.com/longg2/LongBioinformatics). Specifically, you'll need to use:
	* CoreSNPPhylogenetics.sh to create the phylogeny
	* BWAmemMapping.sh for the 16S mapping
	* PoppunkAssigning.sh to update the PopPUNK database
	* ShovillWrapper.sh to create the genome assemblies, and
	* PanGenomeCreation.sh to create the pan-genome

The PopPUNK database we created for this paper can be found at [TpallidumLineages](https://github.com/Phylodynamics-and-AI-for-Public-health/TpallidumLineages).

