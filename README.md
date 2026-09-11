The data and scripts used in the publication titled "Evolutionary expansions of IF2 couple translation initiation to bacterial stress adaptation"

**Data** contains the following:
- genomes used in the study
- raw and analyzed growth curve results
- in vitro experiment results
- predictions for intrinsic disorder
- lengths of the extensions
- nusA gene sequences and their C-terminal region length
- predictions for phase separation
- selection analysis
- temperature and oxygen utilization traits from BacDive database

**Scripts** contains the following:
R scripts: ordered with number
- _01_gtdb_Filtering.R:_ filter the genomes to be used from GTDB
- _02_blast_Filtering.R:_ filter BLAST results to remove duplicates and nonspecific results
- _03_tree_length_properties.R:_ visualiza the species tree and map IF2 extension length, temperature and oxygen utilization predictions
- _04_dNdS_plot.R:_ visualize selection analysis on a representative sequence
- _05_ecoli_IDR.R:_ intrinsic disorder properties shown on E. IF2 predicted whole structure
- _06_idr_Composition.R:_ the frequencies of order and disorder promoting residues in the extension and conserved regions
- _07_preferences.R:_ the distribution of extension length within the groups based on environmental preference (i.e., temperature and oxygen utilization)
- _08_meanDisorder.R:_ the disorder prediction per residue and visualization as violin plots
- _09_meanDP.R:_ the droplet promotion (DP) prediction per residue and visualization as violin plots
- _10_growth_comparisons.R:_ different versions of the same growth curve script specifically used for different comparisons (e.g. truncation of N-termina, N-terminal swap or C-terminal addition, at different temperatures and pHs)
- _11_if2_Types.R:_ tree visualization with IF2 types based on structure prediction and clustering
- _12_invitro.R:_ the visualization of luciferase activity in in vitro transcription and translation system with different IF2 variants
- _13_correlations.R:_ the visualization of correlation of extension lengths with different environmental traits
- _14_nusA_if2.R:_ the visualization of correlation of lenghts of IF2 N-terminal extensions and nusA C-terminal regions
- _FuzDrop_auto.py:_ the automated prediction by submitting each sequence to the server as the software is not available

**Sources** contains the final data to creat plots in excel file, the scripts create final plots and the raw versions of the figure panels

**WT_comparisons** contains growth comparisons between E.coli MG1655 and SL598R strains at different temperature and pH 

**Supplementary Data 1:** The list of IF2 types

**Supplementary Data 2:** The list of organisms used in the dataset with their IF2 N- and C-terminal extension lengths, predicted optimal growth temperatures, and oxygen preference, NusA protein accession list and NusA C-terminal extension lengths

**Supplementary Data 3:** The list of strains, plasmid constructs, and primer sequences used in this study


  
