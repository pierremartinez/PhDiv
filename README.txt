############################################################################
####                                Data                                ####
############################################################################

There are two datasets, with dedicated scripts for analysis in this study:
- sc_inductions : A single-cell RT-qPCR data we produced for EMT, DNA repair
  		  and glycolysis induction and subsequent targeted gene
		  expression analysis. The Biomark_output_table.xls file is
		  the output of the 48.48 chip RT-qPCR on the BiomarkHD,
		  while the emt_dnarep_glyco_induction_sc.csv serves as to
		  provide metadata.
- metadataset   : A dataset constructed from merging 8 public single-cell
  		  RNA-seq datasets. The data from this external sets are not
  		  provided, however the section "Necessary external Data"
		  below details how to download them and where to place
		  them. The AUCell activity scores for all 50 MSigDB
		  hallmark activities in each of the 8 datasets are stored
		  in the "preprocessed" folder as R data files.



#############################################################################
####                               Scripts                               ####
#############################################################################

There is only one script for the SC induction data, in the file
sc_inductions/activity_induction_prediction.R

Scripts are organised for the metadataset analysis are in the
metadataset/scripts directory and have to be run in a particular order
- metadataset_analyses.R, to normalise and integrate all datasets, as well as
  create all leave-one-out cluster and principal component designs.
- *_pca_dist.R and *_clust_dist.R to create set-specific cell-cell divergence
  lists (as R objects) and related plots. Default thresholds have been put to
  2% variance and 8 clusters, but this is customisable by changing the values
  of the thresh and nbest parameters in the script. A commented for loop
  allows to run the code with all thresholds used in the manuscript.
- merge_stats_n_plots.R comprises the analyses across all sets once all
  divergence measures are calculated.
- All plots are created in the metadataset/plots directory



#############################################################################
####                       Necessary external data                       ####
#############################################################################

Download the following files and place them in this folder
Files from the Broad Institute Single Cell Portal require authentication

##**********************##
## Filbin et al dataset ##
##**********************##
PortalK27M_Metadata.vh20180223.txt
PortalK27M_tSNE.vh20180223.txt

Downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP147/single-cell-analysis-in-pediatric-midline-gliomas-with-histone-h3k27m-mutation


##******************##
## Li et al dataset ##
##******************##
GSE81861_CRC_tumor_all_cells_FPKM.csv
GSE81861_series_matrix.txt
GSE81861_CRC_NM_all_cells_FPKM.csv

Downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81861


##**********************##
## Neftel et al dataset ##
##**********************##
IDHwtGBM.Metadata.SS2.txt

Downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma


##***************************************##
## Tirosh et al (melanoma) et al dataset ##
##***************************************##
GSE72056_melanoma_single_cell_revised_v2.txt

Downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056


##******************************************##
## Tirosh et al (oligodendroglioma) dataset ##
##******************************************##
cell_type_assignment_portal.txt
cell_cycle_coordinates_portal2.txt

Downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP12/oligodendroglioma-intra-tumor-heterogeneity


##**************************##
## Venteicher et al dataset ##
##**************************##
IDH_A_cell_type_assignment_portal_v2.txt

Downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP50/single-cell-rna-seq-analysis-of-astrocytoma
