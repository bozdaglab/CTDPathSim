## STAGE 1
## CODE 1
# Created By: Banabithi Bose
# Date Created: 4/23/2020
### EDec for methylation ###
#devtools::install_github("BRL-BCM/EDec")

##### THIS FUNCTION CAN BE RUN FOR 450K AND 27K METHYLATION DATA WITH PROBES

#EDec stage 1 - Deconvolution of methylation profiles

# Run EDec stage 1 with all cell types using the markers_ovr loci

RunDeconvMethyl<-function(Example.Methylation,markers_ovr){
    chosen_markers_ovr <- intersect(row.names(Example.Methylation),markers_ovr)
    #Choosing a good number of cell types
    #set.seed(1)
    stabilityResult <- EDec::estimate_stability(meth_bulk_samples = Example.Methylation,
                                                informative_loci = chosen_markers_ovr,
                                                possible_num_ct = 3:8,
                                                subset_prop = 0.8,
                                                num_subsets = 10,
                                                reps_per_subset = 1,
                                                max_its = 800,
                                                rss_diff_stop = 1e-8)

    #set.seed(1)
    stage1_result_ct =
    EDec::run_edec_stage_1(meth_bulk_samples = Example.Methylation,
                            informative_loci = chosen_markers_ovr,
                            num_cell_types = stabilityResult$most_stable_num_ct)

    return(stage1_result_ct)
}

# ##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
# load("~/CTDPathSim/inst/extdata/Example.Methylation.450K.rda")
# load("~/CTDPathSim/inst/extdata/markers_ovr_450.rda")
#
# stage1_result_ct<-RunDeconvMethyl(Example.Methylation=Example.Methylation.450K,markers_ovr=markers_ovr)
#
# ##RESULT HAS BEEN SAVED IN  ~/CTDPathSim_Bioconductor/ WITH NAME stage1_result_ct_450K.rda
#
# ###Getting Deconv_methylation450.Rda
# Deconv_methylation<-stage1_result_ct$methylation
# Deconv_proportions <-stage1_result_ct$proportions
# ## EXAMPLE Deconv_methylation and Deconv_proportions ARE IN CTDPathSim_Bioconductor EXAMPLE DATASETS
# ## UNDER NAMES Example.Deconv.methylation.450K.rda AND Example.Deconv.proportions.450K.rd
