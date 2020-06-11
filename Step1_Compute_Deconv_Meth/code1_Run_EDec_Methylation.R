## STAGE 1 
## CODE 1
### RUNNING EDec ON TCGA OVARIAN CANCER WITH METHYLATION DATA ####
# Created By: Banabithi Bose
# Date Created: 4/23/2020
### EDec for methylation ###
#devtools::install_github("BRL-BCM/EDec")
library(EDec)
library('statmod')
library("stringr")
library(data.table)
library(sqldf)
library(reshape)
library(reshape2)
library(stringr)
library(janitor)
library(stringr)
library(plyr)
library(readr)
library(gsubfn)
library(pathfindR)
library(parallel)
library(xlsx)
##### FOR 450K DATA
## Reference methylation profile for TCGA patients
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/complete_ref_meth450.Rda")
cn<-read.xlsx(file = "~/PATIENT_CELL_PIPELINE1/ref_names.xlsx",header=T,sheetIndex = 1)
colnames(complete_ref_meth450)<-cn$V1
length(unique(cn$class))

# Create a vector of colors representing the class of each reference
ref_class_colors <- as.factor(as.character(cn$class))
levels(ref_class_colors) <- RColorBrewer::brewer.pal(8,"Pastel1")
ref_class_colors <- as.character(ref_class_colors)

# Create a color gradient to be used in a heatmap of correlations
color_gradient <- colorRampPalette(c("white","steelblue"))

# Compute correlation matrix and draw a heatmap
cors_ref_meth <- cor(as.data.frame(complete_ref_meth450),use = 'pairwise.complete.obs')

################ EDec STAGE 0 ##################
# Selecting marker loci based on comparisons of each class of reference against
# all other samples
## It is important that we select markers which are also present in complete tcga data with no "NA"s

load("~/DRUG_PIPELINE1/DATASETS/TCGA/OV_meth450.Rda")#485577    10
meth450_TCGA<-as.data.frame(meth450_data[complete.cases(meth450_data),])# 365575    10

complete_ref<-complete_ref_meth450[complete.cases(complete_ref_meth450),]
dim(complete_ref_meth450)
dim(complete_ref)#451988     38
x1<-data.table(row.names(meth450_TCGA))
x2<-data.table(row.names(complete_ref))
common_probes<-merge(x1,x2,by="V1")#338134
complete_ref_for_tcga<-complete_ref[common_probes$V1,]#[1] 368041     38

markers_ovr <- 
  EDec::run_edec_stage_0(reference_meth = as.data.frame(complete_ref_for_tcga),
                         reference_classes = as.character(cn$class),
                         max_p_value = 1e-5,
                         num_markers = 500,
                         version = "one.vs.rest")
save(markers_ovr,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/DATAFRAMES/OV_markers_ovr.Rda")


cors_ref_meth_markers_ovr <- cor(as.data.frame(complete_ref_for_tcga)[markers_ovr,])
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/cors_ref_meth_markers_ovr450.pdf")
gplots::heatmap.2(cors_ref_meth_markers_ovr,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(10,10),
                  ColSideColors = ref_class_colors,
                  RowSideColors = ref_class_colors)

dev.off()



#EDec stage 1 - Deconvolution of methylation profiles

# Run EDec stage 1 with all cell types using the markers_ovr loci
#Choosing a good number of cell types
set.seed(1)
stabilityResult <- EDec::estimate_stability(meth_bulk_samples = meth450_TCGA, 
                                            informative_loci = markers_ovr,
                                            possible_num_ct = 3:8,
                                            subset_prop = 0.8,
                                            num_subsets = 10,
                                            reps_per_subset = 1,
                                            max_its = 800,
                                            rss_diff_stop = 1e-8)
stabilityResult$most_stable_num_ct #8

set.seed(1)
stage1_result_ct = 
  EDec::run_edec_stage_1(meth_bulk_samples = meth450_TCGA, 
                         informative_loci = markers_ovr, 
                         num_cell_types = stabilityResult$most_stable_num_ct)

dim(stage1_result_ct$methylation)#395575      8
stage1_result_ct$methylation[1:10,]

stage1_result_ct$proportions[1:10,]

# Compute correlation between estimated methylation profiles,
# and reference methylation profiles
cors_deconv_refs_ct = cor(complete_ref_for_tcga[markers_ovr,],stage1_result_ct$methylation[markers_ovr,])

# Check what references had the highest correlation with each 
# of the estimated methylation profiles 
best_cors = rbind(apply(cors_deconv_refs_ct,2,which.max),
                  apply(cors_deconv_refs_ct,2,max))

best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_ct),
                         ncol=ncol(cors_deconv_refs_ct))
for (i in 1:8){
  best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
}

# Plot correlation matrix
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/OV_cor_matrix_stage1.pdf")
gplots::heatmap.2(cors_deconv_refs_ct,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(4,12),
                  RowSideColors = ref_class_colors,
                  cellnote = best_cor_labels,
                  notecol="black")

dev.off()

## with 3 cell types
set.seed(1)
stage1_result_ct = 
  EDec::run_edec_stage_1(meth_bulk_samples = meth450_TCGA, 
                         informative_loci = markers_ovr, 
                         num_cell_types = 3)

cors_deconv_refs_ct = cor(complete_ref_for_tcga[markers_ovr,],stage1_result_ct$methylation[markers_ovr,])

best_cors = rbind(apply(cors_deconv_refs_ct,2,which.max),
                  apply(cors_deconv_refs_ct,2,max))

best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_ct),
                         ncol=ncol(cors_deconv_refs_ct))
for (i in 1:3){
  best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
}
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/OV_cor_matrix_stage1_stable.pdf")
gplots::heatmap.2(cors_deconv_refs_ct,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(4,12),
                  RowSideColors = ref_class_colors,
                  cellnote = best_cor_labels,
                  notecol="black")

dev.off()

dim(stage1_result_ct$methylation)#395575      3
colnames(stage1_result_ct$methylation) <- c("Cortical_neurons","Vascular_endothelial_cells",
                                             "Adipocytes")

colnames(stage1_result_ct$proportions) <- colnames(stage1_result_ct$methylation)
OV_Deconv_methylation <-stage1_result_ct$methylation
OV_Deconv_proportions <-stage1_result_ct$proportions

save(OV_Deconv_methylation,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_methylation450.Rda")#395575      3
save(OV_Deconv_proportions,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_proportions450.Rda")#10   3

