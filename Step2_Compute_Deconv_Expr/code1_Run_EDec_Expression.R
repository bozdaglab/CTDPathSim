## STAGE 2 
## CODE 1
### FINAL CODE FOR RUNNING EDec ON TCGA OVARIAN CANCER WITH EXPRESSION DATA ####
# Created By: Banabithi Bose
# Date Created: 4/23/2020
### EDec for gene expression ###
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


#### EDec Stage 2

# Run EDec stage 2 for tumor samples

load("~/DRUG_PIPELINE1/DATASETS/TCGA/OV.Rda")## TCGA data
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/Ensemble_hgnc.Rda")
GeneName<-data.table(row.names(RnaSeq_fpkm_data))#56830
RnaSeq_fpkm_data<-data.table(RnaSeq_fpkm_data)
RnaSeq_fpkm_data$ensembl_gene_id<-GeneName$V1
RnaSeq_fpkm_hgnc_data <- merge(Ensemble_hgnc,RnaSeq_fpkm_data,by = "ensembl_gene_id")
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data[,-1]#56537
library(dplyr)
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data %>% distinct(hgnc_symbol, .keep_all = TRUE)#37109

row.names(RnaSeq_fpkm_hgnc_data)<-RnaSeq_fpkm_hgnc_data$hgnc_symbol
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data[,-1]

## Need to find common samples with 450k samples
x1<-data.table(str_sub(colnames(RnaSeq_fpkm_hgnc_data),1,16))
x2<-data.table(str_sub(row.names(OV_Deconv_proportions),1,16))
common_samples<-merge(x1,x2,by="V1")#7
colnames(RnaSeq_fpkm_hgnc_data)<-str_sub(colnames(RnaSeq_fpkm_hgnc_data),1,16)
row.names(OV_Deconv_proportions)<-str_sub(row.names(OV_Deconv_proportions),1,16)

stage2_result_tumors = EDec::run_edec_stage_2(
  gene_exp_bulk_samples = RnaSeq_fpkm_hgnc_data[,common_samples$V1],
  cell_type_props = OV_Deconv_proportions[common_samples$V1,])

OV_Deconv_expression <-stage2_result_tumors$means
#row.names(BRCA_Deconv_expression)<-row.names(RnaSeq_fpkm_hgnc_data)
save(OV_Deconv_expression,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_expression450.Rda")

#########################################  ************  #############################################


##### FOR 27K DATA

## Reference methylation profile for TCGA patients
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/complete_ref_meth27.Rda")
cn<-read.xlsx(file = "~/PATIENT_CELL_PIPELINE1/ref_names.xlsx",header=T,sheetIndex = 1)
colnames(complete_ref_meth27)<-cn$V1
length(unique(cn$class))

# Create a vector of colors representing the class of each reference
ref_class_colors <- as.factor(as.character(cn$class))
levels(ref_class_colors) <- RColorBrewer::brewer.pal(8,"Pastel1")
ref_class_colors <- as.character(ref_class_colors)

# Create a color gradient to be used in a heatmap of correlations
color_gradient <- colorRampPalette(c("white","steelblue"))

# Compute correlation matrix and draw a heatmap
cors_ref_meth <- cor(as.data.frame(complete_ref_meth27),use = 'pairwise.complete.obs')

# pdf("PATIENT_CELL_PIPELINE1/PLOTS/BRCA_corr_complete_ref_meth27.pdf")
# gplots::heatmap.2(cors_ref_meth,
#                   trace="none",
#                   col=color_gradient(10),
#                   breaks=seq(0,1,0.1),
#                   margins=c(10,10),
#                   ColSideColors = ref_class_colors,
#                   RowSideColors = ref_class_colors)
# dev.off()

################ EDec STAGE 0 ##################
# Selecting marker loci based on comparisons of each class of reference against
# all other samples
## It is important that we select markers which are also present in complete tcga data with no "NA"s

load("~/DRUG_PIPELINE1/DATASETS/TCGA/OV_meth27.Rda")#
dim(meth27_data)#27578   613
meth27_TCGA<-as.data.frame(meth27_data[complete.cases(meth27_data),])# 

complete_ref<-complete_ref_meth27[complete.cases(complete_ref_meth27),]
dim(complete_ref_meth27)#25103    38
dim(complete_ref)#24635    38
x1<-data.table(row.names(meth27_TCGA))
x2<-data.table(row.names(complete_ref))
common_probes<-merge(x1,x2,by="V1")#19229
complete_ref_for_tcga<-complete_ref[common_probes$V1,]#[1] 19229    38
dim(complete_ref_for_tcga)
markers_ovr <- 
  EDec::run_edec_stage_0(reference_meth = as.data.frame(complete_ref_for_tcga),
                         reference_classes = as.character(cn$class),
                         max_p_value = 1e-5,
                         num_markers = 500,
                         version = "one.vs.rest")
save(markers_ovr,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/DATAFRAMES/OV_markers_ovr27.Rda")
# Selecting marker loci based on comparisons of between each pair of 
# reference classes
# markers_ep <- 
#   EDec::run_edec_stage_0(reference_meth = as.data.frame(complete_ref_for_tcga),
#                          reference_classes = as.character(cn$class),
#                          max_p_value = 1e-5,
#                          num_markers = 500,
#                          version = "each.pair")
# 
# save(markers_ep,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/DATAFRAMES/BRCA_markers_ep27.Rda")## 498 marker were selected

cors_ref_meth_markers_ovr <- cor(as.data.frame(complete_ref_for_tcga)[markers_ovr,])
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/cors_ref_meth_markers_ovr27.pdf")
gplots::heatmap.2(cors_ref_meth_markers_ovr,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(10,10),
                  ColSideColors = ref_class_colors,
                  RowSideColors = ref_class_colors)

dev.off()


# cors_ref_meth_markers_ep <- cor(as.data.frame(complete_ref_for_tcga)[markers_ep,])
# pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/cors_ref_meth_markers_ep27.pdf")
# gplots::heatmap.2(cors_ref_meth_markers_ep,
#                   trace="none",
#                   col=color_gradient(10),
#                   breaks=seq(0,1,0.1),
#                   margins=c(10,10),
#                   ColSideColors = ref_class_colors,
#                   RowSideColors = ref_class_colors)
# 
# dev.off()


#EDec stage 1 - Deconvolution of methylation profiles

# Run EDec stage 1 with all cell types using the markers_ovr loci
#Choosing a good number of cell types
set.seed(1)
stabilityResult <- EDec::estimate_stability(meth_bulk_samples = meth27_TCGA, 
                                            informative_loci = markers_ovr,
                                            possible_num_ct = 3:8,
                                            subset_prop = 0.8,
                                            num_subsets = 10,
                                            reps_per_subset = 1,
                                            max_its = 800,
                                            rss_diff_stop = 1e-8)
stabilityResult$most_stable_num_ct

set.seed(1)
stage1_result_ct = 
  EDec::run_edec_stage_1(meth_bulk_samples = meth27_TCGA, 
                         informative_loci = markers_ovr, 
                         num_cell_types = 8)

dim(stage1_result_ct$methylation)#21675      8
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
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/OV_cor_matrix_stage1_27.pdf")
gplots::heatmap.2(cors_deconv_refs_ct,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(4,12),
                  RowSideColors = ref_class_colors,
                  cellnote = best_cor_labels,
                  notecol="black")

dev.off()

## With 4 cell types
set.seed(1)
stage1_result_ct = 
  EDec::run_edec_stage_1(meth_bulk_samples = meth27_TCGA, 
                         informative_loci = markers_ovr, 
                         num_cell_types = 4)

dim(stage1_result_ct$methylation)#21675      4
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
for (i in 1:4){
  best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
}

# Plot correlation matrix
pdf("~/PATIENT_CELL_PIPELINE3/OVARIAN/PLOTS/OV_cor_matrix_stage1_stable27.pdf")
gplots::heatmap.2(cors_deconv_refs_ct,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(4,12),
                  RowSideColors = ref_class_colors,
                  cellnote = best_cor_labels,
                  notecol="black")

dev.off()


colnames(stage1_result_ct$methylation) <- c("Monocytes","Adipocytes","Vascular_endothelial_cells","Cortical_neuron")

colnames(stage1_result_ct$proportions) <- colnames(stage1_result_ct$methylation)
OV_Deconv_methylation <-stage1_result_ct$methylation
OV_Deconv_proportions <-stage1_result_ct$proportions

save(OV_Deconv_methylation,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_methylation27.Rda")
save(OV_Deconv_proportions,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_proportions27.Rda")


# Run EDec stage 2 for tumor samples
load("~/DRUG_PIPELINE1/DATASETS/TCGA/OV.Rda")## TCGA data
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/Ensemble_hgnc.Rda")
GeneName<-data.table(row.names(RnaSeq_fpkm_data))#56537
RnaSeq_fpkm_data<-data.table(RnaSeq_fpkm_data)
RnaSeq_fpkm_data$ensembl_gene_id<-GeneName$V1
RnaSeq_fpkm_hgnc_data <- merge(Ensemble_hgnc,RnaSeq_fpkm_data,by = "ensembl_gene_id")
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data[,-1]#56537
library(dplyr)
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data %>% distinct(hgnc_symbol, .keep_all = TRUE)#37109

row.names(RnaSeq_fpkm_hgnc_data)<-RnaSeq_fpkm_hgnc_data$hgnc_symbol
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data[,-1]

## Need to find common samples with 450k samples
x1<-data.table(str_sub(colnames(RnaSeq_fpkm_hgnc_data),1,16))
x2<-data.table(str_sub(row.names(OV_Deconv_proportions),1,16))
common_samples<-merge(x1,x2,by="V1")#367
colnames(RnaSeq_fpkm_hgnc_data)<-str_sub(colnames(RnaSeq_fpkm_hgnc_data),1,16)
row.names(OV_Deconv_proportions)<-str_sub(row.names(OV_Deconv_proportions),1,16)

stage2_result_tumors = EDec::run_edec_stage_2(
  gene_exp_bulk_samples = RnaSeq_fpkm_hgnc_data[,common_samples$V1],
  cell_type_props = OV_Deconv_proportions[common_samples$V1,])

OV_Deconv_expression <-stage2_result_tumors$means

save(OV_Deconv_expression,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_expression27.Rda")

dim(OV_Deconv_expression)#37109     4

