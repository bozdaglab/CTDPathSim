## STAGE 1
## CODE 2
# Created By: Banabithi Bose
# Date Created: 4/23/2020
### EDec for methylation ###
#devtools::install_github("BRL-BCM/EDec")

##THIS CODE WILL PLOT CORRELATION PLOT WITH CELL TYPES USED FOR DECONVOLUTION

##### THIS FUNCTION WILL BE USED AFTER RUNING RunDeconvMethyl FUNCTION FOR THE CORRELATION PLOT

PlotDeconvMethyl<-function(Example.Methylation,markers_ovr,complete_ref_meth,stage1_result_ct,cn){

    chosen_markers_ovr <- intersect(row.names(Example.Methylation),markers_ovr)
    # Compute correlation between estimated methylation profiles,
    # and reference methylation profiles
    cors_deconv_refs_ct = cor(complete_ref_meth[chosen_markers_ovr,],stage1_result_ct$methylation[chosen_markers_ovr,])

    # Check what references had the highest correlation with each
    # of the estimated methylation profiles
    best_cors = rbind(apply(cors_deconv_refs_ct,2,which.max),
                    apply(cors_deconv_refs_ct,2,max))

    best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_ct),
                            ncol=ncol(cors_deconv_refs_ct))

    ## for (i in 1:ncol(stage1_result_ct[["proportions"]])){
    for (i in seq_len(ncol(stage1_result_ct[["proportions"]]))){
        best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
    }

    # Create a vector of colors representing the class of each reference
    ref_class_colors <- as.factor(as.character(cn$class))
    levels(ref_class_colors) <- RColorBrewer::brewer.pal(8,"Pastel1")
    ref_class_colors <- as.character(ref_class_colors)

    # Create a color gradient to be used in a heatmap of correlations
    color_gradient <- colorRampPalette(c("white","steelblue"))
    # Plot correlation matrix
    gplots::heatmap.2(cors_deconv_refs_ct,
                    trace="none",
                    col=color_gradient(10),
                    breaks=seq(0,1,0.1),
                    margins=c(4,4),
                    RowSideColors = ref_class_colors,
                    cellnote = best_cor_labels,
                    notecol="black")


}

# ##RUNNING THE FUNCTION AND SAVING THE PLOT

# load("~/CTDPathSim/inst/extdata/Example.Methylation.450K.rda")
# load("~/CTDPathSim/inst/extdata/complete_ref_meth450.rda")
# load("~/CTDPathSim/inst/extdata/markers_ovr_450.rda")
# load("~/CTDPathSim/inst/extdata/stage1_result_ct.rda")
# # Create a vector of colors representing the class of each reference
# cn<-xlsx::read.xlsx(file = "~/CTDPathSim/inst/extdata/ref_names.xlsx",header=T,sheetIndex = 1)
# pdf("~/CTDPathSim/inst/extdata/_cor_matrix_stage1_450K.pdf")
#
# PlotDeconvMethyl(Example.Methylation=Example.Methylation.450K,markers_ovr=markers_ovr,
#                  complete_ref_meth=complete_ref_meth450,stage1_result_ct=stage1_result_ct,cn=cn)
# # dev.off()
# #
