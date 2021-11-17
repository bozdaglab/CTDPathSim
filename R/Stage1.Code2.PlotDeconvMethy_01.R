## STAGE 1
## CODE 2
# Created By: Banabithi Bose
# Date Created: 4/23/2020
PlotDeconvMethylFun <-
  function(Methylation.Sample,
           Reference.Methylation.Probes,
           Reference.Methylation.CellTypes,
           stage1_result_ct,
           Reference.CellTypes.Names) {
    if ((is(Methylation.Sample, "SummarizedExpirament"))) {
      Methylation.Sample <- assay(Methylation.Sample)
    }
    if ((is(Reference.Methylation.CellTypes, "SummarizedExpirament"))) {
      Reference.Methylation.CellTypes <-
        assay(Reference.Methylation.CellTypes)
    }
    
    if ((is(Reference.CellTypes.Names, "SummarizedExpirament"))) {
      Reference.CellTypes.Names <- assay(Reference.CellTypes.Names)
    }
    
    
    chosen_Reference.Methylation.Probes <-
      intersect(row.names(Methylation.Sample),
                Reference.Methylation.Probes)
    # Compute correlation between estimated methylation profiles,
    # and reference methylation profiles
    cors_deconv_refs_ct = cor(Reference.Methylation.CellTypes[chosen_Reference.Methylation.Probes, ],
                              stage1_result_ct$methylation[chosen_Reference.Methylation.Probes, ])
    
    # Check what references had the highest correlation with each
    # of the estimated methylation profiles
    best_cors = rbind(apply(cors_deconv_refs_ct, 2, which.max),
                      apply(cors_deconv_refs_ct, 2, max))
    
    best_cor_labels = matrix("",
                             nrow = nrow(cors_deconv_refs_ct),
                             ncol = ncol(cors_deconv_refs_ct))
    
    for (i in seq_len(ncol(stage1_result_ct[["proportions"]]))) {
      best_cor_labels[best_cors[1, i], i] = as.character(round(best_cors[2, i], 2))
    }
    
    # Create a vector of colors representing the class of each reference
    ref_class_colors <-
      as.factor(as.character(Reference.CellTypes.Names$class))
    levels(ref_class_colors) <- RColorBrewer::brewer.pal(8, "Pastel1")
    ref_class_colors <- as.character(ref_class_colors)
    
    # Create a color gradient to be used in a heatmap of correlations
    color_gradient <- colorRampPalette(c("white", "steelblue"))
    # Plot correlation matrix
    gplots::heatmap.2(
      cors_deconv_refs_ct,
      trace = "none",
      col = color_gradient(10),
      breaks = seq(0, 1, 0.1),
      margins = c(15, 15),
      RowSideColors = ref_class_colors,
      cellnote = best_cor_labels,
      notecol = "black"
    )
    
  }
