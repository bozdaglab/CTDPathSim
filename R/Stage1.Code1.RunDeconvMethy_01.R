## STAGE 1
## CODE 1
# Created By: Banabithi Bose
# Date Created: 4/23/2020
RunDeconvMethylFun <-
  function(Methylation.Sample,
           Reference.Methylation.Probes) {
    if ((is(Methylation.Sample, "SummarizedExpirament"))) {
      Methylation.Sample <- assay(Methylation.Sample)
    }
    
    
    chosen_Reference.Methylation.Probes <-
      intersect(row.names(Methylation.Sample),
                Reference.Methylation.Probes)
    
    stabilityResult <-
      estimate_stability(
        meth_bulk_samples = Methylation.Sample,
        informative_loci = chosen_Reference.Methylation.Probes,
        possible_num_ct = 3:8,
        subset_prop = 0.8,
        num_subsets = 10,
        reps_per_subset = 1,
        max_its = 800,
        rss_diff_stop = 1e-8
      )
    
    
    stage1_result_ct =
      run_edec_stage_1(
        meth_bulk_samples = Methylation.Sample,
        informative_loci = chosen_Reference.Methylation.Probes,
        num_cell_types = stabilityResult$most_stable_num_ct
      )
    
    stage1_result_ct$proportions <-
      round(stage1_result_ct$proportions, digits = 2)
    return(stage1_result_ct)
  }