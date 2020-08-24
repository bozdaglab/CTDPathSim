##### This is for pathways with cancer genes from patient specific DE genes  #####
## THIS CODE IS FOR FINDING PATHWAYS FOR PATIENTS DE/DM MUTATED GENES
## STAGE 3
## CODE 3
# Created By: Banabithi Bose
# Date Created: 4/28/2020
##### This is for pathways with cancer genes from patient specific DE genes  #####
#
# library('statmod')
# library("stringr")
# library(data.table)
# library(sqldf)
# library(reshape)
# library(reshape2)
# library(stringr)
# library(janitor)
# library(stringr)
# library(plyr)
# library(readr)
# library(gsubfn)
# library(pathfindR)
# library(parallel)

GetPath<-function(CancerGeneList,pat_de_genes){
    common_gene_mut_de <-merge(pat_de_genes,CancerGeneList,by.x = "gene",by.y = "V1")
    reactome<-enrichment(
    input_genes=common_gene_mut_de$gene,
    genes_by_term = pathfindR.data::reactome_genes,
    term_descriptions = pathfindR.data::reactome_descriptions,
    adj_method = "bonferroni",
    enrichment_threshold = 0.05,
    sig_genes_vec=common_gene_mut_de$gene,
    background_genes=unlist(pathfindR.data::reactome_genes))

    if (is.null(nrow(reactome))== TRUE){
        print("no patient reactome")
        pat_reactome<-data.frame("","","")
        colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
        return(pat_reactome)
    }else{
        d1<-data.frame(reactome$Term_Description)
        d2<-data.frame(reactome$adj_p)
        d3<-data.frame(reactome$ID)
        pat_reactome<-cbind(d3,d1,d2)
        colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
        return(pat_reactome)
    }
}

# ##RUNNING THE FUNCTION
#load("~/CTDPathSim/inst/extdata/Example.CancerGeneList.rda")
#load("~/CTDPathSim/inst/extdata/Example.patient1.de.genes.rda")
#x<-GetPath(CancerGeneList,pat_de_genes)
# ### WE NEED TO PROVIDE AN EXAMPLE FOR CELL LINE1 load("~/CTDPathSim_Bioconductor/Example.cellLine1.de.genes.rda")



