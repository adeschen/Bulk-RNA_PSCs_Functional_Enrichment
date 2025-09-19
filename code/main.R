
################################################################################
## Load the required packages
################################################################################
library(gprofiler2)
library(enrichViewNet)
library(ggplot2)
library(stringr)

################################################################################
## Load the required dataset
################################################################################

## Load the data from the differential expression TGFbeta versus control
dataTGF <- read.csv("./data/DEG_TGF_vs_control_from_Mucciolo_TableS1.csv",
                    header=T, stringsAsFactors = F)

## Load the data from the differential expression CM versus control
dataCM <- read.csv("./data/DEG_CM_vs_control_from_Mucciolo_TableS2.csv",
                    header=T, stringsAsFactors = F)

################################################################################
## Fix the gprofiler2 database version for functional enrichment to ensure
## reproducible results
################################################################################

## Force gprofiler2 to use the same database version as in the paper
gprofiler2::set_base_url(url="https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18")

################################################################################
## Functional enrichment for the up-regulated genes in
## the differential expression analysis TGFb vs control
################################################################################

## Retained only the significant up-regulated genes (FDR < 0.05 and logFC > 0)
significant_TGF_UP <-  dataTGF[which(dataTGF$FDR < 0.05 & dataTGF$logFC > 0), ]

## Run a functional enrichment analysis on the up-regulated genes
## The Ensembl ID are used as input
## The organism is mus musculus
## The sources are GO, KEGG, Reactome and Wikipathways
## The GO electronic annotations are excluded (exclude_iea=TRUE)
## The evidence codes are included in the results (evcodes=TRUE)
gostres_TGF_UP <- gost(query=list(TGFbeta_vs_ctl_UP=c(significant_TGF_UP$ensembl_gene_id)),
                        organism="mmusculus",
                        correction_method="g_SCS",
                        sources=c("GO", "KEGG", "REAC", "WP"),
                        user_threshold=0.05,
                        custom_bg=NULL,
                        significant=TRUE,
                        evcodes=TRUE,
                        exclude_iea=TRUE)

saveRDS(object=gostres_TGF_UP, file="results/EnrichmentTGF_Up.RDS")
write.csv(gostres_TGF_UP$results, "results/EnrichmentTGF_Up_summary.csv",
          row.names=FALSE)

## Remove TGF beta related dataset
rm(significant_TGF_UP)

################################################################################
## Functional enrichment for the down-regulated genes in
## the differential expression analysis TGFb vs control
################################################################################

## Retained only the significant up-regulated genes (FDR < 0.05 and logFC < 0)
significant_TGF_DOWN <- dataTGF[which(dataTGF$FDR < 0.05 & dataTGF$logFC < 0), ]

## Run a functional enrichment analysis on the up-regulated genes
## The Ensembl ID are used as input
## The organism is mus musculus
## The sources are GO, KEGG, Reactome and Wikipathways
## The GO electronic annotations are excluded (exclude_iea=TRUE)
## The evidence codes are included in the results (evcodes=TRUE)
gostres_TGF_DOWN <- gost(query=list(TGFbeta_vs_ctl_DOWN=c(significant_TGF_DOWN$ensembl_gene_id)),
                        organism="mmusculus",
                        correction_method="g_SCS",
                        sources=c("GO", "KEGG", "REAC", "WP"),
                        user_threshold=0.05,
                        custom_bg=NULL,
                        significant=TRUE,
                        evcodes=TRUE,
                        exclude_iea=TRUE)

saveRDS(object=gostres_TGF_DOWN, file="results/EnrichmentTGF_Down.RDS")
write.csv(gostres_TGF_DOWN$results, "results/EnrichmentTGF_Down_summary.csv",
            row.names=FALSE)

## Remove TGF beta related dataset
rm(significant_TGF_DOWN)
rm(dataTGF)

################################################################################
## Functional enrichment for the up-regulated genes in
## the differential expression analysis CM vs control
################################################################################

## Retained only the significant up-regulated genes (FDR < 0.05 and logFC > 0)
significant_CM_UP <- dataCM[which(dataCM$FDR < 0.05 & dataCM$logFC > 0), ]

## Run a functional enrichment analysis on the up-regulated genes
## The Ensembl ID are used as input
## The organism is mus musculus
## The sources are GO, KEGG, Reactome and Wikipathways
## The GO electronic annotations are excluded (exclude_iea=TRUE)
## The evidence codes are included in the results (evcodes=TRUE)
gostres_CM_UP <- gost(query=list(CM_vs_ctl_UP=c(significant_CM_UP$ensembl_gene_id)),
                       organism="mmusculus",
                       correction_method="g_SCS",
                       sources=c("GO", "KEGG", "REAC", "WP"),
                       user_threshold=0.05,
                       custom_bg=NULL,
                       significant=TRUE,
                       evcodes=TRUE,
                       exclude_iea=TRUE)

saveRDS(object=gostres_CM_UP, file="results/EnrichmentCM_Up.RDS")
write.csv(gostres_CM_UP$results, "results/EnrichmentCM_Up_summary.csv",
            row.names=FALSE)

rm(significant_CM_UP)

################################################################################
## Functional enrichment for the down-regulated genes in
## the differential expression analysis CM vs control
################################################################################

## Retained only the significant up-regulated genes (FDR < 0.05 and logFC < 0)
significant_CM_DOWN <- dataCM[which(dataCM$FDR < 0.05 & dataCM$logFC < 0), ]

## Run a functional enrichment analysis on the up-regulated genes
## The Ensembl ID are used as input
## The organism is mus musculus
## The sources are GO, KEGG, Reactome and Wikipathways
## The GO electronic annotations are excluded (exclude_iea=TRUE)
## The evidence codes are included in the results (evcodes=TRUE)
gostres_CM_DOWN <- gost(query=list(CM_vs_ctl_DOWN=c(significant_CM_DOWN$ensembl_gene_id)),
                         organism="mmusculus",
                         correction_method="g_SCS",
                         sources=c("GO", "KEGG", "REAC", "WP"),
                         user_threshold = 0.05,
                         custom_bg = NULL,
                         significant=TRUE,
                         evcodes=TRUE,
                         exclude_iea=TRUE)


saveRDS(object=gostres_CM_DOWN, file="results/EnrichmentCM_Down.RDS")
write.csv(gostres_CM_DOWN$results, "results/EnrichmentCM_Down_summary.csv",
            row.names=FALSE)

rm(significant_CM_DOWN)
rm(dataCM)

################################################################################
## Generate enrichment map using results of interest
################################################################################

## Terms of interest from GO:BP that are retained for the graph
## GO:0061564  GO:BP                      axon development
## GO:0007409  GO:BP                          axonogenesis
## GO:0050770  GO:BP            regulation of axonogenesis
## GO:0050772  GO:BP   positive regulation of axonogenesis
## GO:0048675  GO:BP                        axon extension
## GO:0030516  GO:BP          regulation of axon extension
## GO:0045773  GO:BP positive regulation of axon extension
## GO:0008366  GO:BP                     axon ensheathment
## GO:0098930  GO:BP                      axonal transport
## GO:0007411  GO:BP                         axon guidance
selected_terms <- c("GO:0061564", "GO:0007409", "GO:0050770", "GO:0050772",
                        "GO:0048675", "GO:0030516", "GO:0045773",
                        "GO:0008366", "GO:0098930", "GO:0007411")

set.seed(12514)

gg <- enrichViewNet::createEnrichMapMultiComplex(gostObjectList=list(gostres_TGF_UP,
                        gostres_TGF_DOWN, gostres_CM_DOWN),
                queryInfo = data.frame(queryName=c("TGFbeta_vs_ctl_UP",
                            "TGFbeta_vs_ctl_DOWN", "CM_vs_ctl_DOWN"),
                source=c("TERM_ID", "TERM_ID", "TERM_ID"),
                removeRoot=c(TRUE, TRUE, TRUE),
                termIDs=c(paste0(selected_terms, collapse=","),
                            paste0(selected_terms, collapse=","),
                            paste0(selected_terms, collapse=",")),
                groupName=c("TGFbeta Up", "TGFbeta Down", "CM Down"),
                stringsAsFactors = F))

gg <- gg + scale_fill_manual(name="Protocol",
            breaks=c("TGFbeta Up", "TGFbeta Down", "CM Down"),
            values=c("#D55E00", "#0072B2", "#009E73"),
            labels=c(expression(paste("TGF", beta, " up-regulated")),
                        expression(paste("TGF", beta, " down-regulated")),
                        "CM down-regulated")) +
            theme(legend.title = element_text(face="bold"))

pdf(paste0("results/Figure_2C.pdf"))
gg
invisible(dev.off())
