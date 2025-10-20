################################################################################
## Load the required packages
################################################################################
library(stringr)
library(gprofiler2)
library(enrichViewNet)
library(ggplot2)
library(scatterpie)
library(ggtangle)
library(ggrepel)
library(ggnetwork)
library(igraph)

################################################################################
## Load the required dataset
################################################################################

print("Loading dataset")

## Load the data from the differential expression TGFbeta versus control
dataTGF <- read.csv("../data/DEG_TGF_vs_control_from_Mucciolo_TableS1.csv",
                    header=T, stringsAsFactors=FALSE)

## Load the data from the differential expression CM versus control
dataCM <- read.csv("../data/DEG_CM_vs_control_from_Mucciolo_TableS2.csv",
                   header=T, stringsAsFactors=FALSE)


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

print("Enrichmnet on the TGFbeta vs control up-regulated genes")

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

saveRDS(object=gostres_TGF_UP, file="/results/EnrichmentTGF_Up.RDS")
write.csv(gostres_TGF_UP$results, "/results/EnrichmentTGF_Up_summary.csv",
          row.names=FALSE)

## Remove TGF beta related dataset
rm(significant_TGF_UP)

################################################################################
## Functional enrichment for the down-regulated genes in
## the differential expression analysis TGFb vs control
################################################################################

print("Enrichmnet on the TGFbeta vs control down-regulated genes")

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

saveRDS(object=gostres_TGF_DOWN, file="/results/EnrichmentTGF_Down.RDS")
write.csv(gostres_TGF_DOWN$results, "/results/EnrichmentTGF_Down_summary.csv",
          row.names=FALSE)

## Remove TGF beta related dataset
rm(significant_TGF_DOWN)
rm(dataTGF)

################################################################################
## Functional enrichment for the up-regulated genes in
## the differential expression analysis CM vs control
################################################################################

print("Enrichmnet on the CM vs control up-regulated genes")

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

saveRDS(object=gostres_CM_UP, file="/results/EnrichmentCM_Up.RDS")
write.csv(gostres_CM_UP$results, "/results/EnrichmentCM_Up_summary.csv",
          row.names=FALSE)

rm(significant_CM_UP)

################################################################################
## Functional enrichment for the down-regulated genes in
## the differential expression analysis CM vs control
################################################################################

print("Enrichmnet on the CM vs control down-regulated genes")

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
                        user_threshold=0.05,
                        custom_bg=NULL,
                        significant=TRUE,
                        evcodes=TRUE,
                        exclude_iea=TRUE)


saveRDS(object=gostres_CM_DOWN, file="/results/EnrichmentCM_Down.RDS")
write.csv(gostres_CM_DOWN$results, "/results/EnrichmentCM_Down_summary.csv",
          row.names=FALSE)

rm(significant_CM_DOWN)
rm(dataCM)

################################################################################
## Generate enrichment map using results of interest
################################################################################

print("Generate enrichment map")

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

termsIDs <- paste0(selected_terms, collapse=",")

## Generate a enrichment plot in igraph format
## Only TGF beta up, TGF beta down and CM down are used
## Minimum similarity to have a link present is set to 0.2
gg <- enrichViewNet::createEnrichMapMultiComplexAsIgraph(
    gostObjectList=list(gostres_TGF_UP, gostres_TGF_DOWN,
                        gostres_CM_DOWN),
    queryInfo=data.frame(queryName=c("TGFbeta_vs_ctl_UP",
                                       "TGFbeta_vs_ctl_DOWN", "CM_vs_ctl_DOWN"),
                           source=c("TERM_ID", "TERM_ID", "TERM_ID"),
                           removeRoot=c(TRUE, TRUE, TRUE),
                           termIDs=c(termsIDs, termsIDs, termsIDs),
                           groupName=c("TGFbeta Up", "TGFbeta Down", "CM Down"),
                           stringsAsFactors=FALSE),
    similarityCutOff=0.2)

set.seed(612)

## The igraph is transformed into a ggplot with some overspecialization
## First, the line width is function of the similarity coefficient
p<- ggplot2::ggplot(gg, layout=igraph::layout_with_fr) +
    geom_edge(aes(linewidth=weight), color="gray")

## The position of the "axon ensheathment" and "axonal transport" is changed
p$data$y[p$data$y > 0.6]  <- -0.7

p$data$x[p$data$label == "axonal transport"]  <- 11.2
p$data$x[p$data$label == "axon ensheathment"]  <- 9

## Extract information about the pie chart distribution of the nodes
## from the ggplot object (1 column contains all 3 groups)
## Add the information into the ggplot object (1 column per group)
pieInfo <- as.data.frame(do.call(rbind, V(gg)$pie))
colnames(pieInfo) <- c("TGFbeta Up", "TGFbeta Down", "CM Down")
p$data$`TGFbeta Up` <- pieInfo$`TGFbeta Up`
p$data$`TGFbeta Down` <- pieInfo$`TGFbeta Down`
p$data$`CM Down` <- pieInfo$`CM Down`

## Use scatterpie library to add the scatter pies plot into the graph
## The size of the nodes is function of the number of genes in term
## coord_fixed() is used to force a 1:1 ratio
p <- p + geom_scatterpie(aes(x=x, y=y, r=size/500),
            cols=c("TGFbeta Up", "TGFbeta Down", "CM Down"),
            legend_name = "Cluster", color=NA) +
        geom_scatterpie_legend(p$data$size/500,
            breaks=c(0, 0.1125, 0.225, 0.3375, 0.45),
            x=max(p$data$x)+0.9, y=min(p$data$y),
            labeller=function(x) {round(x*500)}, label_position="right") +
    coord_fixed() +
    guides(size="none") +
    guides(linewidth=guide_legend(ncol=3)) +
    scale_linewidth_continuous(name="Similarity coefficient") +
    theme(legend.position=c(0.855, 0.74),
            legend.key.width=unit(4, "mm"),
            legend.key.height=unit(4, "mm"),
            legend.spacing=unit(1, "mm"))

## Update colors for the nodes
p <- p + scale_fill_manual(name="Protocol",
            breaks=c("TGFbeta Up", "TGFbeta Down", "CM Down"),
            values=c("#c84d4c", "#3f78c1", "#28827a"),
            labels=c(expression(paste("TGF", beta, " up-regulated")),
                        expression(paste("TGF", beta, " down-regulated")),
                        "CM down-regulated")) +
        theme(legend.text=element_text(size=13),
                legend.title=element_text(size=14, face="bold"))

## Change position of the labels
p$data$nudge_y <- rep(-0.07, nrow(p$data))
p$data$nudge_y[p$data$label == "axonogenesis"] <- 0.4
p$data$nudge_y[p$data$label == "axon extension"] <- -0.1
p$data$nudge_y[p$data$label == "axon extension"] <- 0.13
p$data$nudge_y[p$data$label == "axon ensheathment"] <- -0.12
p$data$nudge_y[p$data$label == "axon guidance"] <- -0.28
p$data$nudge_y[p$data$label == "axon development"] <- -0.52
p$data$nudge_y[p$data$label == "positive regulation of axonogenesis"] <- 0.225
p$data$nudge_y[p$data$label == "regulation of axonogenesis"] <- -0.05
p$data$nudge_y[p$data$label == "regulation of axon extension"] <- -0.09

p$data$nudge_x<- rep(0, nrow(p$data))
p$data$nudge_x[p$data$label == "axonal transport"] <- -0.21
p$data$nudge_x[p$data$label == "axon extension"] <- 0.02
p$data$nudge_x[p$data$label == "axon guidance"] <- 0.1
p$data$nudge_x[p$data$label == "axon development"] <- 0.12
p$data$nudge_x[p$data$label == "axon ensheathment"] <- -0.1
p$data$nudge_x[p$data$label == "positive regulation of axonogenesis"] <- -0.32
p$data$nudge_x[p$data$label == "regulation of axonogenesis"] <- -0.3
p$data$nudge_x[p$data$label == "regulation of axon extension"] <- -0.04
p$data$nudge_x[p$data$label == "positive regulation of axon extension"] <- -0.2


p$data$label[p$data$label == "regulation of axon extension"] <-
                "regulation of axon\nextension"
p$data$label[p$data$label == "positive regulation of axonogenesis"] <-
                "positive regulation\nof axonogenesis"
p$data$label[p$data$label == "regulation of axonogenesis"] <-
                "regulation of\naxonogenesis"
p$data$label[p$data$label == "axon guidance"] <- "axon\nguidance"
p$data$label[p$data$label == "positive regulation of axon extension"] <-
                "positive regulation\nof axon extension"

## Use ggrepel library to add the labels for the nodes
ggg <- p + geom_text_repel(aes(x=x, y=y, label=label),
                nudge_y=p$data$nudge_y, nudge_x=p$data$nudge_x,
                min.segment.length=6, seed=121, size=5.2, lineheight = 0.7)

## Save graph in pdf
pdf(file="results/Figure_2C.pdf", width=7, height=3.3)
ggg
invisible(dev.off())

## Save graph in svg
svg(filename="results/Figure_2C.svg", width=7, height=3.3)
ggg
invisible(dev.off())