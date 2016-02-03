#' PCA graphs on codon coverage
#' @param data a list of 2 data.frames:
#' one with the number of times each codon type is found in each ORF and
#' one with the number of reads for each codon in each ORF.
#' @param typeData a character string. It is used as title for the PCA.
#' Ex. typeData="codonCoverage"
#' @return a list of length 2:
#' PCA_scores - matrix of the scores on the first 4 principal components.
#' PCA_plots - a list of 5 PCA scatterplots.
#' @examples
#' #How to perform a PCA analysis based on codon coverage
#' #adapted from
#' #http://stackoverflow.com/questions/20260434/test-significance-of-clusters-on-a-pca-plot
#' #either get the codon frequency, coverage, and annotation using a function
#' #such as codonInfo in this package
#' #or create a list of matrices with the above information
#' data(codonDataCtrl)
#' codonData <- codonDataCtrl
#' codonUsage <- codonData[[1]]
#' codonCovMatrix <- codonData[[2]]
#'
#' #keep only genes with a minimum number of reads
#' nbrReadsGene <- apply(codonCovMatrix, 1, sum)
#' ixExpGenes <- which(nbrReadsGene >= 50)
#' codonCovMatrix <- codonCovMatrix[ixExpGenes, ]
#'
#' #get the PCA on the codon coverage
#' codonCovMatrixTransp <- t(codonCovMatrix)
#' rownames(codonCovMatrixTransp) <- colnames(codonCovMatrix)
#' colnames(codonCovMatrixTransp) <- rownames(codonCovMatrix)
#'
#' listPCACodonCoverage <- codonPCA(codonCovMatrixTransp,"codonCoverage")
#' print(listPCACodonCoverage[[2]])
#' #See aditional examples in the pdf manual
#' @export
#'

codonPCA <-
    function(
        data,
        typeData){
    Transcript <- NULL
    PC1 <- NULL
    PC2 <- NULL
    Cluster <- NULL
    PC3 <- NULL
    PC4 <- NULL
    .id <- NULL
    codonID <- NULL
    codon <- NULL

    pca <- stats::prcomp(data, retx=TRUE, center=TRUE)
    #screeplot(pca, type="l", main="Scree Plot")
    scores <- pca$x[, 1:4]

    percVar <- (pca$sdev)^2/sum(pca$sdev^2)
    km <- stats::kmeans(scores, centers=3, nstart=5)

    ggdata <- data.frame(scores, Cluster=km$cluster, Transcript=rownames(data))

    p12 <-
        ggplot2::ggplot(
            ggdata,
            ggplot2::aes_string(x="PC1", y="PC2")
        ) +
        ggplot2::geom_text(
            alpha=.4,
            size=3,
            ggplot2::aes(label=Transcript),
            vjust=2,
            hjust=2
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x=PC1,
                y=PC2,
                color=factor(Cluster)
            ),
            size=5,
            shape=20
        ) +
        ggplot2::ggtitle(typeData)+
        ggplot2::guides(
            color=ggplot2::guide_legend("Cluster"),
            fill=ggplot2::guide_legend("Cluster")
        ) +
        ggplot2::labs(
            x=paste("PC1(", round(percVar[1], 4), ")", collapse="", sep=""),
            y=paste("PC2(", round(percVar[2], 4), ")", collapse="", sep="")
        )
    p13 <-
        ggplot2::ggplot(
            ggdata,
            ggplot2::aes_string(x="PC1", y="PC3")
        ) +
        ggplot2::geom_text(
            alpha=.4,
            size=3,
            ggplot2::aes(label=Transcript),
            vjust=2,
            hjust=1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x=PC1,
                y=PC3,
                color=factor(Cluster)
            ),
            size=5,
            shape=20
        ) +
        ggplot2::ggtitle(typeData) +
        ggplot2::guides(
            color=ggplot2::guide_legend("Cluster"),
            fill=ggplot2::guide_legend("Cluster")
        ) +
        ggplot2::labs(
            x=paste("PC1(", round(percVar[1], 4), ")", collapse="", sep=""),
            y=paste("PC3(", round(percVar[3], 4), ")", collapse="", sep="")
        )
    p14 <-
        ggplot2::ggplot(
            ggdata,
            ggplot2::aes_string(x="PC1", y="PC4")
        ) +
        ggplot2::geom_text(
            alpha=.4,
            size=3,
            ggplot2::aes(label=Transcript),
            vjust=2,
            hjust=1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x=PC1,
                y=PC4,
                color=factor(Cluster)
            ),
            size=5,
            shape=20
        ) +
        ggplot2::ggtitle(typeData) +
        ggplot2::guides(
            color=ggplot2::guide_legend("Cluster"),
            fill=ggplot2::guide_legend("Cluster")
        ) +
        ggplot2::labs(
            x=paste("PC1(", round(percVar[1], 4), ")", collapse="", sep=""),
            y=paste("PC4(", round(percVar[4], 4), ")", collapse="", sep="")
        )
    p23 <-
        ggplot2::ggplot(
            ggdata,
            ggplot2::aes_string(x="PC2", y="PC3")
        ) +
        ggplot2::geom_text(
            alpha=.4,
            size=3,
            ggplot2::aes(label=Transcript),
            vjust=2,
            hjust=1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x=PC2,
                y=PC3,
                color=factor(Cluster)
            ),
            size=5,
            shape=20
        ) +
        ggplot2::ggtitle(typeData)+
        ggplot2::guides(
            color=ggplot2::guide_legend("Cluster"),
            fill=ggplot2::guide_legend("Cluster")
        ) +
        ggplot2::labs(
            x=paste("PC2(", round(percVar[2], 4), ")", collapse="", sep=""),
            y=paste("PC3(", round(percVar[3], 4), ")", collapse="", sep="")
        )
    p24 <-
        ggplot2::ggplot(
            ggdata,
            ggplot2::aes_string(x="PC2", y="PC4")
        ) +
        ggplot2::geom_text(
            alpha=.4,
            size=3,
            ggplot2::aes(label=Transcript),
            vjust=2,
            hjust=1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x=PC2,
                y=PC4,
                color=factor(Cluster)
            ),
            size=5,
            shape=20
        ) +
        ggplot2::ggtitle(typeData)+
        ggplot2::guides(
            color=ggplot2::guide_legend("Cluster"),
            fill=ggplot2::guide_legend("Cluster")
        ) +
        ggplot2::labs(
            x=paste("PC2(", round(percVar[2], 4), ")", collapse="", sep=""),
            y=paste("PC4(", round(percVar[4], 4), ")", collapse="", sep="")
        )
    #put the PCA plots in a list
    listPCAPlots <- list(p12, p13, p14, p23, p24)
    names(listPCAPlots) <- c("PC_1-2", "PC_1-3", "PC_1-4", "PC_2-3", "PC_2-4")

    returnList <- list(scores, listPCAPlots)
    names(returnList) <- c("PCA_scores","PCA_plots")

    return(returnList)
}
