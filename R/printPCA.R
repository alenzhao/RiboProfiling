#' Plots the PCA scatterplots produced by codonPCA function.
#' @param listPCAGraphs a list of 5 PCA ggplot scatterplots.
#' @return a unique plot with the 5 PCA scatterplots.
#' @examples
#' #How to perform a PCA analysis based on codon coverage
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
#' printPCA(listPCACodonCoverage[[2]])
#' @export

printPCA <- function(
    listPCAGraphs){
    stopifnot(inherits(listPCAGraphs, "list"))
    stopifnot(inherits(listPCAGraphs[[1]], "ggplot"))
    suppressWarnings(ggbio::arrangeGrobByParsingLegend(
        listPCAGraphs,
        ncol=3,
        nrow=2,
        legend.idx=1
        ))
}