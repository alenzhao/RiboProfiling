#' Graphs of sample read counts (quality assesment)
#'
#' @param listCounts a list of data.frame objects.
#' It contains the counts on the genomic features.
#' Each data.frame in the list should have the same number of columns.
#' @param ixCounts a numeric (a vector of integers).
#' It contains the index of the columns containing counts in the dataFrame.
#' @param log2Bool a numeric, either 0 or 1.
#' 0 (default) for no log2 transformation and 1 for log2 transformation.
#' @return A list of pairs and boxplots between the counts data in each data.frame.
#' @examples
#' #read the BAM file into a GAlignments object using
#' #GenomicAlignments::readGAlignments
#' #the GAlignments object should be similar to ctrlGAlignments
#' data(ctrlGAlignments)
#' aln <- ctrlGAlignments
#'
#' #transform the GAlignments object into a GRanges object (faster processing)
#' alnGRanges <- readsToStartOrEnd(aln, what="start")
#'
#' #make a txdb object containing the annotations for the specified species.
#' #In this case hg19.
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' #Please make sure that seqnames of txdb correspond to
#' #the seqnames of the alignment files ("chr" particle)
#' #if not rename the txdb seqlevels
#' #renameSeqlevels(txdb, sub("chr", "",seqlevels(txdb)))
#' #get the flanking region around the promoter of the best expressed CDSs
#'
#' #get all CDSs by transcript
#' cds <- GenomicFeatures::cdsBy(txdb,by="tx",use.names=TRUE)
#'
#' #get all exons by transcript
#' exonGRanges <- GenomicFeatures::exonsBy(txdb,by="tx",use.names=TRUE)
#'
#' #get the per transcript relative position of start and end codons
#' cdsPosTransc <- orfRelativePos(cds, exonGRanges)
#'
#' #compute the counts on the different features after applying
#' #the specified shift value on the read start along the transcript
#' countsData <-
#'    countShiftReads(
#'          exonGRanges[names(cdsPosTransc)],
#'          cdsPosTransc,
#'          alnGRanges,
#'          -14
#'      )
#'
#' #now make the plots
#' listCountsPlots <- countsPlot(
#'    list(countsData[[1]]),
#'    grep("_counts$", colnames(countsData[[1]])),
#'    1
#' )
#' listCountsPlots
#' @export
#' @import ggplot2




countsPlot <-
    function(
        listCounts,
        ixCounts,
        log2Bool){
    stopifnot(inherits(listCounts,"list"))
    elemList <- length(listCounts)
    #check that the dimensions of all element in the list is the same
    if(elemList > 1){
        dimMat <- sapply(listCounts,ncol)
        if(length(which(dimMat != dimMat[1])) > 0){
            stop("All elements in listCounts don't have the same dimensions!\n")
        }
    }

    if(missing(log2Bool)){
        log2Bool <- 0
    }

    #intra list element plots
    if(length(ixCounts) > 1){
        #pairs on the counts in each data.frame
        myListPairs <- list(length(elemList))
        for(ixEl in 1:elemList){
            tmpCounts <- listCounts[[ixEl]][, ixCounts]
            if(log2Bool == 1){
                tmpCounts <- apply(tmpCounts, 2, as.numeric)
                tmpCounts <- log2(tmpCounts + 1)
                tmpCounts <- data.frame(tmpCounts)
                yText <- "log2_counts"
            }
            else{
                yTest <- "counts"
            }
            pPairs <- funcPlotPairs(tmpCounts, names(listCounts)[ixEl])
            myListPairs[[ixEl]] <- pPairs
        }

    }

    #overall boxplot of counts
    countGgplot <- funcDataFrameGGPlot2(listCounts, ixCounts, log2Bool)
    figBoxplot <- funcBoxplot(countGgplot, "", yText, "")
    return(list(myListPairs, figBoxplot))

}
