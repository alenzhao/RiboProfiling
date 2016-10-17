#' Plots the summarized coverage in a specified range (e.g. around TSS)
#' for the specified match sizes
#'
#' @param covSummarized a list of GRanges objects.
#' For each matchSize a GRanges object of the summarized coverage.
#' @examples
#' #read the BAM file into a GAlignments object using
#' #GenomicAlignments::readGAlignments
#' #the GAlignments object should be similar to ctrlGAlignments
#' data(ctrlGAlignments)
#' aln <- ctrlGAlignments
#' #transform the GAlignments object into a GRanges object (faster processing)
#' alnGRanges <- readsToStartOrEnd(aln, what="start")
#' #make a txdb object containing the annotations for the specified species.
#' #In this case hg19.
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' #Please make sure that seqnames of txdb correspond to
#' #the seqnames of the alignment files ("chr" particle)
#' #if not rename the txdb seqlevels
#' #renameSeqlevels(txdb, sub("chr", "",seqlevels(txdb)))
#' #get the flanking region around the promoter of the best expressed CDSs
#' oneBinRanges <- aroundPromoter(txdb, alnGRanges)
#' #the read start coverage around the TSS as a percentage for all match sizes.
#' covSummarized <- readStartCov(alnGRanges, oneBinRanges, matchSize="all",
#' c(-20,20), "aroundTSS", charPerc="perc")
#' trackPlotTSS <- plotSummarizedCov(covSummarized)
#' print(trackPlotTSS)
#' @return a ggplot2 plot of read coverage in interval
#' @export
#' @importFrom ggbio tracks
#' @importFrom grid unit
#'

plotSummarizedCov <-
    function(
        covSummarized){

    if(!inherits(covSummarized, "list")){
        stop("The covSummarized object is not a list!\n")
    }
    else{
        if(!is(covSummarized[[1]], "GRanges")){
            stop("The covSummarized object is not a list of GRanges objects!\n")
        }
    }

    listPlotSum <- lapply(covSummarized, function(iSumCov){
        maxPeak <- max(iSumCov$values);
        maxPeakPos <- start(iSumCov)[which(iSumCov$values == maxPeak)];
        if(maxPeak <= 100){yLab <- "% of reads"}
        else{yLab <- "nbr reads"}
        iPlot <-
            ggplot(iSumCov, ggplot2::aes(start, values)) +
            ggplot2::geom_point() +
            ggplot2::geom_line() +
            ggplot2::xlab("") +
            ggplot2::ylab(yLab);
        iPlot <-
            iPlot +
            ggplot2::geom_point(
                data=data.frame(start=maxPeakPos, values=maxPeak),
                ggplot2::aes(x=start, y=values),
                colour="brown2",
                size=3
            );
        iPlot <-
            iPlot +
            ggplot2::geom_text(
                data=data.frame(start=maxPeakPos, values=maxPeak),
                ggplot2::aes(label=start),
                colour="brown2",
                hjust=-1.5,
                fontface="bold"
            )
    })
    names(listPlotSum) <- names(covSummarized)

    #try to adapt the size of the track and that of the image
    trackPlotTSS <- ggbio::tracks(
        listPlotSum,
        heights=rep(grid::unit(4, "cm"), length(listPlotSum)),
        main.height=grid::unit(1.75 * length(listPlotSum), "cm")
    )

    return(trackPlotTSS)

}
