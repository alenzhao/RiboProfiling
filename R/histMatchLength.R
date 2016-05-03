#' Histogram of match length distribution of reads.
#'
#' @param aln A GAlignments object of the BAM mapping file.
#' @param log10Transf A boolean. Either 0 (default) or 1 (log10).
#' @param titleHist a character. The main title for the histogram.
#' Default - none.
#' @return A list with 2 elements.
#' The first element:
#' a data.frame of the number of counts per match length distribution.
#' The second element in the list:
#' a ggplot2 histogram of the match length distribution.
#' @examples
#' #starting from a GAlignment object
#' data(ctrlGAlignments)
#' aln <- ctrlGAlignments
#' #no log10 scaling
#' matchLenDistr <- histMatchLength(aln, 0)
#' #to plot the histogram
#' matchLenDistr[[2]]
#' @export
#'
histMatchLength <-
    function(
        aln,
        log10Transf=0,
        titleHist){
    # CHECK if the objects are the correct class
    if(missing(log10Transf)){
        log10Transf == 0
    }
    stopifnot(log10Transf == 0 || log10Transf == 1)
    stopifnot(is(aln, "GAlignments"))

    if(missing(titleHist)){
        titleHist <- ""
    }
    else{
        if(titleHist != ""){
            tmpName <- unlist(strsplit(titleHist, "/"))
            titleHist <- tmpName[length(tmpName)]
        }
    }

    matchSize <- NULL

    #get match length
    vecMLength <-
        GenomicAlignments::cigarWidthAlongReferenceSpace(
            GenomicAlignments::cigar(aln),
            N.regions.removed=TRUE
        )

    #contingency table
    countsMLength <- table(vecMLength)

    #ggplot2 dataFrame format
    dataGGPMLength <- cbind(names(countsMLength), countsMLength)
    dataGGPMLength <- data.frame(dataGGPMLength)
    colnames(dataGGPMLength) <- c("matchSize", "counts")
    dataGGPMLength$counts <- as.numeric(as.character(dataGGPMLength$counts))

    if(log10Transf == 0){
        #barplot of nbr of reads per match size
        barPlotMatchSize <-
            ggplot2::ggplot(dataGGPMLength, ggplot2::aes(matchSize, counts)) +
            ggplot2::geom_bar(
                stat="identity",
                position="dodge",
                color="steelblue",
                fill="steelblue"
            ) +
            ggplot2::labs(x="Match size", y="Number of reads",fill=NULL)+
            ggplot2::theme(
                axis.title.y=ggplot2::element_text(
                    family="serif",
                    size=15,
                    face="bold"
                )
            ) +
            ggplot2::theme(
                axis.title.x=ggplot2::element_text(
                    family="serif",
                    size=15,
                    face="bold"
                )
            ) +
            ggplot2::ggtitle(titleHist) +
            ggplot2::theme(
                axis.text.x=ggplot2::element_text(angle=90, hjust=1)
            )
    }else{
        #barplot of nbr of reads per match size
        barPlotMatchSize <-
            ggplot2::ggplot(
                dataGGPMLength,
                ggplot2::aes(matchSize, log10(counts))
            ) +
            ggplot2::geom_bar(
                stat="identity",
                position="dodge",
                color="steelblue",
                fill="steelblue"
            ) +
            ggplot2::labs(
                x="Match size",
                y="log10 Number of reads",
                fill=NULL
            ) +
            ggplot2::theme(
                axis.title.y=ggplot2::element_text(
                    family="serif",
                    size=15,
                    face="bold"
                )
            ) +
            ggplot2::theme(
                axis.title.x=ggplot2::element_text(
                    family="serif",
                    size=15,
                    face="bold"
                )
            ) +
            ggplot2::ggtitle(titleHist)+
            ggplot2::theme(
                axis.text.x=ggplot2::element_text(angle=90, hjust=1)
            )
    }


    return(list(dataGGPMLength,barPlotMatchSize))

}
