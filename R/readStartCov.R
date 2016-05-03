#' Read start coverage around the TSS on the predifined CDSs
#'
#' @param alnGRanges A GRanges object containing the alignment information.
#' In order to improve the performance transform the GAlignments BAM object
#' into a GRanges object containing cigar match size as metadata.
#' @param oneBinRanges A GRanges object.
#' Transform the gene GRangesList into one big GRanges object.
#' Add the info on the cds_id.
#' @param matchSize either "all" or a vector of read match sizes.
#' If matchSize <- "all", then all the reads are used to compute the coverage.
#' If the matchSize is a vector of read match sizes,
#' the summarized coverage is reported per match size and for the sum up.
#' @param fixedInterval a numeric vector with the extremities of the interval.
#' Ex. fixedInterval <- c(-20,20) or fixedInterval <- c(0,40) ...
#' @param renameChr a character object.
#' It contains the name to be given to the new summarized coverage interval.
#' Ex. renameChr <- "aroundTSS" the summarized region around the TSS.
#' @param charPerc a character object. Either "perc" (the default) or "sum"
#' for percentage of counts per position or the sum of counts per position.
#  Ex. charPerc <- "perc"
#' @return a list of GRanges objects (for each matchSize chosen).
#' It contains the summarized coverage for the specified read match sizes.
#' @examples
#' #read the BAM into a GAlignments object using
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
#' #renameSeqlevels(txdb, sub("chr", "", seqlevels(txdb)))
#'
#' #get the flanking region around the promoter of the best expressed CDSs
#' oneBinRanges <- aroundPromoter(txdb, alnGRanges, percBestExpressed=0.01)
#'
#' #the coverage in the TSS flanking region for the summarized read match sizes
#' listPromoterCov <- readStartCov(
#'      alnGRanges,
#'      oneBinRanges,
#'      matchSize="all",
#'      fixedInterval=c(-20, 20),
#'      renameChr="aroundTSS",
#'      charPerc="perc"
#' )
#' @export
#' @import GenomicRanges

readStartCov <-
    function(
        alnGRanges,
        oneBinRanges,
        matchSize="all",
        fixedInterval,
        renameChr,
        charPerc="perc"){
    if(!is(alnGRanges, "GRanges")){
        stop("The alignment object is not a GRanges object!\n")
    }
    if(!inherits(fixedInterval, "numeric")){
        stop("The parameter fixedInterval should be a numeric vector!\n")
    }
    #the fixedInterval parameter should be a vector of exactly 2 values
    #containing the extremities of the interval.
    if(length(fixedInterval) != 2){
        stop("The fixedInterval must contain only the 2 interval borders!\n")
    }

    if((charPerc != "perc" && charPerc != "sum") || missing(charPerc)){
        warning("charPerc parameter is invalid. Set to default perc value!\n")
        charPerc <- "perc"
    }
    if(!is(oneBinRanges, "GRanges")){
        stop("The 1bp region object (oneBinRanges) is not a GRanges object!\n")
    }
    lenFlankRegion <- unique(table(names(oneBinRanges)))
    if(length(lenFlankRegion) > 1){
        stop(
            paste(
                "The oneBinRanges object has multiple fixedInterval sizes ",
                lenFlankRegion, "!.\n", sep=" "
            )
        )
    }
    if(lenFlankRegion != fixedInterval[2] - fixedInterval[1] + 1){
        stop(
            paste(
                "The oneBinRanges parameter should contain 1bp ranges equal to:",
                fixedInterval[2] - fixedInterval[1] + 1,
                ".\n",
                sep=""
            )
        )
    }

    #get the size of the match. Also works for spliced alignments
    vecMLength <- alnGRanges$score
    rangeOfMatch <- sort(unique(vecMLength))

    if(!identical(matchSize, "all")){
        if(inherits(matchSize, "numeric") || inherits(matchSize, "integer")){
            #check if all demanded match sizes are found in the alignment object
            specDefMatch <- setdiff(matchSize, rangeOfMatch)
            if(length(specDefMatch) > 0){
                stop(
                    paste(
                        "The following defined match sizes: ",
                        paste(specDefMatch, collapse=", "),
                        " are not part of the specified alignment object! ","
                        Please select match sizes in this range: ",
                        paste(rangeOfMatch, collapse=", "),
                        "\n",
                        sep=""
                    )
                )
            }
            else{
                #first do a coverage for the summed match sizes
                ixMatchAln <- which(!is.na(match(alnGRanges$score, matchSize)))
                globCovSumm=
                    readStartCov1Aln(
                        alnGRanges[ixMatchAln],
                        oneBinRanges,
                        fixedInterval,
                        renameChr,
                        charPerc
                    )
                #for each match size substract the alignment subset
                listCovSumPerMatch <- lapply(matchSize, function(iMatch){
                    ixMatchAln=which(alnGRanges$score == iMatch);
                    tmpAln=alnGRanges[ixMatchAln];
                    covSummarized=
                        readStartCov1Aln(
                            tmpAln,
                            oneBinRanges,
                            fixedInterval,
                            renameChr,
                            charPerc
                        )
                })

                listCovSum <- append(list(globCovSumm), listCovSumPerMatch)
                names(listCovSum) <-
                    c("sumUp", paste(matchSize, "_match", sep=""))
            }
        }
        else{
            warning("The matchSize parameter is invalid. Default: all!\n")
            matchSize <- "all"
        }
    }
    if(identical(matchSize, "all")){
        covSummarized <-
            readStartCov1Aln(
                alnGRanges,
                oneBinRanges,
                fixedInterval,
                renameChr,
                charPerc
            )
        listCovSum <- list(covSummarized)
        names(listCovSum) <- "sumUp"
    }

    return(listCovSum)

}
