
####### redefine the start/ end values of a GRanges object
# rangeObj  : a GRanges object.
# newValues : a vector of the same size as rangeObj
# startOrEnd: which values should be replaced: either "start", "end", or "both"
redefineStartEnd <-
    function(
        rangeObj,
        newValues,
        startOrEnd){
    if(length(rangeObj) != length(newValues)){
        stop("The rangeObj and your new values don't have the same length!\n")
    }
    else{
        if(startOrEnd == "start" || startOrEnd == "both"){
            start(rangeObj) <- newValues
        }
        if(startOrEnd == "end" || startOrEnd == "both"){
            end(rangeObj) <- newValues
        }
    }
    return(rangeObj)
}


readStartCov1Aln <-
    function(
        alnGRanges,
        oneBinRanges,
        fixedInterval,
        renameChr,
        charPerc){
    oneBinReadStartCov <- oneBinRanges
    vecReadStartCov <-
        suppressWarnings(findOverlaps(oneBinRanges, alnGRanges))
    values(oneBinReadStartCov) =
        S4Vectors::DataFrame(
            values=IRanges::as.table(vecReadStartCov),
            idSeq=oneBinRanges$idSeq
        )
    #regroup per gene name into one object:GRangesList
    reduceRSC <-
        reduce(oneBinReadStartCov, with.revmap=TRUE)
    revmap <-
        mcols(reduceRSC)$revmap
    listReadStartCov <-
        relist(oneBinReadStartCov[unlist(revmap)], revmap)

    #normalize the previous ranges to the specified interval: 0-100 or 0 center
    #for the reverse strand reverse the fixedInterval extremes
    ixReverse <- which(unlist(runValue(strand(listReadStartCov))) == "-")
    if(length(ixReverse) > 0){
        listReadStartCovNormPlus <-
            normRange(listReadStartCov[-ixReverse], fixedInterval)
        listReadStartCovNormRev <-
            normRange(listReadStartCov[ixReverse], -1 * fixedInterval)
        listReadStartCovNorm <-
            c(listReadStartCovNormPlus, listReadStartCovNormRev)
    }
    else{
        listReadStartCovNorm <-
            normRange(listReadStartCov[ixReverse], fixedInterval)
    }

    #get the summarized coverage
    covSummarized <- getSummarizedCov(
        listReadStartCovNorm,
        charPerc,
        seqNameRedef=renameChr)

    return(covSummarized)

}


####### normalizes the GRanges sizes to the fixed interval size
####### and report the associated counts values
# listTSSReadStartCov : a list of GRanges object.
#                       : it contains the nbr of counts in 1bin ranges
#                       : example:
# listTSSReadStartCov <- computeReadStartCov(posByChrom,flankPromoterUniq)
# fixedInterval         : a vector containing the min and max of the interval
#                       : the range of the fixed interval
#                       : to which all the ranges should be normalized.
#                       : example: c(0,100) or c(-20,20)
#                       : Warning! a non-integer interval of a 0-1 interval
#                       : will be automatically transformed to an integer range.
#                       : Instead of 0-1, use 0-100.
# renameChr             : a character object
#                       : it provides a name for the fixed interval.
#                       : The new unique seqnames.
#                       : example: CDSLength or aroundTSS
#it returns a list of GRanges objects
normRange <-
    function(
        listTSSReadStartCov,
        fixedInterval){
    iStart <- fixedInterval[1]
    iEnd <- fixedInterval[2]

    #check that the length of ranges in the listTSSReadStartCov object
    #is the same for all objects
    testLength <-
        sapply(listTSSReadStartCov, length)
    testLengthList <- unique(testLength)
    if(length(testLengthList) != 1){
        ixLengthPb <- which(testLength != length(iStart:iEnd))
        if(length(ixLengthPb)/length(testLength) >= 0.7){
            stop(paste("More than 70% of the coding sequences (idSeq) have",
                       "multiple gene annontations!\n", sep=""))
        }
        else{
            warning("idSeqs with multiple gene ids have been eliminated!\n")
            listTSSReadStartCov <- listTSSReadStartCov[-ixLengthPb]
            if(length(iStart:iEnd) <=1 ){
                warning("The fixedInterval parameter is not correct. \n")
            }

            #the new ranges would start at iStart, end at iEnd and have a length
            #equal to length of the old range here I take advantage that my old
            #ranges is made up of 1bp bins
            newIntervalPos <-
                seq(iStart, iEnd, length.out=length(listTSSReadStartCov[[1]]))

            listTSSReadStartCovRedefined <-
                sapply(
                    listTSSReadStartCov,
                    redefineStartEnd,
                    newValues=newIntervalPos,
                    startOrEnd="both"
                )
            return(listTSSReadStartCovRedefined)
        }

    }
    else{
        #and the same length as the fixed interval
        if(testLengthList != length(iStart:iEnd)){
            stop("The fixedInterval has different length from the GRanges objects!\n")
        }
        else{
            if(length(iStart:iEnd) <= 1){
                warning("The fixedInterval parameter is not correct. \n")
            }

            #the new ranges would start at iStart, end at iEnd and have a length
            #equal to length of the old range here I take advantage that my old
            #ranges is made up of 1bp bins
            newIntervalPos <-
                seq(
                    iStart,
                    iEnd,
                    length.out=length(listTSSReadStartCov[[1]])
                )

            listTSSReadStartCovRedefined <-
                sapply(
                    listTSSReadStartCov,
                    redefineStartEnd,
                    newValues=newIntervalPos,
                    startOrEnd="both"
                )
            return(listTSSReadStartCovRedefined)
        }
    }


}

####### apply the reference rage to the summarized coverage of all GRanges
# listCovNorm          : a list of GRanges object.
#                       : list of read counts in exactly the same 1bin ranges
#                       : example:
# listTSSReadStartCov <- computeReadStartCov(posByChrom,flankPromoterUniq)
# listCovNorm=normRange(listTSSReadStartCov,c((-1) * 20, 20), "TSS")
# charPerc              : a string
#                       : "perc" (the default) or "sum"
#                       : example: charPerc <- "perc"
#it returns a GRanges object: the sum or the percentage of coverage per position
getSummarizedCov <-
    function(
        listCovNorm,
        charPerc="perc",
        seqNameRedef){
    #check that the length of ranges in all the objects of the list
    testLengthList <- unique(sapply(listCovNorm, length))
    if(length(testLengthList) != 1){
        stop("Objects in the listCovNorm list have the different ranges!\n")
    }

    if(charPerc != "perc" && charPerc != "sum"){
        warning("charPerc parameter invalid: default value = perc!\n")
        charPerc <- "perc"
    }

    #now get the sum and the % of read counts per each position in the range
    #check the number of objects found.
    #If less than 10 send a warning message
    nbrObj <- length(listCovNorm)
    if(nbrObj < 1){
        stop("No range has been provided for coverage plot!\n")
    }
    else{
        if(nbrObj >= 1 && nbrObj < 10){
            warning(paste("Only ", nbrObj, " ranges analyzed!\n", sep=""))
        }
        if(seqNameRedef == "" || missing(seqNameRedef) || is.na(seqNameRedef)){
            warning("Invalid seqNameRedef parameter -> No renaming.\n")
        }

        #sum the read Start counts over all the TSS
        listReadStartCov <- lapply(listCovNorm,
                                 function(TSSCounts){
                                     TSSCounts = sort(TSSCounts);
                                     TSSCounts$values
                                 })
        sumReadStartPerPos <- Reduce("+", listReadStartCov)
        if(charPerc == "perc"){
            percReadStartPerPos <-
                sumReadStartPerPos/sum(sumReadStartPerPos)*100
            TSSFlankPercReadStart <-
                GRanges(
                    seqnames=seqNameRedef,
                    ranges(sort(listCovNorm[[1]])),
                    values=percReadStartPerPos
                )
        }
        else{
            TSSFlankPercReadStart <-
                GRanges(
                    seqnames=seqNameRedef,
                    ranges(sort(listCovNorm[[1]])),
                    values=sumReadStartPerPos
                )
        }


    }

    return(TSSFlankPercReadStart)

}

####### prepare the feature object and shift it accordingly
# transcGRangesList    : a GRangesList object.
#                       : it contains the GRanges of exons for each transcript
#                       : example:
# transcGRangesList=exonsBy(txdb,by="tx",use.names=T)
# shiftValue            : numeric
#                       : offset of read starts on transcript.
#                       : output of plotSummarizedCov function
#                       : example: shift=-14
# returns a list of integers, for transcripts sizes superior to the shiftValue
applyShiftFeature <-
    function(
        transcGRangesList,
        shiftValue){
    #check parameter validity
    if(missing(shiftValue) || !inherits(shiftValue, "numeric")){
        shiftValue <- 0
        warning("Incorrect shiftValue parameter! No shift is performed!\n")
    }
    if(!is(transcGRangesList, "GRangesList")){
        stop(
            paste(
                "transcGRangesList parameter is of class ",
                class(transcGRangesList),
                " instead of GRangesList!\n",
                sep=""
            )
        )
    }

    transcWidth <-
        GenomicFeatures::transcriptWidths(
            start(transcGRangesList),
            end(transcGRangesList)
        )
    absShiftVal <- abs(shiftValue)

    #if width of transcript is smaller than the absolute shiftValue
    #eliminate the transcript
    ixSmallTransc <- which(transcWidth <= absShiftVal)
    if(length(ixSmallTransc) > 0){
        transcGRangesList <- transcGRangesList[-ixSmallTransc]
        transcWidth <-
            GenomicFeatures::transcriptWidths(
                start(transcGRangesList),
                end(transcGRangesList)
            )
    }

    #now if the shiftValue is positive, the start of the transcript is shifted
    if(shiftValue > 0){
        usefulRangeOnTransc <-
            cbind(
                startT=rep(absShiftVal+1, length(transcGRangesList)),
                endT=transcWidth
            )
    }
    #else, it is the end of the transcript that we shift
    else{
        usefulRangeOnTransc <- cbind(startT=1, endT=transcWidth - absShiftVal)
    }
    listeUsefulRanges <-
        lapply(seq_len(length(transcGRangesList)),
               function(ixTransc){
                   usefulRangeOnTransc[ixTransc, 1]:usefulRangeOnTransc[ixTransc, 2]}
        )
    #for the remaining positions in the transcript, make 1bp bins of the genomic positions
    shiftedTransc <-
        GenomicFeatures::transcriptLocs2refLocs(
            listeUsefulRanges,
            start(transcGRangesList),
            end(transcGRangesList),
            as.character(S4Vectors::runValue(strand(transcGRangesList))),
            decreasing.rank.on.minus.strand=TRUE
        )
    names(shiftedTransc) <- names(transcGRangesList)
    return(shiftedTransc)

}



naTozeroRle=function(rleObject){
    #check rleObject class
    ixNA=which(is.na(S4Vectors::runValue(rleObject)))
    if(length(ixNA)>0)
    {
        S4Vectors::runValue(rleObject)[ixNA]=0;
        S4Vectors::runLength(rleObject)[ixNA]=0
    }
    return(rleObject)
}


#make the ggplot2 format
funcDataFrameGGPlot2 <-
    function(
        listData,
        columnIndex,
        logBool){
    dataGgplot <- data.frame()
    for(ixData in seq_len(length(listData))){
        data <- listData[[ixData]]
        for(i in seq_len(length(columnIndex))){
            if(logBool == 1){
                value <-
                    log2(as.numeric(as.character(data[,columnIndex[i]])) + 1)
            }
            else{
                value <- as.numeric(as.character(data[,columnIndex[i]]))
            }
            #if no names for the list
            if(is.null(names(listData)[ixData]) || is.na(names(listData)[ixData])){
                sampleName <- paste("data", ixData, sep="")
            }
            else{
                sampleName <- names(listData)[ixData]
            }
            tmpData <-
                cbind(
                    rep(sampleName, nrow(data)),
                    rep(colnames(data)[columnIndex[i]], nrow(data)),
                    value
                )
            if(ixData == 1 && i == 1){
                dataGgplot <- tmpData
            }
            else{
                dataGgplot <- rbind(dataGgplot, tmpData)
            }
        }
    }

    colnames(dataGgplot) <- c("sample", "type", "value")
    dataGgplot <- as.data.frame(dataGgplot)
    dataGgplot$value <- as.numeric(as.character(dataGgplot$value))

    return(dataGgplot)
}


#ggplot boxplot
funcBoxplot <-
    function(
        dataGgplot,
        xText,
        yText,
        titleText){
    type=NULL
    value=NULL
    x=NULL


    pBoxplot <- ggplot2::ggplot(data = dataGgplot, ggplot2::aes(x=type, y=value))
    if(length(which(colnames(dataGgplot) == "sample")) == 1){
        pBoxplot <- pBoxplot +
            ggplot2::geom_boxplot(ggplot2::aes(fill = sample))
        pBoxplot <- pBoxplot +
            ggplot2::geom_point(
                ggplot2::aes(y = value, group = sample),
                position=ggplot2::position_dodge(width=0.75)
            )
    }
    else{
        pBoxplot <- pBoxplot +
            ggplot2::geom_boxplot(ggplot2::aes(fill = type))
        pBoxplot <- pBoxplot +
            ggplot2::geom_point(
                ggplot2::aes(y = value, group = type),
                position =ggplot2::position_dodge(width=0.75)
            )
    }
    pBoxplot <- pBoxplot +
        ggplot2::xlab(xText) + ggplot2::ylab(yText) +
        ggplot2::ggtitle(titleText)
    pBoxplot <- pBoxplot +
        ggplot2::guides(fill=ggplot2::guide_legend(title="sample")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, hjust=1))
    pBoxplot
    return(pBoxplot)
}


#http://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
makePairs <-
    function(data){
    x <- NULL
    y <- NULL

    grid <- expand.grid(x=1:ncol(data), y=1:ncol(data))
    gridUpper <- subset(grid, x < y)
    upper <- do.call("rbind", lapply(seq_len(NROW(gridUpper)), function(i){
        xcol <- gridUpper[i, "x"]
        ycol <- gridUpper[i, "y"]
        data.frame(xvar=names(data)[ycol], yvar=names(data)[xcol],
                   x=data[,xcol], y=data[,ycol], data)
    }))
    upper$xvar <- factor(upper$xvar, levels=names(data))
    upper$yvar <- factor(upper$yvar, levels=names(data))

    gridLower<- subset(grid,x>y)
    lower <- do.call("rbind", lapply(seq_len(NROW(gridLower)), function(i){
        xcol <- gridLower[i, "x"]
        ycol <- gridLower[i, "y"]
        data.frame(
            xvar=names(data)[ycol],
            yvar=names(data)[xcol],
            x=(max(data[,xcol])+min(data[,xcol]))/2,
            y=(max(data[,ycol])+min(data[,ycol]))/2,
            labs=round(cor(data[,xcol],data[,ycol]),3)
        )
    }))


    densities <- do.call("rbind", lapply(1:ncol(data), function(i){
        data.frame(xvar=names(data)[i], yvar=names(data)[i], x=data[, i])
    }))
    list(upper=upper, lower=lower, densities=densities)
}

#this function is an adaptation of the ggplot pairs plot
#it takes as input a matrix or data.frame and makes the pairs ggplot
funcPlotPairs <-
    function(
        data,
        title){
    ..scaled.. <- NULL
    x <- NULL
    y <- NULL
    labs <- NULL

    gridData <- makePairs(data)

    upperGrid <- data.frame(gridData$upper)

    pPairs <- ggplot2::ggplot(upperGrid, ggplot2::aes_string(x="x", y="y")) +
        ggplot2::facet_grid(xvar ~ yvar, scales="free") +
        ggplot2::geom_point(color="#6495ED") +
        ggplot2::stat_density(
            ggplot2::aes(x=x, y=..scaled.. * diff(range(x)) + min(x)),
            data=gridData$densities, position="identity",
            colour="grey20", geom="line", lwd=1) +
        ggplot2::geom_abline(
            ggplot2::aes(color="Diagonal", fill="Diagonal", intercept = 0),
            lwd=1) +
        ggplot2::geom_smooth(
            ggplot2::aes(x=x, y=y, color = "Linear_regression", fill = "Linear_regression"),
            method = "lm", lwd=1) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title.x=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank()) +
        ggplot2::geom_text(
            data=gridData$lower,
            ggplot2::aes(x=x, y=y, label=labs, group=NULL),
            fontface=2,
            colour="steelblue4") +
        ggplot2::scale_fill_manual(
            name="Lines",
            values=c('Diagonal'='white', 'Linear_regression'='white')) +
        ggplot2::scale_colour_manual(
            name="Lines",
            values=c('Diagonal'='brown2', 'Linear_regression'='black'),
            guide='legend') +
        ggplot2::guides(
            colour=ggplot2::guide_legend(override.aes=list(linetype=c(1, 1)))) +
        ggplot2::ggtitle(title)

    return(pPairs)
}



getCodons <-
    function(
        seqChar){
    oneCharSplit <- sapply(seqChar, strsplit, "")[[1]]
    codonSeq <-
        paste0(
            oneCharSplit[c(TRUE, FALSE, FALSE)],
            oneCharSplit[c(FALSE, TRUE, FALSE)],
            oneCharSplit[c(FALSE, FALSE, TRUE)]
        )
    codonsInSeq <- cbind(seq_len(length(codonSeq)), codonSeq)
    colnames(codonsInSeq) <- c("codonID", "codon")
    codonsInSeq[,1] <- as.numeric(as.character(codonsInSeq[, 1]))
    return(codonsInSeq)
}

