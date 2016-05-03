#' Starting from a BAM file path: quality plots, shift ribosome position,
#' coverage on multiple transcript features and on codons.
#'
#'
#' @importMethodsFrom GenomeInfoDb "seqlevels<-"
#' @importMethodsFrom S4Vectors "runValue<-"
#' @importFrom Rsamtools BamFile
#' @param listeInputBamFile A character path or a vector of paths to the
#'   ribo-seq BAM file(s). If multiple BAM files are provided, they should come
#'   from the same genome alignment.
#' @param paramScanBAM NULL or ScanBamParam object specifying what fields and
#' which records are imported. Default value is NULL.
#' @param genomeName a character object containing the name of the genome used
#'   for the alignment BAM file. The name should be one of the UCSC ensGene
#'   list: ucscGenomes()[ , "db"]. Ex. "hg19" or "mm10". This parameter is used
#'   to build the TxDb object.
#' @param txdb a TxDb object containing the annotations info to comfront with
#'   the alignment files. Either genomeName or txdb parameters should be
#'   provided.
#' @param percBestExpressed numeric [between 0 and 1]. The percentage of best
#'   expressed CDSs on which to plot the coverage around the TSS. Necessary if
#'   the shiftValue parameter must be estimated. Default value 0.03 (3\%).
#' @param flankSize an integer. How many bp left and right of the TSS should the
#'   coverage be performed?
#' @param offsetStartEnd a character object. Either "start" (the default) or "end"
#' for read start or read end to define the offset.
#' @param listShiftValue a vector of integer. It should have the same length as
#'   the inputBamFile vector. The numeric value for shifting ranges of reads on
#'   genomic features when computing coverage. Set this parameter to 0 if no
#'   shift should be performed. If this parameter is missing, the shiftValue is
#'   computed based on the maximum peak of read start coverage around the TSS. A
#'   plot is produced to illustrate this estimation.
#' @return A list of list for each BAM file in the inputBamFile list. For each
#'   BAM file 2 objects are returned: one data.frame with info on the genomic
#'   features and the corresponding coverage column, and one list of per ORF
#'   codon coverage.
#' @examples
#' #the txdb object can be given as parameter or not.
#' #If it is not specified, a txdb object is build from UCSC.
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' #in this example only one BAM file is treated.
#' #However, multiple BAM files can be analyzed together.
#' myFile <- system.file("extdata", "ctrl_sample.bam", package="RiboProfiling")
#' listeInputBam <- c(myFile)
#'
#' #when running this function it is important that chromosome names
#' #in UCSC and your BAM correspond: the "chr" particle
#' covData <- riboSeqFromBAM(listeInputBam, txdb=txdb, listShiftValue=c(-14))
#' @export
#' @import GenomicAlignments
#' @import rtracklayer
#'
riboSeqFromBAM <-
    function(
        listeInputBamFile,
        paramScanBAM,
        genomeName,
        txdb,
        percBestExpressed,
        flankSize,
        offsetStartEnd,
        listShiftValue){

    listCounts <- list(length=length(listeInputBamFile))
    listInfo <- list(length=length(listeInputBamFile))
    listCodonCounts <- list(length=length(listeInputBamFile))

    #check the offsetStartEnd parameter validity
    if(missing(offsetStartEnd)){
        offsetStartEnd <- "start"
    }
    if((offsetStartEnd != "start" && offsetStartEnd != "end")){
        warning("offsetStartEnd parameter is invalid. Set to default start value!\n")
        offsetStartEnd <- "start"
    }

    #check for paramScanBAM object
    if(missing(paramScanBAM)){
        paramScanBAM <- NULL
    }
    if(!is(paramScanBAM, "ScanBamParam")){
        warning("paramScanBAM parameter is not a ScanBamParam object. Set to default NULL value!\n")
        paramScanBAM <- NULL
    }

    if(missing(txdb)){
        message("\nGet UCSC ensGene annotations.\n")
        if(missing(genomeName)){
            stop("Please provide either a txdb or a genomeName parameter!\n")
        }
        else{
            ####### get the UCSC ensGene annotations
            #for the annotations, check the available genomes
            listUCSCAvGenomes <- rtracklayer::ucscGenomes()[ , "db"]
            if(!(genomeName %in% listUCSCAvGenomes)){
                stop(
                    paste(
                        "The genome ",
                        genomeName,
                        " is not present in the UCSC database! Available genomes: \n",
                        paste(listUCSCAvGenomes, collapse="\t"),
                        "\n",sep=""
                    )
                )
            }
            #gene annotations from UCSC
            txdb <- suppressWarnings(
                GenomicFeatures::makeTxDbFromUCSC(
                    genome=genomeName,
                    tablename="ensGene",
                    url="http://genome-euro.ucsc.edu/cgi-bin/"
                ))
        }

    }

    #first get only protein_coding transcripts
    #listTranscriptsPerGene=transcriptsBy(txdb,by="gene")
    #grouped by transcript
    cds <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
    #it groups all described exons by gene
    cdsAll <- GenomicFeatures::cdsBy(txdb, by="gene")
    exonGRanges <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)

    cdsPosTransc <- orfRelativePos(cds, exonGRanges)

    if(missing(listShiftValue)){
        globBoolShift <- 0
    }else{
        if(length(listeInputBamFile) != length(listShiftValue)){
            stop("listShiftValue should have the same length as listInputBamFile!\n")
        }
        else{
            globBoolShift <- 1
        }
    }
    namesInputBF <- vector(length=length(listeInputBamFile))

    #BAM file treatment
    for(ixFile in seq_len(length(listeInputBamFile))){
        if(inherits(listeInputBamFile[[ixFile]], "character")){
            message(paste(listeInputBamFile[[ixFile]], "\n", sep=""))
            inputBamFile <- listeInputBamFile[[ixFile]]
        }
        if(is(listeInputBamFile[[ixFile]], "BamFile")){
            message(paste(listeInputBamFile[[ixFile]]$path,"\n",sep=""))
            inputBamFile=paste(listeInputBamFile[[ixFile]]$path,"",sep="")
        }

        tmpName <- unlist(strsplit(inputBamFile, "/"))
        namesInputBF[ixFile] <- tmpName[length(tmpName)]

        #read the BAM file
        message("\nRead alignment file\n")
        aln <- GenomicAlignments::readGAlignments(inputBamFile, param=paramScanBAM)
        #transform the GAlignments vs GRanges (containg the read start info)
        #also add the info on the match size of the read
        alnGRanges <- readsToStartOrEnd(aln, what=offsetStartEnd)

        if(length(unique(GenomeInfoDb::seqnames(alnGRanges))) <= 0){
            stop("Coverage is null for all chromosomes. Check seqnames!\n")
        }

        message("\nPlot match length distribution.\n")
        #plot the match length distribution in the alignment file in number of reads
        dataGGPMLength <- histMatchLength(aln, 0, namesInputBF[ixFile])
        dataGGPMLength[[2]]
        print(dataGGPMLength[[2]])

        #if txdb and the BAM have differences in terms of "chr" name annotation
        seqLvlsInters <- length(
            intersect(
                GenomeInfoDb::seqlevels(aln),
                GenomeInfoDb::seqlevels(txdb)
            ))
        myCond1 <-
            seqLvlsInters / length(S4Vectors::runValue(GenomeInfoDb::seqnames(aln))) * 100 <= 80
        myCond2 <-
            seqLvlsInters / length(GenomeInfoDb::seqlevels(txdb)) * 100 <= 80
        if(myCond1 && myCond2){
            warning("Differences in seqlevels between the txdb and the BAM! \n")
            message(
                paste(
                    "# The BAM file has the following seqlevels: \n",
                    paste(GenomeInfoDb::seqlevels(aln), sep="\t"),
                    "\n",
                    sep="")
            )
            message(
                paste(
                    "# The transcript database file has the following seqlevels: \n",
                    paste(GenomeInfoDb::seqlevels(txdb), sep="\t"),
                    "\n",
                    sep="")
            )
        }


        message("\nPlot read start coverage around TSS.\n")
        #first get the chromosomes present in the alignment file
        seqInAlignment <- unique(seqnames(alnGRanges))
        if(length(seqInAlignment) <= 0){
            stop("Coverage is null for all chromosomes. Check seqnames!\n")
        }

        #keep only CDS on the chromosomes for which there are reads
        #cdsAll_chr=cdsAll[seqnames(cdsAll) %in% seqInAlignment]
        seqlevels(cdsAll, force=TRUE) <- as.character(seqInAlignment)
        cdsByGene <- cdsAll[S4Vectors::elementNROWS(cdsAll) != 0]


        ####### plot the coverage x bp left and right from the TSS
        #get coverage on CDS
        countsPCGenesAllExons <-
            GenomicAlignments::summarizeOverlaps(
                cdsByGene,
                alnGRanges,
                mode="IntersectionNotEmpty"
            )
        #then select the percBestExpressed\% best expressed genes
        if(missing(percBestExpressed)){
            percBestExpressed <- 0.03
        }
        if(!inherits(percBestExpressed, "numeric") || percBestExpressed < 0){
            warning("Invalid percBestExpressed parameter. Default 0.03!\n")
            percBestExpressed <- 0.03
        }

        vecCountsPerGene <- SummarizedExperiment::assays(countsPCGenesAllExons)$counts
        quantCounts <- quantile(vecCountsPerGene, 1-percBestExpressed)
        if(quantCounts <= 0 || missing(quantCounts)){
            stop("No gene had counts overlapping the CDS!\n")
        }

        indexGenesBestExpressed <- which(vecCountsPerGene >= quantCounts)
        geneNamesBestExpressed <-
            rownames(countsPCGenesAllExons)[indexGenesBestExpressed]

        #coverage n bases left and right from the start codon
        cdsByBestExpressed <-
            cdsByGene[names(cdsByGene) %in% geneNamesBestExpressed]

        #choose the longest transcript per gene (very rapid)
        maxWidthTranscBestExpGene <- max(width(cdsByBestExpressed))
        cdsByBestExprLongTransc <- cdsByBestExpressed[width(cdsByBestExpressed) == maxWidthTranscBestExpGene, ]
        #range around the promoters of the best expressed, longest transcripts
        if(missing(flankSize)){
            flankSize <- 20
        }
        if(!inherits(flankSize, "numeric") || flankSize < 0){
            warning("Invalid flankSize parameter. Default 20!\n")
            flankSize <- 20
        }

        flankPromoter <-
            GenomicRanges::promoters(
                cdsByBestExprLongTransc,
                flankSize,
                flankSize + 1
            )
        #keep only unique TSS regions
        #for those genes with multiple cds per gene keep only the first
        flankPromoterUniq <-
            S4Vectors::endoapply(flankPromoter, function(ixTSS){ ixTSS[1] })

        #transform the gene ranges into one GRanges object and add cds_id.
        oneBinRanges <-
            unlist(endoapply(
                flankPromoterUniq,
                function(ixTSS){
                    tmpOneBinGRanges=unlist(tile(ixTSS, n=width(ixTSS)));
                    values(tmpOneBinGRanges)=S4Vectors::DataFrame(
                        values=0,
                        idSeq=rep(ixTSS$cds_id[1], length(tmpOneBinGRanges))
                    );
                    tmpOneBinGRanges}))

        #get a list of coverage in the TSS flanking region
        listCovSummarized <-
            readStartCov(
                alnGRanges,
                oneBinRanges,
                matchSize="all",
                c(-flankSize, flankSize),
                "aroundTSS",
                charPerc="perc"
            )
        #plot the coverage of read starts around the TSS +- flankSize
        trackPlotTSS <- plotSummarizedCov(listCovSummarized)
        print(trackPlotTSS)
        #if no shiftValue is defined it is estimated as the position of the maximum
        #coverage value around the TSS
        boolShift <- 0
        if(globBoolShift == 0){
            while(boolShift == 0){
                shiftValue <- start(listCovSummarized[[1]])[which(listCovSummarized[[1]]$values == max(listCovSummarized[[1]]$values))] - 1
                message(
                    paste(
                        "# The following shift value will be applied: ",
                        shiftValue,
                        "\n",
                        sep="")
                )
                boolShift <- 1
            }

        }
        else{
            shiftValue <- listShiftValue[ixFile]
        }

        message("\nCount read start coverage on the CDS, 5pUTR, and 3pUTR.\n")
        ####### count the read start coverage on the different genomic features
        listInfo[[ixFile]] <-
            countShiftReads(
                exonGRanges[names(cdsPosTransc)],
                cdsPosTransc,
                alnGRanges,
                shiftValue
            )
        listCounts[[ixFile]] <- listInfo[[ixFile]][[1]]
        listCodonCounts[[ixFile]] <- listInfo[[ixFile]][[2]]

    }

    names(listCounts) <- namesInputBF
    names(listInfo) <- namesInputBF
    names(listCodonCounts) <- namesInputBF

    listCountsPlots <-
        countsPlot(listCounts, grep("_counts$", colnames(listCounts[[1]])), 1)
    invisible(utils::capture.output(
        suppressWarnings(
            suppressMessages(
                print(listCountsPlots))
            )
        )
        )

    return(listInfo)
}
