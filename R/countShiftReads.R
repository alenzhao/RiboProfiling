#' Apply an offset on the read start along the transcript and
#' returns the coverage on the 5pUTR, CDS, 3pUTR,
#' as well as a matrix of codon coverage per ORF.
#'
#' @param exonGRanges a GRangesList.
#' It contains the exon coordinates grouped by transcript.
#' @param cdsPosTransc a list.
#' It contains the relative positions of the start and end of the ORFs.
#' The transcript names in exonGRanges and cdsPosTransc should be the same.
#' @param alnGRanges A GRanges object containing the alignment information.
#' In order to improve the performance the GAlignments BAM object
#' should be transformed into a GRanges object with cigar match size metadata.
#' @param shiftValue integer.
#' The offset for recalibrating reads on transcripts when computing coverage.
#' The default value for this parameter is 0, no offset should be performed.
#' @param motifSize an integer. The number of nucleotides in each motif
#' on which to compute coverage and usage. Either 3, 6, or 9.
#' Default 3 nucleotides (codon).
#' @return a list with 2 objects.
#' The first object in the list is a data.frame containing:
#' information on ORFs (names, chromosomal position, length)
#' as well as the counts on the 5pUTR, CDS and 3pUTR once the offset is applied.
#' The second object in the list is a list in itself.
#' It contains: for each ORF in the cdsPosTransc,
#' for each codon the sum of read starts covering the 3 codon nucleotides.
#' For motifs of size 6 nucleotides, the motif coverage is computed only
#' for the first codon in the motif, considered as the codon in the P-site.
#' For motifs of size 9 nucleotides, the motif coverage is computed only
#' for the second codon in the motif, considered as the codon in the P-site.
#' This per codon coverage does not contain information on the codon type,
#' just its position in the ORF and its coverage.
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
#' #renameSeqlevels(txdb, sub("chr", "", seqlevels(txdb)))
#'
#' #get all CDSs by transcript
#' cds <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
#'
#' #get all exons by transcript
#' exonGRanges <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
#' #get the per transcript relative position of start and end codons
#' cdsPosTransc <- orfRelativePos(cds, exonGRanges)
#' #compute the counts on the different features after applying
#' #the specified shift value on the read start along the transcript
#' countsData <- countShiftReads(exonGRanges[names(cdsPosTransc)], cdsPosTransc,
#'            alnGRanges, -14)
#' @export

countShiftReads <-
    function(
        exonGRanges,
        cdsPosTransc,
        alnGRanges,
        shiftValue,
        motifSize){
    #check cdsPosTransc no NA values, no missing, all integer
    if(missing(cdsPosTransc)){
        stop("Missing cdsPosTransc parameter!\n")
    }
    if(length(exonGRanges) != length(cdsPosTransc)){
        stop("Different lengths for exonGRanges and cdsPosTransc parameters!\n")
    }
    myCondNA <- which(is.na(unlist(cdsPosTransc)) | is.null(unlist(cdsPosTransc)))
    if(length(myCondNA) > 0){
        stop("Non-null, non-NA values for the cdsPosTransc parameter!\n")
    }

    #check parameter validity
    if(missing(shiftValue) || !inherits(shiftValue, "numeric")){
        shiftValue <- 0
        warning("Incorrect shiftValue parameter! No shift is performed!\n")
    }

    if(!is(exonGRanges, "GRangesList")){
        stop(
            paste(
                "exonGRanges parameter is of class ",
                class(exonGRanges),
                " instead of GRangesList!\n",
                sep=""
            )
        )
    }
    if(!is(alnGRanges, "GRanges")){
        stop(
            paste(
                "alnGRanges parameter is of class ",
                class(alnGRanges),
                " instead of GRanges!\n",
                sep=""
            )
        )
    }

    if(missing(motifSize) ||
       !is(motifSize, "numeric") ||
       motifSize %% 1 != 0 ||
       motifSize <= 0 ||
       !(motifSize %in% c(3, 6, 9))){
        warning("Param motifSize should be an integer! Accepted values 3, 6 or 9. Default value is 3.\n")
        motifSize <- 3
    }


    exonGRangesRestrict <- exonGRanges[names(cdsPosTransc)]

    if(length(exonGRangesRestrict) <= 5){
        stop("Less than 5 common transcripts btw exonGRanges and cdsPosTransc!\n")
    }
    else{
        if(length(exonGRangesRestrict) <= 10){
            warning("Less than 10 common transcripts between exonGRanges and cdsPosTransc!\n")
        }

    }



    transcWidth <-
        GenomicFeatures::transcriptWidths(
            start(exonGRangesRestrict),
            end(exonGRangesRestrict)
        )
    absShiftVal <- abs(shiftValue)

    #if width of transcript is smaller than the absolute shiftValue,
    #eliminate the transcript
    ixSmallTransc <- which(transcWidth <= absShiftVal)
    if(length(ixSmallTransc) > 0){
        #put a stop if too many small
        transcBig <- exonGRangesRestrict[-ixSmallTransc]
        cdsPosTRanscBig <- cdsPosTransc[-ixSmallTransc]
    }
    else{
        transcBig <- exonGRangesRestrict
    }

    #find read start overlap on the unshifted transcripts
    overlapReads <- suppressWarnings(findOverlaps(alnGRanges, transcBig))
    startOverlapReads <-
        split(
            start(alnGRanges[queryHits(overlapReads)]),
            factor(subjectHits(overlapReads))
        )
    overlapReadsRle <- sapply(startOverlapReads, S4Vectors::Rle)

    #now keep only those transcripts that have reads on them
    #put a warning if nothing left
    transcWithReads <- transcBig[as.numeric(names(overlapReadsRle))]
    cdsPosTranscWithReads <- cdsPosTransc[as.numeric(names(overlapReadsRle))]
    newTranscWidth <-
        GenomicFeatures::transcriptWidths(
            start(transcWithReads),
            end(transcWithReads)
        )


    #apply the shiftValues
    cdsPosTranscShifted <- do.call(rbind, cdsPosTranscWithReads) + shiftValue

    #maximum between the transcript start (1) and the shifted CDS start and minimum between the transcript end and the shifted CDS end
    listeRangesCDS <-
        lapply(seq_len(NROW(cdsPosTranscShifted)),
               function(ixTransc){
                   max(1, cdsPosTranscShifted[ixTransc, 1]):min(newTranscWidth[ixTransc], cdsPosTranscShifted[ixTransc, 2])
                   })

    #for the shifted 5UTR
    listeRanges5UTR <-
        lapply(seq_len(NROW(cdsPosTranscShifted)),
               function(ixTransc){
                   if((cdsPosTranscShifted[ixTransc, 1] - 1) < 1){
                       0
                   }
                   else{
                       max(1, shiftValue):(cdsPosTranscShifted[ixTransc, 1] - 1)
                   }
               }
        )
    #for the shifted 3UTR
    listeRanges3UTR <-
        lapply(seq_len(NROW(cdsPosTranscShifted)),
               function(ixTransc){
                   if((cdsPosTranscShifted[ixTransc, 2] + 1) > newTranscWidth[ixTransc]){
                       0
                   }
                   else{
                       (cdsPosTranscShifted[ixTransc, 2] + 1):min(newTranscWidth[ixTransc], newTranscWidth[ixTransc] + shiftValue)
                   }
               }
        )


    #transform transcripts into 1bp bin ranges
    binTransc <- applyShiftFeature(transcWithReads, 0)

    #additional info
    strandInfo <- S4Vectors::runValue(strand(transcWithReads))

    #find reads mapped on the shifted transcripts
    shiftedTranscMatches <- lapply(seq_len(length(binTransc)), function(ixTransc){
        if(strandInfo[[ixTransc]] == "-"){
            binTranscVal <- sort(binTransc[[ixTransc]], decreasing=TRUE)
        }
        else{
            binTranscVal <- sort(binTransc[[ixTransc]])
        };
        matchedReadsTransc <-
            match(sort(overlapReadsRle[[ixTransc]]), binTranscVal);
        matchedReadsCDS <-
            naTozeroRle(match(matchedReadsTransc, listeRangesCDS[[ixTransc]]));
        matchedReads5UTR <-
            naTozeroRle(match(matchedReadsTransc, listeRanges5UTR[[ixTransc]]));
        matchedReads3UTR <-
            naTozeroRle(match(matchedReadsTransc, listeRanges3UTR[[ixTransc]]));
        #count the number of reads per codon
        if(length(matchedReadsCDS) > 0){
            # if(motifSize <= 3){
            allCodonCounts <-
                aggregate(
                    S4Vectors::runLength(matchedReadsCDS),
                    by=list(ceiling(S4Vectors::runValue(matchedReadsCDS) / 3)),
                    FUN=sum
                )
            if(motifSize <= 3){
                myCodonCounts <- allCodonCounts
            }
            else{
                if(motifSize == 6){
                    myCodonCounts <-
                        allCodonCounts[1:(nrow(allCodonCounts)-1), ]
                }
                else{
                    if(motifSize == 9){
                        myCodonCounts <-
                            allCodonCounts[2:(nrow(allCodonCounts)-1), ]
                    }
                }
            }
#                 myCodonCounts <-
#                     aggregate(
#                         S4Vectors::runLength(matchedReadsCDS),
#                         by=list(ceiling(S4Vectors::runValue(matchedReadsCDS) / motifSize)),
#                         FUN=sum
#                     )
#             }
#             #for motifs >3 make overlapping windows shifted with 3 nucleotides
#             else{
#                 motifStart <- seq.int(
#                     from=1,
#                     to=max(runValue(matchedReadsCDS)),
#                     by=3)
#                 motifEnd <- motifStart + motifSize - 1L
#                 motifIRanges <- IRanges(start=motifStart, end=motifEnd)
#                 myRleToIRanges <- IRanges(
#                     start=rep(
#                         runValue(matchedReadsCDS),
#                         runLength(matchedReadsCDS)),
#                     width=1)
#                 myMotifCounts <- countOverlaps(motifIRanges, myRleToIRanges)
#                 myCodonCounts <- data.frame(cbind(1:length(myMotifCounts), myMotifCounts))
#
#             }
        }
        else{
            nbrCodons <- ceiling(length(listeRangesCDS[[ixTransc]]) / motifSize)
            myCodonCounts <- data.frame(cbind(1:nbrCodons, rep(0, nbrCodons)))
        }
        names(myCodonCounts) <- c("codonID","nbrReads");
        list(
            c(
                sum(S4Vectors::runLength(matchedReadsCDS)),
                sum(S4Vectors::runLength(matchedReads5UTR)),
                sum(S4Vectors::runLength(matchedReads3UTR))
            ),
            myCodonCounts
        )
    })

    names(shiftedTranscMatches) <- names(transcWithReads)
    countsFeatures <- do.call(rbind,lapply(shiftedTranscMatches, `[[`, 1))
    colnames(countsFeatures) <- c("CDS_counts", "fiveUTR_counts", "threeUTR_counts")
    rownames(countsFeatures) <- names(shiftedTranscMatches)

    #additional info
    #strandInfo=S4Vectors::runValue(strand(transc_withReads))
    chrInfo <- S4Vectors::runValue(seqnames(transcWithReads))
    startInfo <- min(start(transcWithReads))
    endInfo <- max(end(transcWithReads))
    cdsInfo <- do.call(rbind, cdsPosTranscWithReads)
    cdsLength <- cdsInfo[, 2] - cdsInfo[, 1] + 1
    cdsStart <- cdsInfo[, 1]
    cdsEnd <- cdsInfo[, 2]

    countsData <-
        cbind(
            as.character(rownames(countsFeatures)),
            as.character(unlist(chrInfo)),
            as.character(unlist(strandInfo)),
            startInfo,
            endInfo,
            newTranscWidth,
            cdsStart,
            cdsEnd,
            cdsLength,
            countsFeatures
        )
    colnames(countsData) <-
        c(
            "gene",
            "chr",
            "strand",
            "transc_genomic_start",
            "transc_genomic_end",
            "transc_length",
            "orf_start",
            "orf_end",
            "orf_length",
            colnames(countsFeatures)
        )

    codonReadCoverage <- lapply(shiftedTranscMatches, `[[`, 2)
    names(codonReadCoverage) <- names(shiftedTranscMatches)
    return(list(as.data.frame(countsData), codonReadCoverage))
}
