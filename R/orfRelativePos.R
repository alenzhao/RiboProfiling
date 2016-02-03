#' Relative position of the start and stop codon along the transcript
#'
#' @param cdsTransc a GRangesList.
#' It contains the CDS coordinates grouped by transcript.
#' @param exonGRanges a GRangesList.
#' It contains the exon coordinates grouped by transcript.
#' @return a list.
#' A list of relative positions of the start and end of ORFs.
#' @examples
#' #make a txdb object containing the annotations for the specified species.
#' #In this case hg19.
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' #get all CDSs by transcript
#' cds <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
#'
#' #get all exons by transcript
#' exonGRanges <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
#'
#' #retrieve the positions of start and end codons relative to the transcript
#' cdsPosTransc <- orfRelativePos(cds, exonGRanges)
#' @export

orfRelativePos <-
    function(
        cdsTransc,
        exonGRanges){
    if(!is(cdsTransc, "GRangesList") || !is(exonGRanges, "GRangesList")){
        stop("Parameters cdsTransc and exonGRanges should be GRangesList.\n")
    }
    #keep only transcripts with CDS
    transcWithCDS <- exonGRanges[names(cdsTransc)]
    #transform the GRanges of the transcripts in 1bp bins
    binTransc <- applyShiftFeature(transcWithCDS, 0)

    strandCDS <- as.character(S4Vectors::runValue(strand(cdsTransc)))

    #now two analyses on the plus and minus strand separately
    ## plus strand
    ixPlus <- which(strandCDS != "-")
    cdsTranscPlus <- cdsTransc[ixPlus]
    shiftedTranscCdsPlus <- binTransc[ixPlus]

    #find the start and end codon genomic positions
    startCodonPlus <- min(start(cdsTranscPlus))
    endCodonPlus <- max(end(cdsTranscPlus))

    cdsPosTranscPlus <- lapply(seq_len(length(shiftedTranscCdsPlus)), function(ixCDS){
        ixStart <- which(shiftedTranscCdsPlus[[ixCDS]] == startCodonPlus[[ixCDS]]);
        ixEnd <- which(shiftedTranscCdsPlus[[ixCDS]] == endCodonPlus[[ixCDS]]);
        c(ixStart, ixEnd)})
    names(cdsPosTranscPlus) <- names(shiftedTranscCdsPlus)

    ## minus strand
    ixMinus <- which(strandCDS == "-")
    cdsTranscMinus <- cdsTransc[ixMinus]
    shiftedTranscCdsMinus <- binTransc[ixMinus]

    #find the start and end codon genomic positions
    startCodonMinus <- max(end(cdsTranscMinus))
    endCodonMinus <- min(start(cdsTranscMinus))

    cdsPosTranscMinus<-
        lapply(
            seq_len(length(shiftedTranscCdsMinus)),
            function(ixCDS){
                transc1bins <- rev(sort(shiftedTranscCdsMinus[[ixCDS]]));
                ixStart <- which(transc1bins == startCodonMinus[[ixCDS]]);
                ixEnd <- which(transc1bins == endCodonMinus[[ixCDS]]);
                c(ixStart, ixEnd)
            }
        )
    names(cdsPosTranscMinus) <- names(shiftedTranscCdsMinus)

    cdsPosTransc <- append(cdsPosTranscPlus, cdsPosTranscMinus)

    return(cdsPosTransc)
}
