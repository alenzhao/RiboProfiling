#' Transforms the reads in GAlignments to Read Start Positions in GRanges
#'
#' @param aln A GAlignments object of the BAM mapping file.
#' @return A GRanges object containing the read start genomic positions.
#' For the reverse strand the read position on the genome is the end of the read
#'
#' @examples
#' #read the BAM file into a GAlignments object using readGAlignments
#' #the GAlignments object should be similar to ctrlGAlignments object
#' data(ctrlGAlignments)
#' aln <- ctrlGAlignments
#' #transform the GAlignments object into a GRanges object (faster processing)
#' alnGRanges <- readsToReadStart(aln)
#' @export
#' @import methods BiocGenerics S4Vectors IRanges GenomicAlignments


readsToReadStart <-
    function(
        aln){
    #check that the aln object is a GAlignments
    stopifnot(is(aln, "GAlignments"))

    #find the indexes of the reads on the reverse strand
    ixMinusStrand <- which(strand(aln) == "-")

    startRead <- IRanges(start=start(aln), width=1)
    #if the read is on the reverse strand, the start position is the end
    startRead[ixMinusStrand] <- IRanges(start=end(aln[ixMinusStrand]), width=1)

    alnGRanges <-
        GRanges(
            seqnames(aln),
            startRead,
            strand=strand(aln),
            score=cigarWidthAlongReferenceSpace(
                cigar(aln), N.regions.removed=TRUE)
        )

    return(alnGRanges)
}
