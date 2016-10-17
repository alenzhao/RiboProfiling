#' Reads in GAlignments converted to either Read Start (5') or End (3') Positions
#'
#' @param aln A GAlignments object of the BAM mapping file.
#' @param what A character object. Either "start" (the default) or "end"
#' for read start or read end.
#  Ex. what <- "end"
#' @return A GRanges object containing either the read start or end
#' genomic positions.
#' @examples
#' #read the BAM file into a GAlignments object using
#' #GenomicAlignments::readGAlignments
#' #the GAlignments object should be similar to ctrlGAlignments object
#' data(ctrlGAlignments)
#' aln <- ctrlGAlignments
#' #transform the GAlignments object into a GRanges object (faster processing)
#' alnGRanges <- readsToStartOrEnd(aln, what = "end")
#' @export


readsToStartOrEnd <-
    function(
        aln,
        what){
    #check that the aln object is a GAlignments
    stopifnot(is(aln, "GAlignments"))

    if((what != "start" && what != "end") || missing(what)){
        warning("what parameter is invalid. Set to default start value!\n")
        what <- "start"
    }

    #find the indexes of the reads on the reverse strand
    ixMinusStrand <- which(BiocGenerics::strand(aln) == "-")

    if(what == "start"){
        borderRead <- IRanges::IRanges(start=BiocGenerics::start(aln), width=1)
        #if the read is on the reverse strand, the start position is the end
        borderRead[ixMinusStrand] <- IRanges::IRanges(start=BiocGenerics::end(aln[ixMinusStrand]), width=1)
    }
    if(what == "end"){
        borderRead <- IRanges::IRanges(end=BiocGenerics::end(aln), width=1)
        #if the read is on the reverse strand, the end position is the start
        borderRead[ixMinusStrand] <- IRanges::IRanges(end=BiocGenerics::start(aln[ixMinusStrand]), width=1)

    }

    alnGRanges <-
        GenomicRanges::GRanges(
            GenomeInfoDb::seqnames(aln),
            borderRead,
            strand=BiocGenerics::strand(aln),
            score=GenomicAlignments::cigarWidthAlongReferenceSpace(
                GenomicAlignments::cigar(aln), N.regions.removed=TRUE)
        )

    return(alnGRanges)
}