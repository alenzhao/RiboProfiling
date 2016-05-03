#' Returns the flank size around the TSS for the x \% CDSs
#'
#' @param txdb a TxDb object containing the annotations info to intersect with
#'   the alignment files.
#' @param alnGRanges A GRanges object containing the alignment information. In
#'   order to improve the performance of this function the GAlignments BAM
#'   object should be transformed into a GRanges object containing the cigar
#'   match size information as metadata.
#' @param percBestExpressed a numeric [between 0 and 1]. The percentage of the
#'   best expressed CDSs on which to plot the coverage around the TSS. Necessary
#'   if the shiftValue parameter must be estimated. Default value 0.03 (3\%).
#' @param flankSize a numeric positive integer. How many bp left and right of
#'   the TSS should the coverage be performed? Ex. flankSize <- 20
#' @return A GRanges object containing the 1 bp ranges for the selected CDSs in
#'   the TSS defined flanking region.
#' @examples
#' #read the BAM into a GAlignments object using
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
#' @export
#' @import IRanges


aroundPromoter <-
    function(
        txdb,
        alnGRanges,
        percBestExpressed,
        flankSize){
    #check the parameter validity
    if(!is(txdb, "TxDb")){
        stop("The annotation parameter (txdb) is not a valid TxDb object!\n")
    }
    if(!is(alnGRanges, "GRanges")){
        stop("The alignment object (alnGRanges) is not a GRanges object!\n")
    }

    ##### GET THE BEST x\% covered CDSs
    #first use cdsByGene to get the counts on all exons for each gene

    #first get the chromosomes present in the alignment file
    seqInAlignment <- unique(seqnames(alnGRanges))
    if(length(seqInAlignment) <= 0){
        stop("There are no chromosomes with reads on them.
         No seqnames in the alnGRanges object.\n")
    }
    #it groups all described CDS exons by gene
    cdsAll <- GenomicFeatures::cdsBy(txdb, by="gene")

    #keep only CDS on the chromosomes for which there are reads
    #cdsAllChr=cdsAll[seqnames(cdsAll) %in% seqInAlignment]
    seqlevels(cdsAll, force=TRUE) <- as.character(seqInAlignment)
    cdsByGene <- cdsAll[S4Vectors::elementNROWS(cdsAll) != 0]

    #choose the best expressed CDSs
    #first get the CDS coverage
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
        warning(paste("The percBestExpressed parameter is not numeric or <0.",
                      " It has been set to the default ",
                      percBestExpressed, "!\n", sep=""))
        percBestExpressed <- 0.03
    }

    vecCountsPerGene <- SummarizedExperiment::assays(countsPCGenesAllExons)$counts
    quantCounts <- quantile(vecCountsPerGene[which(vecCountsPerGene > 0)], 1-percBestExpressed)
    if(quantCounts <= 0){
        stop("No gene had counts overlapping the CDS!\n")
    }

    indexGenesBestExpressed <- which(vecCountsPerGene >= quantCounts)
    geneNamesBestExpressed <-
        rownames(countsPCGenesAllExons)[indexGenesBestExpressed]
    cdsByBestExpressed <-
        cdsByGene[names(cdsByGene) %in% geneNamesBestExpressed]

    #Get the coverage n bases left and right from the start codon

    #choose the longest transcript per gene (very rapid)
    maxWidthTranscBestExpGene <- max(width(cdsByBestExpressed))
    cdsByBestExprLongTransc <-
        cdsByBestExpressed[width(cdsByBestExpressed) == maxWidthTranscBestExpGene, ]

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

    #transform the gene ranges into one GRanges object with info on the cds_id.
    oneBinRanges <-
        unlist(S4Vectors::endoapply(flankPromoterUniq, function(ixTSS){
        tmpOneBinGRanges <- unlist(tile(ixTSS, n=width(ixTSS)));
        values(tmpOneBinGRanges)=
            S4Vectors::DataFrame(
                values=0,
                idSeq=rep(ixTSS$cds_id[1], length(tmpOneBinGRanges))
            );
        tmpOneBinGRanges}))

    return(oneBinRanges)
}
