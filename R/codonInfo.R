#' Associates the read counts on codons with the codon type for each ORF.
#'
#' @param listReadsCodon a list of data.frame objects.
#' It contains the number of reads per codon in a CDS.
#' @param genomeSeq a BSgenome object.
#' It contains the full genome sequences for the organism.
#' @param orfCoord a GRangesList.
#' The coordinates of the ORFs on the genome.
#' @return a list of 2 data.frame objects:
#' one with the number of times each codon type is found in each ORF and
#' one with the number of reads for each codon type in each ORF.
#' @examples
#' #for each codon in each ORF get the read coverage
#' #parameter listReadsCodon can be returned by the riboSeqFromBam function
#' #it corresponts to the 2nd element in the list returned by riboSeqFromBam
#' data(codonIndexCovCtrl)
#' listReadsCodon <- codonIndexCovCtrl
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' #get the names of the ORFs
#' #grouped by transcript
#' cds <- cdsBy(txdb, use.names=TRUE)
#' orfCoord <- cds[names(cds) %in% names(listReadsCodon)]
#'
#' #get the genome, please check that the genome has the same seqlevels
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genomeSeq <- BSgenome.Hsapiens.UCSC.hg19
#' #if not rename it
#' #gSeq <- renameSeqlevels(genomeSeq,
#' #sub("chr", "", seqlevels(genomeSeq)))
#'
#' #codon frequency, coverage, and annotation
#' codonData <- codonInfo(listReadsCodon, genomeSeq, orfCoord)
#' @export
#' @import S4Vectors Biostrings GenomicFeatures
#' @importFrom reshape cast
#' @importFrom plyr ldply


codonInfo <-
    function(
        listReadsCodon,
        genomeSeq,
        orfCoord){
    if(!is(genomeSeq, "BSgenome")){
        stop("Param genomeSeq should be of BSgenome class!\n")
    }
    if(!is(orfCoord, "GRangesList")){
        stop("The ORF coordinates (orfCoord) are not a GRangesList class!\n")
    }

    #extract the orf sequence
    cdsSeqs <- extractTranscriptSeqs(genomeSeq, orfCoord)
    orfNames <- names(cdsSeqs)

    codonUsage <- trinucleotideFrequency(cdsSeqs, step=3)
    codonUsage <- data.frame(codonUsage)
    rownames(codonUsage) <- orfNames
    #codonUsage <- ldply(codonUsage)

    dataListReadsCodonID <- ldply(.data=listReadsCodon)
    # system.time(ldply(.data=listReadsCodon))
    # user  system elapsed
    # 16.217   0.001  16.202
    codonTypeID <- sapply(lapply(cdsSeqs, as.character), getCodons)
    # system.time(sapply(lapply(cdsSeqs,as.character),getCodons))
    # user  system elapsed
    # 139.369   0.016 139.255

    codonTypeID <- ldply(codonTypeID)
    #   system.time(ldply(codonTypeID))
    #   user  system elapsed
    #   515.160  21.022 535.669

    codonTypeCov <- merge(codonTypeID, dataListReadsCodonID)
    #   system.time(merge(codonTypeID, dataListReadsCodonID))
    #   user  system elapsed
    #   73.968   2.368  76.265
    codonTypeCoverage <- aggregate(nbrReads~ .id + codon, codonTypeCov, sum)
    codonTypeCoverage <- codonTypeCoverage[order(codonTypeCoverage$.id), ]

    codonTypeCoverageFrame<-
        cast(codonTypeCoverage, .id ~ codon, value="nbrReads")
    codonCovMatrix <- codonTypeCoverageFrame[, 2:ncol(codonTypeCoverageFrame)]
    rownames(codonCovMatrix) <- codonTypeCoverageFrame$.id
    codonCovMatrix[is.na(codonCovMatrix)] <- 0

    #order the genes according to codonUsage order
    codonCovMatrix<-
        codonCovMatrix[match(rownames(codonUsage), rownames(codonCovMatrix)), ]

    return(list(codonUsage, codonCovMatrix))
}

