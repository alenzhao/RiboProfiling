#' Associates the read counts on codons with the codon type for each ORF.
#'
#' @importFrom Biostrings reverseComplement
#' @param listReadsCodon a list of data.frame objects.
#' It contains the number of reads per codon in a CDS.
#' @param genomeSeq a BSgenome object.
#' It contains the full genome sequences for the organism.
#' @param orfCoord a GRangesList.
#' The coordinates of the ORFs on the genome.
#' @param motifSize an integer. The number of nucleotides in each motif
#' on which to compute coverage and usage. Either 3, 6, or 9.
#' Default 3 nucleotides (codon).
#' Attention! For long motifs, the function can be quite slow!!
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
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' #get the names of the ORFs
#' #grouped by transcript
#' cds <- GenomicFeatures::cdsBy(txdb, use.names=TRUE)
#' orfCoord <- cds[names(cds) %in% names(listReadsCodon)]
#'
#' #get the genome, please check that the genome has the same seqlevels
#' genomeSeq <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' #if not rename it
#' #gSeq <- GenomeInfoDb::renameSeqlevels(genomeSeq,
#' #sub("chr", "", GenomeInfoDb::seqlevels(genomeSeq)))
#'
#' #codon frequency, coverage, and annotation
#' codonData <- codonInfo(listReadsCodon, genomeSeq, orfCoord)
#' @export
#' @import S4Vectors
#' @import GenomicFeatures
#' @importFrom reshape2 dcast
#' @importFrom plyr ldply
#' @importFrom data.table data.table
#' @import sqldf


codonInfo <-
    function(
        listReadsCodon,
        genomeSeq,
        orfCoord,
        motifSize){
    if(!is(genomeSeq, "BSgenome")){
        stop("Param genomeSeq should be of BSgenome class!\n")
    }
    if(!is(orfCoord, "GRangesList")){
        stop("The ORF coordinates (orfCoord) are not a GRangesList class!\n")
    }
    if(missing(motifSize) ||
       !is(motifSize, "numeric") ||
       motifSize %% 1 != 0 ||
       motifSize <= 0 ||
       !(motifSize %in% c(3, 6, 9))){
        warning("Param motifSize should be an integer! Accepted values 3, 6 or 9. Default value is 3.\n")
        motifSize <- 3
    }
    .id <- NULL
    codonID <- NULL
    codon <- NULL
    stepSize <- motifSize

    #extract the orf sequence
    cdsSeqs <- GenomicFeatures::extractTranscriptSeqs(genomeSeq, orfCoord)
    orfNames <- names(cdsSeqs)

    #if motifSize >3 then create counts for overlapping 1 codon shifted motifs
    if(motifSize > 3){
        stepSize = 3
    }

    #codonUsage <- Biostrings::trinucleotideFrequency(cdsSeqs, step=3)
    #I use with.labels=F for rapidity
    codonUsage <- Biostrings::oligonucleotideFrequency(
        cdsSeqs,
        width=motifSize,
        step=stepSize,
        with.labels=F
        )
    #I launch with labels only on the the first sequence, to get the patterns
    testCodonUsage <- Biostrings::oligonucleotideFrequency(
        cdsSeqs[[1]],
        width=motifSize,
        step=stepSize
    )
    colnames(codonUsage) <- names(testCodonUsage)
    rownames(codonUsage) <- orfNames
    #codonUsage <- codonUsage[,colSums(codonUsage) !=0]
    codonUsage <- data.frame(codonUsage)

    #codonUsage <- ldply(codonUsage)

    dataListReadsCodonID <- plyr::ldply(.data=listReadsCodon)
    # system.time(ldply(.data=listReadsCodon))
    # user  system elapsed
    # 16.217   0.001  16.202
    #### ????here group the info: if 3 codons, than paste the 3 codons together and sum their info
    #something of the type
    codonTypeID <- lapply(as.character(cdsSeqs), function(x){
        getCodons(x, sizeMotif=motifSize, stepSize)
        })
    # system.time(sapply(lapply(cdsSeqs,as.character),getCodons))
    # user  system elapsed
    # 139.369   0.016 139.255

    codonTypeID <- plyr::ldply(codonTypeID)
    #   system.time(ldply(codonTypeID))
    #   user  system elapsed
    #   515.160  21.022 535.669
    ### ???? also on codonTypeID paste every 3 codons in a cds
    ### ddply(codonTypeID, ".id", summarise, seq=unname(tapply(codon, (seq_along(codon)-1) %/% 3, paste,collapse="")))
    codonTypeID$codonID <- as.numeric(as.character(codonTypeID$codonID))

    #merging is much faster on data.table then with merge
    dtCodonTypeID <- data.table(codonTypeID, key=c(".id", "codonID"))
    dtReads <- data.table(dataListReadsCodonID, key=c(".id", "codonID"))
    codonTypeCov <- data.frame(
        dtCodonTypeID[dtReads, list(.id, codonID, codon, nbrReads=dtReads$nbrReads)]
        )

#     codonTypeCov <- suppressMessages(plyr::join(codonTypeID, dataListReadsCodonID))
#     #   system.time(merge(codonTypeID, dataListReadsCodonID))
#     #   user  system elapsed
#     #   73.968   2.368  76.265
#    codonTypeCoverage <- aggregate(nbrReads~ .id + codon, codonTypeCov, sum)
    colnames(codonTypeCov)[1] <- "gene"
    codonTypeCoverage <- sqldf::sqldf(c("create index ix on codonTypeCov(gene, codon)",
                         "select gene, codon, sum(nbrReads) from main.codonTypeCov group by gene, codon"))
    codonTypeCoverage <- na.omit(codonTypeCoverage[order(codonTypeCoverage$gene), ])
    colnames(codonTypeCoverage) <- c(
        ".id", "codon", "nbrReads"
    )

#     codonTypeCoverageFrame <-
#         reshape::cast(codonTypeCoverage, .id ~ codon, value="nbrReads")
    codonTypeCoverageFrame <-
        reshape2::dcast(codonTypeCoverage, .id ~ codon, value.var="nbrReads")
    codonCovMatrix <- codonTypeCoverageFrame[, 2:ncol(codonTypeCoverageFrame)]
    rownames(codonCovMatrix) <- codonTypeCoverageFrame$.id
    #codonCovMatrix[is.na(codonCovMatrix)] <- 0
    for (j in seq_len(ncol(codonCovMatrix)))
        data.table::set(codonCovMatrix, which(is.na(codonCovMatrix[[j]])),j,0)
    codonCovMatrix <- as.data.frame(lapply(codonCovMatrix,as.numeric))
    rownames(codonCovMatrix) <- codonTypeCoverageFrame$.id

    #order the genes according to codonUsage order
    codonCovMatrix <-
        codonCovMatrix[match(rownames(codonUsage), rownames(codonCovMatrix)), ]

    return(list(codonUsage, codonCovMatrix))
}

