#' RiboProfiling.
#'
#' @name RiboProfiling
#' @docType package
#' @import BiocGenerics
NULL

#' Ribosome profiling data on chr1 in human primary BJ fibroblasts control data: PMID: 23594524.
#'
#' A dataset containing the alignment information on chromosome 1 from the control BAM.
#' The data object is a GAlignments object containing 3,504,859 hg19 mapped reads.
#'
#' @docType data
#' @keywords datasets
#' @name ctrlGAlignments
#' @usage data(ctrlGAlignments)
#' @format A GAlignments object with 3,504,859 reads.
#' @return the GAlignments object of reads on chr 1
NULL

#' Codon frequency and coverage in ORFs on chromosome 1, for dataset ctrlGAlignments
#'
#' A list of 2 data.frame objects:
#' one with the number of times each codon type is found in each ORF and
#' one with the number of reads for each codon type in each ORF.
#'
#' @docType data
#' @keywords datasets
#' @name codonDataCtrl
#' @usage data(codonDataCtrl)
#' @format A list of 2 lists.
#' @return A list of 2 lists.
NULL

#' The read coverage for each codon in ORFs on chromosome 1, for dataset ctrlGAlignments
#'
#' A list containing the number of reads for each codon in each ORF.
#' Codons are reported on their index in the ORF and no information is available about their type/sequence.
#'
#' @docType data
#' @keywords datasets
#' @name codonIndexCovCtrl
#' @usage data(codonIndexCovCtrl)
#' @format A list of 2 columns dataframes.
#' @return A list of 2 columns dataframes.
NULL



