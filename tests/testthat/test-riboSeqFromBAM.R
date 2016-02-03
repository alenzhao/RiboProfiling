library(testthat)
library(RiboProfiling)

myBamFile <- system.file("extdata", "ctrl_sample.bam", package="RiboProfiling")
listeInputBam <- c(myBamFile)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
covData <-
    suppressMessages(
        suppressWarnings(
            riboSeqFromBAM(
                listeInputBam,
                txdb=txdb,
                listShiftValue=c(-14)
            )
        )
    )

test_that("riboSeqFromBam output type",{
    expect_is(covData,"list")
    expect_equal(length(covData),1)
    expect_is(covData[[1]],"list")
    expect_equal(length(covData[[1]]),2)
    expect_is(covData[[1]][[1]],"data.frame")
    expect_is(covData[[1]][[2]],"list")
}
)