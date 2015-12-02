library(testthat)
library(RiboProfiling)

data(ctrlGAlignments)
aln <- ctrlGAlignments
alnGRanges <- readsToReadStart(aln)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

oneBinRanges <-
    aroundPromoter(
        txdb,
        alnGRanges
    )

listPromoterCov <- readStartCov(
    alnGRanges,
    oneBinRanges,
    matchSize="all",
    fixedInterval=c(-20, 20),
    renameChr="aroundTSS",
    charPerc="perc"
)

test_that("readStartCov output type",{
    expect_is(listPromoterCov,"list")
    expect_is(listPromoterCov[[1]],"GRanges")
}
)

test_that("readStartCov output dim",{
    expect_equal(length(listPromoterCov[[1]]),41)
}
)
