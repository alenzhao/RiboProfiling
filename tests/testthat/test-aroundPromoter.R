library(testthat)
library(RiboProfiling)

data(ctrlGAlignments)
aln <- ctrlGAlignments
alnGRanges <- readsToStartOrEnd(aln, what="start")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

oneBinRanges <-
    aroundPromoter(
        txdb,
        alnGRanges
    )

test_that("aroundPromoter output type",{
  expect_is(oneBinRanges,"GRanges")
}
)

test_that("aroundPromoter output length positive",{
  expect_more_than(length(oneBinRanges),0)
}
)
