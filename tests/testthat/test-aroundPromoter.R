library(testthat)
library(RiboProfiling)

data(ctrlGAlignments)
aln <- ctrlGAlignments
alnGRanges <- readsToReadStart(aln)

oneBinRanges <-
    aroundPromoter(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
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
