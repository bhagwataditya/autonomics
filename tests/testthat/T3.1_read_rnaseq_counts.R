require(testthat)
context('read_rnaseq_counts: GSE161731')

file <- download_data("billing19.rnacounts.txt")

test_that('read_rnaseq_counts(file, cpm = FALSE)', {
    object <- read_rnaseq_counts(file, cpm = FALSE)
    expect_s4_class(object, 'SummarizedExperiment') 
    expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2counts')
})

test_that('read_rnaseq_counts(file)', {
    object <- read_rnaseq_counts(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2cpm')
})


# Moved to longtest_3_read_rnaseq_counts
# --------------------------------------
#
# test_that('read_rnaseq_counts(file, pca = TRUE)', {
#     object <- read_rnaseq_counts(file, pca = TRUE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true(  scorenames('pca', 'sample_id', 1) %in% svars(object))
#     expect_true(loadingnames('pca', 'sample_id', 1) %in% fvars(object))
#     expect_true(  methodname('pca', 'sample_id')    %in% names(S4Vectors::metadata(object)))
# })
# 
# test_that("read_rnaseq_counts(file, fit = 'limma')", {
#     object <- read_rnaseq_counts(file, fit = 'limma')
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
# })
