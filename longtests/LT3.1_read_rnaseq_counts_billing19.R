require(testthat)

file <- download_data("billing19.rnacounts.txt")

test_that('read_rnaseq_counts(file, pca = TRUE)', {
    object <- read_rnaseq_counts(file, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(  scorenames('pca', 'sample_id', 1) %in% svars(object))
    expect_true(loadingnames('pca', 'sample_id', 1) %in% fvars(object))
    expect_true(  methodname('pca', 'sample_id')    %in% names(S4Vectors::metadata(object)))
})

test_that("read_rnaseq_counts(file, fit = 'limma')", {
    object <- read_rnaseq_counts(file, fit = 'limma')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})

