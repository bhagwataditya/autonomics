context('read_rnaseq_counts: GSE161731')
    file <- download_data("billing19.rnacounts.txt")

test_that('read_rnaseq_counts(file, cpm=FALSE)', {
    object <- read_rnaseq_counts(file, cpm = FALSE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment') 
    expect_identical(SummarizedExperiment::assayNames(object), 'log2counts')
})

test_that('read_rnaseq_counts(file)', {
    object <- read_rnaseq_counts(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2cpm')
})

test_that('read_rnaseq_counts(file, pca=TRUE)', {
    object <- read_rnaseq_counts(file, pca = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('pca1' %in% svars(object))
    expect_true('pca1' %in% fvars(object))
    expect_true('pca'  %in% names(S4Vectors::metadata(object)))
})

test_that("read_rnaseq_counts(file, fit='limma')", {
    object <- read_rnaseq_counts(file, fit = 'limma', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(S4Vectors::metadata(object)))
})

