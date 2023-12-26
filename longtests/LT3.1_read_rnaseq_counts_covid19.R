require(testthat)
context('read_rnaseq_counts: GSE161731')
    basedir <- tools::R_user_dir('autonomics', 'cache')
    subdir  <- file.path(basedir, 'GSE161731')
    if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir)
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')

test_that('read_rnaseq_counts(file, cpm = FALSE)', {
    object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', formula = ~cohort, cpm = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('log2counts' %in% assayNames(object))
})

test_that('read_rnaseq_counts(file)', {
    object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', formula = ~cohort)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2cpm')
})

test_that('read_rnaseq_counts(file, pca = TRUE)', {
    object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', formula = ~cohort, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca' )))
})

# test_that("read_rnaseq_counts(file, fit = 'limma')", {
#     object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', formula = ~0+cohort, 
#                                  fit = 'limma', contrasts = c('COVID.19-Bacterial'))
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('limma' %in% names(S4Vectors::metadata(object)))
# })

