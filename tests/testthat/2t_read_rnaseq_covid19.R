context('read_rnaseq_counts: GSE161731')
    basedir <- '~/autonomicscache/datasets'
    subdir  <- '~/autonomicscache/datasets/GSE161731'
    if (!dir.exists(subdir)){
        GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir) }
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')

test_that('read_rnaseq_counts(file, cpm=FALSE)', {
    object <- read_rnaseq_counts(file, sfile = sfile, sfileby = 'rna_id', 
                 subgroupvar = 'cohort', cpm = FALSE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment') 
    expect_identical(assayNames(object), 'log2counts')
})

test_that('read_rnaseq_counts(file)', {
    object <- read_rnaseq_counts(file, sfile = sfile, sfileby = 'rna_id', 
                 subgroupvar = 'cohort', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_identical(assayNames(object)[1], 'log2cpm')
})

test_that('read_rnaseq_counts(file, pca=TRUE)', {
    object <- read_rnaseq_counts(file, sfile = sfile, sfileby = 'rna_id',
                 subgroupvar='cohort', pca = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('pca1' %in% svars(object))
    expect_true('pca1' %in% fvars(object))
    expect_true('pca'  %in% names(S4Vectors::metadata(object)))
})

test_that("read_rnaseq_counts(file, fit='limma')", {
    object <- read_rnaseq_counts(file, sfile = sfile, sfileby = 'rna_id', 
                subgroupvar='cohort', fit = 'limma', 
                contrastdefs = c('COVID.19-Bacterial'), plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(S4Vectors::metadata(object)))
})

