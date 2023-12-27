require(testthat)

#============================================================================
#                                                                           #
             context(" read_rnaseq_counts: billing19 ")                     #
#                                                                           #
#============================================================================


test_that(" read_rnaseq_counts: billing19, cpm = FALSE ", {
    # Read
        file <- download_data("billing19.rnacounts.txt")
        object <- read_rnaseq_counts(file, cpm = FALSE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment') 
        expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2counts')
})


test_that(" read_rnaseq_counts: billing19 ", {
    # Read
        file <- download_data("billing19.rnacounts.txt")
        object <- read_rnaseq_counts(file)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2cpm')
})


test_that(" read_rnaseq_counts: billing19, pca = TRUE ", {
    # Read
        file <- download_data("billing19.rnacounts.txt")
        object <- read_rnaseq_counts(file, pca = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(  scorenames('pca', 'sample_id', 1) %in% svars(object))
        expect_true(loadingnames('pca', 'sample_id', 1) %in% fvars(object))
        expect_true(  methodname('pca', 'sample_id')    %in% names(S4Vectors::metadata(object)))
})


test_that(" read_rnaseq_counts: billing19, fit = 'limma' ", {
    # Read
        file <- download_data("billing19.rnacounts.txt")
        object <- read_rnaseq_counts(file, fit = 'limma')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})


#============================================================================
#                                                                           #
             context(" read_rnaseq_counts: covid19 ")                       #
#                                                                           #
#============================================================================


test_that(" read_rnaseq_counts: covid19 ", {
    # Read
        basedir <- tools::R_user_dir('autonomics', 'cache')
        subdir  <- file.path(basedir, 'GSE161731')
        if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir = basedir)
        file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
        sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
        object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', formula = ~cohort)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_identical(SummarizedExperiment::assayNames(object)[1], 'log2cpm')
})


test_that(" read_rnaseq_counts: covid19, cpm = FALSE ", {
    # Read
        basedir <- tools::R_user_dir('autonomics', 'cache')
        subdir  <- file.path(basedir, 'GSE161731')
        if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir = basedir)
        file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
        sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
        object <- read_rnaseq_counts(
                        file, sfile = sfile, by.y = 'rna_id', formula = ~cohort, cpm = FALSE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2counts' %in% assayNames(object))
})


test_that(" read_rnaseq_counts: file, pca = TRUE ", {
    # Read
        basedir <- tools::R_user_dir('autonomics', 'cache')
        subdir  <- file.path(basedir, 'GSE161731')
        if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir = basedir)
        file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
        sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
        object <- read_rnaseq_counts(
            file, sfile = sfile, by.y = 'rna_id', formula = ~cohort, pca = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca' )))
})

test_that(" read_rnaseq_counts: file, fit = 'limma' ", {
    # Read
        basedir <- tools::R_user_dir('autonomics', 'cache')
        subdir  <- file.path(basedir, 'GSE161731')
        if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir = basedir)
        file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
        sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
        object <- read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id', 
                         formula = ~cohort, fit = 'limma')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})

