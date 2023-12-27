require(testthat)
if (requireNamespace('Rsubread', quietly = TRUE)){
    context('read_rnaseq_bams')
    file <- download_data('billing16.bam.zip')
    test_that('read_rnaseq_bams(file, cpm = FALSE)', {
        object <- read_rnaseq_bams(file, paired = TRUE, genome = 'hg38', cpm = FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2counts' %in% SummarizedExperiment::assayNames(object))
    })

    test_that('read_rnaseq_bams(file)', {
        object <- read_rnaseq_bams(file, paired = TRUE, genome = 'hg38')
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2cpm' %in% SummarizedExperiment::assayNames(object))
    })

    test_that('read_rnaseq_bams(file, pca = TRUE)', {
        object <- read_rnaseq_bams(file, paired = TRUE, genome = 'hg38', pca = TRUE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2cpm' %in% SummarizedExperiment::assayNames(object))
        expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca' )))
    })

    test_that("read_rnaseq_counts(file, fit = 'limma')", {
        object <- read_rnaseq_bams(file, paired = TRUE, genome = 'hg38', fit = 'limma')
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
    })
}
