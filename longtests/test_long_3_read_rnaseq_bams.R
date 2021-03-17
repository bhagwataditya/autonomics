if (requireNamespace('Rsubread', quietly = TRUE)){
    context('read_rnaseq_bams')
    file <- download_data('billing16.bam.zip')
    test_that('read_rnaseq_bams(file, cpm=FALSE)', {
        object <- read_rnaseq_bams(
                    file, paired=TRUE, genome='hg38', cpm=FALSE, plot=FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_identical(SummarizedExperiment::assayNames(object),'log2counts')
    })

    test_that('read_rnaseq_bams(file)', {
        object <- read_rnaseq_bams(
                    file, paired=TRUE, genome='hg38', plot=FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_identical(SummarizedExperiment::assayNames(object),'log2cpm')
    })

    test_that('read_rnaseq_bams(file, pca = TRUE)', {
        object <- read_rnaseq_bams(
                    file, paired=TRUE, genome='hg38', pca=TRUE, plot=FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_identical(SummarizedExperiment::assayNames(object),'log2cpm')
        expect_true('pca1' %in% svars(object))
        expect_true('pca1' %in% fvars(object))
        expect_true('pca'  %in% names(S4Vectors::metadata(object)))
    })

    test_that("read_rnaseq_counts(file, fit='limma')", {
        object <- read_rnaseq_bams(
                    file, paired=TRUE, genome='hg38', fit='limma', plot=FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('limma' %in% names(S4Vectors::metadata(object)))
    })
}
