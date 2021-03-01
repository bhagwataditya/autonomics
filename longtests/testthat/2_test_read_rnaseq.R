context('read_rnaseq_counts: GSE161731')
    basedir <- '~/autonomicscache/datasets'
    subdir  <- '~/autonomicscache/datasets/GSE161731'
    if (!dir.exists(subdir)){
        GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir) }
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')

test_that('read_rnaseq_counts(file, cpm=FALSE)', {
    object <- read_rnaseq_counts(
        file, sfile = sfile, sfileby ='rna_id', subgroupvar = 'gender', 
        cpm = FALSE)
    expect_s4_class(object, 'SummarizedExperiment') 
})

test_that('read_rnaseq_counts(file, cpm=TRUE)', {
    object <- read_rnaseq_counts(
        file, sfile = sfile, sfileby ='rna_id', subgroupvar = 'gender')
    expect_s4_class(object, 'SummarizedExperiment')
})
    # pca
        msg <- 'read_rnaseq_counts(pca=TRUE)'
        object <- read_rnaseq_counts(file, sfile = sfile, sfileby='rna_id',
                                     subgroupvar='gender', pca=TRUE)
        test_that(msg, {
                  expect_s4_class(object, 'SummarizedExperiment')
                  expect_true('pca1' %in% svars(object))
                  expect_true('pca1' %in% fvars(object))
                  expect_true('pca'  %in% names(S4Vectors::metadata(object)))})
    # limma
        msg <- 'read_rnaseq_counts(limma=TRUE)'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby='rna_id', subgroupvar='gender', 
            limma = TRUE)
        test_that(msg, {
            expect_s4_class(object, 'SummarizedExperiment')
            expect_true('limma' %in% names(S4Vectors::metadata(object)))})
    # limma + block        
        msg <- 'read_rnaseq_counts(limma=TRUE, block=blockvar)'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby='rna_id', subgroupvar='gender', 
            limma=TRUE, block='subject_id')
        test_that(msg, {
            expect_s4_class(object, 'SummarizedExperiment')
            expect_true('limma' %in% names(S4Vectors::metadata(object)))})
        
# read_rnaseq_bams: billing16 
    # download
        context('read_rnaseq_bams: billing16')
        dir <- download_data("billing16.bam.zip")
    # read
        object <- read_rnaseq_bams(dir, paired=TRUE, genome="hg38")
        msg <- 'read_rnaseq_counts("GSE161731_counts.csv.gz")'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # external gtf
        # gtffile <- download_gtf("Homo sapiens")
        # object <- read_rnaseq_bams(dir, paired=TRUE, genome=gtffile)


