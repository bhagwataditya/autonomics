context('read_rnaseq_bams: billing16')
dir <- download_data("billing16.bam.zip")
    object <- read_rnaseq_bams(
        dir, paired=TRUE, genome="hg38", pca=FALSE, lmfit=FALSE, plot=FALSE)
    msg <- 'read_rnaseq_counts("GSE161731_counts.csv.gz") works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # external gtf
        # gtffile <- download_gtf("Homo sapiens")
        # object <- read_rnaseq_bams(dir, paired=TRUE, genome=gtffile)

    
context('read_rnaseq_counts: GSE161731')
basedir <- '~/autonomicscache/datasets'
subdir  <- '~/autonomicscache/datasets/GSE161731'
if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir)
file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
    # read
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby ='rna_id', subgroupvar = 'gender', 
            cpm = FALSE, voom = FALSE, pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_rnaseq_counts("GSE161731_counts.csv.gz") works'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # cpm
        msg <- 'read_rnaseq_counts(cpm=TRUE) works'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby ='rna_id', subgroupvar = 'gender', 
            voom = FALSE,  pca = FALSE, lmfit = FALSE, plot = FALSE)
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # voom
        msg <- 'read_rnaseq_counts(voom=TRUE) works'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby ='rna_id', subgroupvar = 'gender', 
            pca = FALSE, lmfit = FALSE, plot = FALSE)
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # block        
        msg <- 'read_rnaseq_counts(block=blockvar) works'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby='rna_id', subgroupvar='gender', 
            block='subject_id', pca=FALSE, lmfit=FALSE, plot=FALSE)
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # pca
        msg <- 'read_rnaseq_counts(pca=TRUE) works'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby='rna_id', subgroupvar='gender', 
            pca=TRUE, lmfit=FALSE, plot=FALSE)
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # lmfit
        msg <- 'read_rnaseq_counts(lmfit=TRUE) works'
        object <- read_rnaseq_counts(
            file, sfile = sfile, sfileby='rna_id', subgroupvar='gender',
            pca=FALSE, plot=FALSE)
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # plot
        msg <- 'read_rnaseq_counts(lmfit=TRUE) works'
        object <- read_rnaseq_counts(
                file, sfile = sfile, sfileby='rna_id', subgroupvar='gender')
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
