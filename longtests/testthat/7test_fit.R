require(magrittr)
context('fit'
        )
# UNPAIRED

    msg <- 'fit: fukuda20.proteingroups)'
    file <- download_data('fukuda20.proteingroups.txt') 
    object <- read_proteingroups(file, plot = FALSE)
    object %<>% fit_wilcoxon(plot=FALSE)
    object %<>% fit_lm(      plot=FALSE)
    object %<>% fit_limma(   plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    
    msg <- 'fit: billing19.proteingroups'
    file <- download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(file, select_subgroups = select, plot = FALSE)
    exprs(object) %<>% na_to_zero()
    object %<>% fit_wilcoxon(plot=FALSE)
    object %<>% fit_lm(      plot=FALSE)
    object %<>% fit_limma(   plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    
    msg <- 'fit: halama18.metabolon'
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    object %<>% impute_systematic_nondetects(plot=FALSE)
    object %<>% fit_wilcoxon(plot=FALSE)
    object %<>% fit_lm(      plot=FALSE)
    object %<>% fit_limma(   plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    
    msg <- 'fit: billing19.rnacounts'
    file <- download_data('billing19.rnacounts.txt')
    object <- read_rnaseq_counts(file, voom=TRUE, plot=FALSE)
    object %<>% fit_wilcoxon(plot=FALSE)
    object %<>% fit_lm(      plot=FALSE)
    object %<>% fit_limma(   plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))

# PAIRED
    
    msg <- 'fit: atkin18.somascan'
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, plot=FALSE)
    object %<>% fit_wilcoxon(plot=FALSE)
    object %<>% fit_lm(      plot=FALSE)
    object %<>% fit_limma()
    
    object %<>% fit_wilcoxon(block = 'Subject_ID', plot=FALSE)
    object %<>% fit_limma(   block = 'Subject_ID', plot=FALSE)
    object %<>% fit_lme(     block = 'Subject_ID', plot=FALSE)
    
    object %<>% subtract_controls(block='Subject_ID', subgroup = 'subgroup', refgroup='t0')
    object %<>% fit_lm(formula = value ~ 0 + subgroup)
    object %<>% fit_limma(contrastdefs = c('t1', 't2', 't3'))
    
#' pca(object, plot=TRUE, color=SET)

    msg <- 'fit: GSE161731.rnacounts'
    basedir <- '~/autonomicscache/datasets'
    subdir  <- '~/autonomicscache/datasets/GSE161731'
    if (!dir.exists(subdir)){
        GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir) }
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
    object <- read_rnaseq_counts(
        file, sfile = sfile, sfileby='rna_id', subgroupvar='gender', 
        block='subject_id', voom=TRUE)
    object %<>% fit_lm(      plot = FALSE)
    object %<>% fit_limma(   plot = FALSE)
    object %<>% fit_wilcoxon(plot = FALSE)
    object %<>% fit_limma(block='subject_id', plot = FALSE)
    #object %<>% fit_lme(  block='subject_id', plot = FALSE) # convergence error
    #object %<>% fit_lmer( block='subject_id', plot = FALSE)  # keeps running
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))

# L2R
    
    object %<>% subtract_controls()
            
context('summarize_fit')

# RNASEQCOUNTS
    # ~ 0 + subgroup | weights
        msg <- 'summarize_fit("billing19.rnacounts.txt")'
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file, fit='limma', plot=FALSE)
        fitdt <- summarize_fit(object, fit='limma')
        test_that(msg, expect_s3_class(fitdt, 'data.table'))

    # ~ 0 + subgroups
        msg <- 'summarize_fit("billing19.rnacounts.txt" no voom)'
        weights(object) <- NULL
        object %<>% fit_limma(plot=FALSE)
        fitdt <- summarize_fit(object, fit = 'limma')
        test_that(msg, expect_s3_class(fitdt, 'data.table'))

# METABOLON
    # ~ 0 + subgroup
        msg <- 'summarize_fit("atkin18.metabolon.xlsx")'
        file <- download_data('atkin18.metabolon.xlsx')
        object <- read_metabolon(file, limma=TRUE, block='SUB', plot=FALSE)
        fitdt <- summarize_fit(object, fit = 'limma')
        test_that(msg, expect_s3_class(fitdt, 'data.table'))

   # ~ 0 + subgroup | block
        msg <- 'summarize_fit("atkin18.metabolon.xlsx, block")'
        object %<>% fit_limma(plot=FALSE, block = 'SUB')
        fitdt <- summarize_fit(object, fit = 'limma')
        test_that(msg, expect_s3_class(fitdt, 'data.table'))

   # ~ 0 + subgroup + t2d | block
        msg <- 'summarize_fit("atkin18.metabolon.xlsx, block, T2D")'
        object %<>% fit_limma(formula=~0+subgroup+T2D, block='SUB', plot=FALSE)
        fitdt <- summarize_fit(object, fit = 'limma')
        test_that(msg, expect_s3_class(fitdt, 'data.table'))

context('plot_contrastogram')
    # subgroup vector
    msg <- 'plot_contrastogram("billing19.proteingroups")'
    file <-  download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(
                file, select_subgroups = select, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object, curve=0.8), NA))
    
    # subgroup vector
    msg <- 'plot_contrastogram("fukuda20.proteingroups.txt")'
    file <-  download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object), NA))
    
    # Ratios: self-contrasts
    msg <- 'plot_contrastogram("billing16.proteingroups.txt")'
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_proteingroups(
               file, invert_subgroups=invert, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object), NA))

        
context('plot_volcano')
    # proteingroup group ratios
    msg  <- 'plot_volcano("billing16.proteingroups.txt")'
    file <- download_data("billing16.proteingroups.txt")
    inv <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_proteingroups(
        file, invert_subgroups=inv, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))
    
    # metabolon intensities: complex design
    msg  <- 'plot_volcano("halama18.metabolon.xlsx")'
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object, ntop=0), 'gg'))

    # proteingroup internalstandard ratios
    msg  <- 'plot_volcano("billing19.proteingroups.txt")'
    file <-  download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(
                 file, select_subgroups = select, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))

