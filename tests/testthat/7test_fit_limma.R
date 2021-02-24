# Prepare minimal full-feature datset: subgroup, block, weights
    require(magrittr)
    basedir <- '~/autonomicscache/datasets'
    subdir  <- '~/autonomicscache/datasets/GSE161731'
    if (!dir.exists(subdir)){
        GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir) }
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
    object0 <- .read_rnaseq_counts(file, sfile=sfile, sfileby='rna_id')
    object0 %<>% rm_singleton_samples('subject_id')
    object0 %<>% filter_samples(cohort == 'COVID-19', verbose=TRUE)
    complete_subjects <- data.table(sdata(object0))[, 
        .SD[all(c('early', 'late') %in% time_since_onset)], 
        by='subject_id']$subject_id
    object0 %<>% filter_samples(subject_id %in% complete_subjects)
    object0 %<>% filter_samples(time_since_onset != 'middle')
    object0 %<>% preprocess_rnaseq_counts(subgroupvar = 'subject_id')
    object <- subtract_baseline(object0, 
                    block='subject_id', subgroupvar='time_since_onset')
    # object0 %<>% pca()
    # biplot(object0, color=subject_id,group=subject_id, shape=time_since_onset)
    # object %<>% pca()
    # biplot(object, color=subject_id, group=subject_id, shape=time_since_onset)

    
#=====================================
    
context('fit_limma: ~ 1')       # late - early : 32 down, 75 up
    
    sumexp_has_limma <- function(object){
        is(object, 'SummarizedExperiment') &
        is.array(limma(object)) &
        (nrow(limma(object))==nrow(object)) & 
        (ncol(limma(object))==length(unlist(contrastdefs(object)))) }
    
    test_that(  "subgroup = 'subgroup1'", {
        object %<>% fit_limma()
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==32)
        expect_true(summarize_fit(object, 'limma')$nup  ==75)
    })

    test_that(  "subgroup = NULL", {
        object$subgroup <- NULL
        object %<>% fit_limma()
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==32)
        expect_true(summarize_fit(object, 'limma')$nup  ==75)
    })
    
    test_that(  "formula=~1", {
        object$subgroup <- NULL
        object %<>% fit_limma(formula=~1)
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==32)
        expect_true(summarize_fit(object, 'limma')$nup  ==75)
    })

#=============================================================================
    
context('fit_limma: ~ (0 +) time_since_onset')  # 23 down, 90 up
    
    test_that(  "object$subgroup", {
        object <- object0
        object$subgroup <- object$time_since_onset
        object %<>% fit_limma()
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==23)
        expect_true(summarize_fit(object, 'limma')$nup  ==90)
    })

    test_that(  "subgroupvar = 'time_since_onset'", {
        object$subgroup <- NULL
        object %<>% fit_limma(subgroupvar = 'time_since_onset')
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==23)
        expect_true(summarize_fit(object, 'limma')$nup  ==90)
    })
    
    test_that(  "formula = ~ 0 + time_since_onset", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ 0 + time_since_onset)
        expect_true(sumexp_has_limma(object))
    })
    
    test_that(  "formula = ~ time_since_onset", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ time_since_onset)
        expect_true(sumexp_has_limma(object))
    })

#=============================================================================
    
context('fit_limma: ~ subject_id + subgroup')  # 90 down, 147 up
    
    test_that(  "formula = ~ subject_id + time_since_onset", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ subject_id + time_since_onset, 
                              contrastdefs = 'late')
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==90)
        expect_true(summarize_fit(object, 'limma')$nup  ==147)
    })
    
    test_that(  "formula = ~ time_since_onset + subject_id", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ time_since_onset + subject_id)
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==90)
        expect_true(summarize_fit(object, 'limma')$nup  ==147)
    })
    
    test_that(  "formula = ~ 0 + subject_id + time_since_onset", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ 0 + subject_id + time_since_onset, 
                              contrastdefs = 'late')
        expect_true(sumexp_has_limma(object))
        expect_true(summarize_fit(object, 'limma')$ndown==90)
        expect_true(summarize_fit(object, 'limma')$nup  ==147)
    })
    
    test_that(  "formula = ~ 0 + time_since_onset + subject_id", {
        object <- object0       # needs to be avoided !
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ 0 + time_since_onset + subject_id)
        expect_true(summarize_fit(object, 'limma')$ndown<90)
        expect_true(summarize_fit(object, 'limma')$nup  <147)
        expect_true(sumexp_has_limma(object))
    })

        
#=============================================================================
    
context('fit_limma: ~ subject_id  |  subgroup')  # 53 down, 127 up
    
    test_that(  "formula = ~ time_since_onset | subject_id", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~time_since_onset, block='subject_id')
        expect_true(summarize_fit(object, 'limma')$ndown== 53)
        expect_true(summarize_fit(object, 'limma')$nup  ==127)
        expect_true(sumexp_has_limma(object))
    })
    
    test_that(  "formula = ~ 0 + time_since_onset | subject_id", {
        object <- object0
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~0+time_since_onset, block='subject_id')
        expect_true(summarize_fit(object, 'limma')$ndown==53)
        expect_true(summarize_fit(object, 'limma')$nup  ==127)
        expect_true(sumexp_has_limma(object))
    })
    
#=============================================================================
    
context('fit_limma: ~ subject_id  + time_since_onset')  # 0 down, 1 up
# TODO: need to add weights
# Is colnames / row - rownames /col correct in .lmx?

    #test_that(  "formula = ~ time_since_onset + subject_id", {
    #    object <- object0
    #    object$subgroup <- NULL
    #    object %<>% fit_lm(formula = ~subject_id + time_since_onset)
    #    expect_true(summarize_fit(object, 'limma')$ndown== 0)
    #    expect_true(summarize_fit(object, 'limma')$nup  ==1)
    #    expect_true(sumexp_has_limma(object))
    #})
    
    #test_that(  "formula = ~ subject_id + time_since_onset", {
    #    object <- object0
    #    object$subgroup <- NULL
    #    object %<>% fit_lm(formula = ~time_since_onset+subject_id)
    #    expect_true(summarize_fit(object, 'limma')$ndown== 0)
    #    expect_true(summarize_fit(object, 'limma')$nup  ==1)
    #    expect_true(sumexp_has_limma(object))
    #})
    
    
    

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

