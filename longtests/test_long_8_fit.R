sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    is.array(S4Vectors::metadata(object)[[fit]]) &
    (nrow(S4Vectors::metadata(object)[[fit]])==nrow(object)) 
}

#==============================================================================

context('fit: GSE161731')

    # Prepare minimal full-feature datset: subgroup, block, weights
    require(magrittr)
    basedir <- rappdirs::user_cache_dir(appname = 'autonomics')
    subdir  <- file.path(basedir, 'GSE161731')
    if (!dir.exists(subdir)){
        GEOquery::getGEOSuppFiles("GSE161731",baseDir=basedir) }
    file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
    sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
    object <- .read_rnaseq_counts(file, sfile=sfile, sfileby='rna_id')
    object %<>% rm_singleton_samples('subject_id')
    object %<>% filter_samples(cohort == 'COVID-19', verbose=TRUE)
    object %<>% filter_samples(time_since_onset != 'middle')
    complete_subjects <- data.table::data.table(sdata(object))[, 
        .SD[all(c('early', 'late') %in% time_since_onset)], 
        by='subject_id']$subject_id
    object %<>% filter_samples(subject_id %in% complete_subjects, verbose=TRUE)
    # object %<>% pca()
    # biplot(object, color = subject_id, group = subject_id, label = sample_id, 
    #       shape = time_since_onset)
    object %<>% filter_samples(!sample_id %in% c(
        'DU18-02S0011620', 'DU18-02S0011625', 'DU18-02S0011619'), verbose=TRUE)
    object %<>% preprocess_rnaseq_counts(formula=~subject_id+time_since_onset)
    # object %<>% pca()
    # biplot(object, color=subject_id,group=subject_id, shape=time_since_onset)

    test_that(  "fit_limma: ~ time_since_onset", {             # 0 own, 35 up
        # defaults
        object$subgroup <- object$time_since_onset
        object %<>% fit_limma()
        expect_true(sumexp_contains_fit(object))
        ndown <- summarize_fit(object, 'limma')$ndown
        nup   <- summarize_fit(object, 'limma')$nup
        # subgroupvar
        object %<>% fit_limma(subgroupvar = 'time_since_onset')
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
        # formula without intercept
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ 0 + time_since_onset)
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
        # fomula with intercept
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ time_since_onset)
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
    })

    
    test_that("fit_limma: ~ subject_id + time_since_onset", { # 10 down, 54 up
        # ~ subject_id + time_since_onset
        object %<>% fit_limma(formula = ~ subject_id + time_since_onset, 
                              contrasts = 'late')
        expect_true(sumexp_contains_fit(object))
        ndown <- summarize_fit(object, 'limma')$ndown==10
        nup   <- summarize_fit(object, 'limma')$nup  ==54
        # ~ time_since_onset + subject_id
        object %<>% fit_limma(formula = ~ time_since_onset + subject_id)
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
        # ~ 0 + subject_id + time_since_onset
        object %<>% fit_limma(formula = ~ 0 + subject_id + time_since_onset, 
                              contrasts = 'late')
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
        # 0 + time_since_onset + subject_id
        object$subgroup <- NULL
        object %<>% fit_limma(formula = ~ 0 + time_since_onset + subject_id)
        expect_true(sumexp_contains_fit(object)) # is different: 12 (not 10)
                                                 #               52 (not 54)
    }) # https://stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html

        
    # a bit slow    
    #test_that(  "fit_limma: ~ time_since_onset | subject_id", {
    #    # ~ time_since_onset | subject_id
    #    object %<>% fit_limma(formula=~0+time_since_onset, block='subject_id')
    #    expect_true(sumexp_contains_fit(object))
    #    expect_true(summarize_fit(object, 'limma')$ndown==8)
    #    expect_true(summarize_fit(object, 'limma')$nup  ==86)
    #    
    #    # ~ 0 + time_since_onset | subject_id
    #    object %<>% fit_limma(formula = ~time_since_onset, block='subject_id')
    #    expect_true(summarize_fit(object, 'limma')$ndown== 8)
    #    expect_true(summarize_fit(object, 'limma')$nup  == 85)
    #    expect_true(sumexp_contains_fit(object))
    #})

    # fails
    # lme: nlminb problem, convergence error code = 1 singular convergence (7)
    # lmer: keeps running
    # test_that(  "fit_lme(r): formula = ~ 0 + time_since_onset | subject_id", {
        #object$subgroup <- NULL
        #object %<>% fit_lme(formula = ~time_since_onset, block='subject_id')
        #object %<>% fit_lmer(formula = ~time_since_onset, block='subject_id')
        #expect_true(sumexp_contains_fit(object))
        #expect_true(summarize_fit(object, 'limma')$ndown==8)
        #expect_true(summarize_fit(object, 'limma')$nup  ==86)
    # })
    
    
    test_that("fit_limma: diff ~ 1", { # 10 down, 54 up
        object %<>% subtract_baseline(
                        block='subject_id', subgroupvar='time_since_onset')
        # object %<>% pca()
        # biplot(object, 
        #       color=subject_id, group=subject_id, shape=time_since_onset)
        # 'subgroup1'
        object %<>% fit_limma()
        expect_true(sumexp_contains_fit(object))
        ndown <- summarize_fit(object, 'limma')$ndown
        nup <- summarize_fit(object, 'limma')$nup
        # NULL subgroup
        object$subgroup <- NULL
        object %<>% fit_limma()
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
        # ~ 1
        object %<>% fit_limma(formula=~1)
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(object, 'limma')$ndown==ndown)
        expect_true(summarize_fit(object, 'limma')$nup  ==nup)
    })

    # a bit slow    
    #test_that(  "lm: formula = ~1", {
    #    object$subgroup <- NULL
    #    object %<>% fit_lm(formula=~1)
    #    expect_true(sumexp_contains_fit(object))
    #    expect_true(summarize_fit(object, 'limma')$ndown==0)
    #    expect_true(summarize_fit(object, 'limma')$nup  ==0)
    #})


#==============================================================================

# UNPAIRED: wilcoxon generally fails
    
    test_that(  "fit: billing19.proteingroups", {
        file <- download_data('billing19.proteingroups.txt')
        select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
        select %<>% paste0('_STD')
        object <- read_maxquant_proteingroups(
                    file, select_subgroups = select, plot = FALSE)
        values(object) %<>% na_to_zero()
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })
    
    test_that(  "fit: halama18.metabolon", {
        file <- download_data('halama18.metabolon.xlsx')
        object <- read_metabolon(file, plot = FALSE)
        object %<>% impute_consistent_nas(subgroup = Group, plot=FALSE)
        #object %<>% fit_wilcoxon(subgroupvar = 'Group')
        #expect_true(sumexp_contains_fit(object, 'wilcoxon'))
        #object %<>% fit_lm(subgroupvar = 'Group')
        #expect_true(sumexp_contains_fit(object, 'lm'))
        object %<>% fit_limma(object, subgroupvar = "Group")
        expect_true(sumexp_contains_fit(object, 'limma'))
    })

context('plot_contrastogram')
    # subgroup vector
    msg <- 'plot_contrastogram("billing19.proteingroups")'
    file <-  download_data('billing19.proteingroups.txt')
    subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    subgroups %<>% paste0('_STD')
    object <- read_maxquant_proteingroups(
                file, subgroups = subgroups, fit = 'limma', plot = FALSE)
    test_that(msg, expect_error(
        plot_contrastogram(object, subgroupvar = 'subgroup', curve = 0.8), NA))
    
    # Ratios: self-contrasts
    msg <- 'plot_contrastogram("billing16.proteingroups.txt")'
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_maxquant_proteingroups(
               file, invert = invert, fit = 'limma', plot = FALSE)
    test_that(msg, expect_error(
        plot_contrastogram(object, subgroupvar = 'subgroup'), NA))

        
context('plot_volcano')
    # proteingroup group ratios
    msg  <- 'plot_volcano("billing16.proteingroups.txt")'
    file <- download_data("billing16.proteingroups.txt")
    invert <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_maxquant_proteingroups(
        file, invert = invert, fit = 'limma', plot = FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))
    
    # metabolon intensities: complex design
    msg  <- 'plot_volcano("halama18.metabolon.xlsx")'
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, fit = 'limma', plot = FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object, ntop = 0), 'gg'))

    # proteingroup internalstandard ratios
    msg  <- 'plot_volcano("billing19.proteingroups.txt")'
    file <-  download_data('billing19.proteingroups.txt')
    subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    subgroups %<>% paste0('_STD')
    object <- read_maxquant_proteingroups(
                 file, subgroups = subgroups, fit = 'limma', plot = FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))

