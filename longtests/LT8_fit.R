require(testthat)

sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    is.array(S4Vectors::metadata(object)[[fit]]) &
    (nrow(S4Vectors::metadata(object)[[fit]])==nrow(object)) 
}

#==============================================================================

context('fit: GSE161731')

    # Prepare minimal full-feature datset: subgroup, block, weights
        require(magrittr)
        basedir <- tools::R_user_dir('autonomics', 'cache')
        subdir  <- file.path(basedir, 'GSE161731')
        if (!dir.exists(subdir))  GEOquery::getGEOSuppFiles("GSE161731", baseDir = basedir)
        file  <- paste0(subdir,'/GSE161731_counts.csv.gz')
        sfile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
        object <- .read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id')
        object %<>% rm_singleton_samples('subject_id')
        object %<>% filter_samples(cohort == 'COVID-19', verbose = TRUE)
        object %<>% filter_samples(time_since_onset != 'middle')
        complete_subjects <- sdt(object)[, .SD[all(c('early', 'late') %in% time_since_onset)], by='subject_id']$subject_id
        object %<>% filter_samples(subject_id %in% complete_subjects, verbose = TRUE)
        # object %<>% pca()
        # biplot(object, color = 'subject_id', group = 'subject_id', label = 'sample_id', shape = 'time_since_onset')
        object %<>% filter_samples(!sample_id %in% c('DU18-02S0011620', 'DU18-02S0011625', 'DU18-02S0011619'), verbose = TRUE)
        object$subject_id %<>% make.names()
        object %<>% preprocess_rnaseq_counts(formula = ~ time_since_onset, block = 'subject_id')
        # object %<>% pca()
        # biplot(object, color=subject_id,group=subject_id, shape=time_since_onset)

    test_that(  "fit_limma( formula = ~ time_since_onset)", {             # 0 down, 37 up
        # ~ time_since_onset
            object %<>% fit_limma(formula =     ~ time_since_onset, block = 'subject_id')                            # 37 up
            object %<>% fit_limma(formula = ~ 0 + time_since_onset, block = 'subject_id', contrasts = 'late-early')  # 37 up
    })

    
                                                 # is different: 12 (not 10)
                                                 #               52 (not 54)
    # https://stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html

        

    # fails
    # lme: nlminb problem, convergence error code = 1 singular convergence (7)
    # lmer: keeps running
    test_that(  "fit_lme(r): formula = ~ 0 + time_since_onset | subject_id", {
        object %<>% fit_lme(formula = ~time_since_onset, block = 'subject_id')
        #object %<>% fit_lmer(formula = ~time_since_onset, block='subject_id')
        #expect_true(sumexp_contains_fit(object))
        #expect_true(summarize_fit(fdt(object), 'limma')$ndown==8)
        #expect_true(summarize_fit(fdt(object), 'limma')$nup  ==86)
    })
    
    
    test_that("fit_limma: diff ~ 1", { # 10 down, 54 up
        object %<>% subtract_baseline(block = 'subject_id', subgroupvar = 'time_since_onset')
        # object %<>% pca()
        # biplot(object, 
        #       color=subject_id, group=subject_id, shape=time_since_onset)
        # 'subgroup1'
        object %<>% fit_limma()
        ndown <- summarize_fit(fdt(object), 'limma', 'late')$downp
        nup   <- summarize_fit(fdt(object), 'limma', 'late')$upp
        # NULL subgroup
        object$subgroup <- NULL
        object %<>% fit_limma()
        expect_true(summarize_fit(fdt(object), 'limma', 'Intercept')$downp==ndown)
        expect_true(summarize_fit(fdt(object), 'limma', 'Intercept')$upp  ==nup)
        # ~ 1
        object %<>% fit_limma(formula=~1)
        expect_true(sumexp_contains_fit(object))
        expect_true(summarize_fit(fdt(object), 'limma')$ndown==ndown)
        expect_true(summarize_fit(fdt(object), 'limma')$nup  ==nup)
    })

    # a bit slow    
    #test_that(  "lm: formula = ~1", {
    #    object$subgroup <- NULL
    #    object %<>% fit_lm(formula=~1)
    #    expect_true(sumexp_contains_fit(object))
    #    expect_true(summarize_fit(fdt(object), 'limma')$ndown==0)
    #    expect_true(summarize_fit(fdt(object), 'limma')$nup  ==0)
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
    




