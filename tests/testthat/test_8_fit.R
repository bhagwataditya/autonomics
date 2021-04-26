sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    any(stri_detect_fixed(names(fdata(object)), paste0('~', fit)))
}

# UNPAIRED: wilcoxon generally fails
    
    test_that(  "fit: fukuda20", {
        file <- download_data('fukuda20.proteingroups.txt') 
        object <- read_proteingroups(file, plot = FALSE)
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })

    test_that('fit: billing19.rnacounts', {
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file, voom=TRUE, plot=FALSE)
        #expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        #expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        object %<>% fit_limma()
        expect_true(sumexp_contains_fit(object,    'limma'))
        
        object %<>% fit_limma(weightvar = NULL)
        expect_true(sumexp_contains_fit(object,    'limma'))
        expect_s3_class(summarize_fit(object, fit='limma'), 'data.table')
    })
    
# PAIRED: wilcoxon does work
    test_that(  "fit: atkin18.somascan", {
        file <- download_data('atkin18.somascan.adat')
        object <- read_somascan(file, plot=FALSE)
        
        object %<>% fit_wilcoxon(subgroupvar='SampleGroup')
        expect_true(sumexp_contains_fit(object, 'wilcoxon'))
        
        object %<>% fit_lm(      subgroupvar='SampleGroup')
        expect_true(sumexp_contains_fit(object, 'lm'))
        
        object %<>% fit_limma(   subgroupvar='SampleGroup')
        expect_true(sumexp_contains_fit(object, 'limma'))
        
        object %<>% fit_wilcoxon(subgroupvar='SampleGroup', block='Subject_ID')
        expect_true(sumexp_contains_fit(object, 'wilcoxon'))
        
        object %<>% fit_limma(   subgroupvar='SampleGroup', block='Subject_ID')
        expect_true(sumexp_contains_fit(object, 'limma'))
        
        #object %<>% fit_lme(     subgroupvar='SampleGroup', block='Subject_ID')
        #expect_true(sumexp_contains_fit(object, 'lme'))
        
        #object %<>% fit_lmer(    subgroupvar='SampleGroup', block='Subject_ID')
        #expect_true(sumexp_contains_fit(object, 'lmer'))
        
        object %<>% subtract_differences(
                        block='Subject_ID', subgroupvar ='SampleGroup')
        object %<>% fit_lm(   formula = ~ 0 + SampleGroup)
        expect_true(sumexp_contains_fit(object, 'lm'))
        
        object %<>% fit_limma(formula = ~ 0 + SampleGroup,
                            contrastdefs = c('t1_t0', 't2_t1', 't3_t2'))
        expect_true(sumexp_contains_fit(object, 'limma'))
    })
    
# METABOLON
    test_that(  "fit: atkin18.metabolon", {
        # ~ 0 + Group
        file <- download_data('atkin18.metabolon.xlsx')
        object <- read_metabolon(file, plot=FALSE)
        suppressWarnings(object %<>% fit_limma(subgroupvar = 'Group'))
        expect_true(sumexp_contains_fit(object, 'limma'))
        fitdt <- summarize_fit(object, fit = 'limma')
        expect_s3_class(fitdt, 'data.table')

        # ~ 0 + Group | block
        suppressWarnings(object %<>% fit_limma(
                                        subgroupvar = 'Group', block = 'SUB'))
        expect_true(sumexp_contains_fit(object, 'limma'))
        fitdt <- summarize_fit(object, fit = 'limma')
        expect_s3_class(fitdt, 'data.table')

        # ~ 0 + Group + t2d | block
        object %<>% fit_limma(formula=~0+Group+T2D, block='SUB', plot=FALSE)
        fitdt <- summarize_fit(object, fit = 'limma')
        expect_s3_class(fitdt, 'data.table')
    })
    
if (requireNamespace('diagram', quietly = TRUE)){
    context('plot_contrastogram')
    msg <- 'plot_contrastogram("fukuda20.proteingroups.txt")'
    file <-  download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, fit='limma', plot=FALSE)
    test_that(msg, expect_error(
        plot_contrastogram(object, subgroupvar = 'subgroup'), NA))
}

