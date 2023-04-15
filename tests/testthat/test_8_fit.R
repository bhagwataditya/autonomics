sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    any(stri_detect_fixed(fvars(object), paste0('~', fit)))
}

# MASSPEC
    test_that(  "fit: fukuda20", {
        file <- download_data('fukuda20.proteingroups.txt') 
        object <- read_maxquant_proteingroups(file, plot = FALSE)
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })

    test_that('fit: billing19.rnacounts', {
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file)
        #expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        #expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        object %<>% fit_limma()
        expect_true(sumexp_contains_fit(object,    'limma'))
        
        object %<>% fit_limma(weightvar = NULL)
        expect_true(sumexp_contains_fit(object,    'limma'))
        expect_s3_class(summarize_fit(fdt(object), fit = 'limma'), 'data.table')
    })
    
# SOMASCAN
    test_that(  "fit: atkin18.somascan", {
        file <- download_data('atkin18.somascan.adat')
        object <- read_somascan(file)
        
        object %<>% fit_limma( ~subgroup + T2D,  block = 'Subject_ID', coding = 'reference')
        expect_true(sumexp_contains_fit(object, 'limma'))
        
        object %<>% fit_lm(~ subgroup + T2D, coding = 'reference')
        expect_true(sumexp_contains_fit(object, 'lm'))
        
        object %<>% fit_lme(~ subgroup + T2D, block = 'Subject_ID', coding = 'reference')
        expect_true(sumexp_contains_fit(object, 'lme'))
        
        object %<>% fit_lmer(~ subgroup + T2D, block = 'Subject_ID', coding = 'reference')
        expect_true(sumexp_contains_fit(object, 'lmer'))
        
        object %<>% fit_wilcoxon(~ subgroup, block = 'Subject_ID', coding = 'reference')
        expect_true(sumexp_contains_fit(object, 'wilcoxon'))
        
        object %<>% subtract_differences(block = 'Subject_ID', subgroupvar ='subgroup')
        object %<>% fit_lm(formula = ~ 0 + subgroup)
        expect_true(sumexp_contains_fit(object, 'lm'))
        
        object %<>% fit_limma(formula = ~ 0 + subgroup)
        expect_true(sumexp_contains_fit(object, 'limma'))
    })
    
# METABOLON
    test_that(  "fit: atkin18.metabolon", {
        
        # limma
        object %<>% fit_limma(formula = ~ subgroup + T2D, block = 'SUB', coding = 'reference')
        fitdt <- summarize_fit(fdt(object), fit = 'limma')
        expect_s3_class(fitdt, 'data.table')
        
        # lm
        object %<>% fit_lm(formula = ~ subgroup + T2D, coding = 'reference')
        fitdt <- summarize_fit(fdt(object), fit = 'lm')
        expect_s3_class(fitdt, 'data.table')
        
        # lme
        object %<>% fit_lme(formula = ~ subgroup + T2D, block = 'SUB', coding = 'reference')
        fitdt <- summarize_fit(fdt(object), fit = 'lme')
        expect_s3_class(fitdt, 'data.table')

        # lmer
        object %<>% fit_lmer(formula = ~ subgroup + T2D, block = 'SUB', coding = 'reference')
        fitdt <- summarize_fit(fdt(object), fit = 'lmer')
        expect_s3_class(fitdt, 'data.table')

        # wilcoxon        
        object %<>% fit_wilcoxon(formula = ~ subgroup, block = 'SUB', coding = 'reference')
        fitdt <- summarize_fit(fdt(object), fit = 'wilcoxon')
        expect_s3_class(fitdt, 'data.table')
    })
    
    
if (requireNamespace('diagram', quietly = TRUE)){
    context('plot_contrastogram')
    msg <- 'plot_contrastogram("fukuda20.proteingroups.txt")'
    file <-  download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file)
    expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
}

