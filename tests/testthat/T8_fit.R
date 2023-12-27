require(testthat)
sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    any(stri_detect_fixed(fvars(object), paste0('~', fit)))
}

#============================================================================
#                                                                           #
             context(" fit ")                                               #
#                                                                           #
#============================================================================


# MASSPEC
    test_that(  "fit: fukuda20", {
        file <- download_data('fukuda20.proteingroups.txt') 
        object <- read_maxquant_proteingroups(file)
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })

    test_that('fit: billing19.rnacounts', {
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file)
        expect_true(sumexp_contains_fit(fit_limma(object),                   'limma'))
        expect_true(sumexp_contains_fit(fit_limma(object, weightvar = NULL), 'limma'))
        # expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))         # slow
        # expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))   # slow
    })
    
# SOMASCAN
    test_that(  "fit: atkin.somascan", {
        
        # Original
            file <- download_data('atkin.somascan.adat')
            object <- read_somascan(file)
            expect_true(sumexp_contains_fit(
                    fit_limma(object, ~ Time + Diabetes, block = 'Subject'), 'limma'))
            
            expect_true(sumexp_contains_fit(
                       fit_lm(object, ~ Time + Diabetes                   ), 'lm'))
            
            expect_true(sumexp_contains_fit(
                      fit_lme(object, ~ Time + Diabetes, block = 'Subject'), 'lme'))
            
            expect_true(sumexp_contains_fit(
                     fit_lmer(object, ~ Time + Diabetes, block = 'Subject'), 'lmer'))
            
            expect_true(sumexp_contains_fit(
                 fit_wilcoxon(object, ~ Time,            block = 'Subject'), 'wilcoxon'))
            
        # Subtracted
            object %<>% subtract_differences(block = 'Subject', subgroupvar ='Time')
            expect_true(sumexp_contains_fit(
                       fit_lm(object, formula = ~ 0 + Time), 'lm'))
            expect_true(sumexp_contains_fit(
                    fit_limma(object, formula = ~ 0 + Time), 'limma'))
    })
    
# METABOLON
    test_that(  "fit: atkin.metabolon", {
        # read
            file <- download_data('atkin.metabolon.xlsx')
            object <- read_metabolon(file)
        # limma
            object %<>% fit_limma(formula = ~ Time + Diabetes, block = 'Subject')
            fitdt <- summarize_fit(fdt(object), fit = 'limma')
            expect_s3_class(fitdt, 'data.table')
        # lm
            object %<>% fit_lm(formula = ~ Time + Diabetes)
            fitdt <- summarize_fit(fdt(object), fit = 'lm')
            expect_s3_class(fitdt, 'data.table')
        # lme
            object %<>% fit_lme(formula = ~ Time + Diabetes, block = 'Subject')
            fitdt <- summarize_fit(fdt(object), fit = 'lme')
            expect_s3_class(fitdt, 'data.table')
        # lmer
            suppressWarnings(object %<>% fit_lmer(formula = ~ Time + Diabetes, block = 'Subject'))
            fitdt <- summarize_fit(fdt(object), fit = 'lmer')
            expect_s3_class(fitdt, 'data.table')
        # wilcoxon        
            object %<>% fit_wilcoxon(formula = ~ subgroup, block = 'Subject')
            fitdt <- summarize_fit(fdt(object), fit = 'wilcoxon')
            expect_s3_class(fitdt, 'data.table')
    })

        
#============================================================================
#                                                                           #
             context(" plot_contrastogram ")                                #
#                                                                           #
#============================================================================
    
if (requireNamespace('diagram', quietly = TRUE)){
    msg <- 'plot_contrastogram("fukuda20.proteingroups.txt")'
    file <-  download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file)
    expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
}

