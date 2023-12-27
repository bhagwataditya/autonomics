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


    test_that( " fit: billing19.rnacounts ", {
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file)
        expect_true(sumexp_contains_fit(fit_limma(object),                   'limma'))
        expect_true(sumexp_contains_fit(fit_limma(object, weightvar = NULL), 'limma'))
        # expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))         # slow
        # expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))   # slow
    })
    
    test_that( " fit: fukuda20.proteingroups ", {
        file <- download_data('fukuda20.proteingroups.txt') 
        object <- read_maxquant_proteingroups(file)
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })

    test_that( " fit: atkin.somascan ", {
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
                               fit_lm(object, ~ 0 + Time), 'lm'))
            expect_true(sumexp_contains_fit(
                            fit_limma(object, ~ 0 + Time), 'limma'))
    })
    
    test_that( " fit: atkin.metabolon ", {
        # read
            file <- download_data('atkin.metabolon.xlsx')
            object <- read_metabolon(file)
        # test
            expect_true(sumexp_contains_fit(
                              fit_limma(object, ~ Time + Diabetes, block = 'Subject'),  'limma'))
            expect_true(sumexp_contains_fit(
                                 fit_lm(object, ~ Time + Diabetes),                     'lm'))
            expect_true(sumexp_contains_fit(
                                fit_lme(object, ~ Time + Diabetes, block = 'Subject'),  'lme'))
            expect_true(sumexp_contains_fit(
              suppressWarnings(fit_lmer(object, ~ Time + Diabetes, block = 'Subject')), 'lmer'))
            expect_true(sumexp_contains_fit(
                           fit_wilcoxon(object, ~ subgroup,        block = 'Subject'),  'wilcoxon'))
    })
    
    test_that(  " fit: halama18.metabolon ", {
        file <- download_data('halama18.metabolon.xlsx')
        object <- read_metabolon(file)
        expect_true(sumexp_contains_fit(   fit_limma(object), 'limma'))
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(      fit_lm(object), 'lm'))
    })
    

        
#============================================================================
#                                                                           #
             context(" plot_contrastogram ")                                #
#                                                                           #
#============================================================================
    
if (requireNamespace('diagram', quietly = TRUE)){
    test_that(' plot_contrastogram: fukuda20.proteingroups', {
        # Read
            file <-  download_data('fukuda20.proteingroups.txt')
            object <- read_maxquant_proteingroups(file)
        # Test
            expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
})}

    
if (requireNamespace('diagram', quietly = TRUE)){
    test_that(' plot_contrastogram: billing19.proteingroups', {
        # Read
            file <-  download_data('billing19.proteingroups.txt')
            subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
            subgroups %<>% paste0('_STD')
            object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'limma')
        # Test
            expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup', curve = 0.8))
})}


if (requireNamespace('diagram', quietly = TRUE)){
    test_that(' plot_contrastogram: billing16.proteingroups', {  # ratios : self-contrasts
        # Read
            file <- download_data('billing16.proteingroups.txt')
            invert <- c('EM_E', 'BM_E', 'BM_EM')
            object <- read_maxquant_proteingroups(
                        file, invert = invert, fit = 'limma', formula = ~ 0 + subgroup)
        # Test
            expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
})}

