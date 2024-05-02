sumexp_contains_fit <- function(object, fit = 'limma'){
    is(object, 'SummarizedExperiment') &
    any(stri_detect_fixed(fvars(object), paste0('~', fit)))
}

#============================================================================
#                                                                           #
#            context(" fit ")                                               #
#                                                                           #
#============================================================================


    test_that( " fit: billing19.rnacounts ", {
        file <- system.file('extdata/billing19.rnacounts.txt', package = 'autonomics')
        object <- read_rnaseq_counts(file)
        expect_true(sumexp_contains_fit(fit_limma(object),                   'limma'))
        expect_true(sumexp_contains_fit(fit_limma(object, weightvar = NULL), 'limma'))
        # expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))         # slow
        # expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))   # slow
    })


    test_that(  " fit: billing19.proteingroups ", {
        # original
            file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
            select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
            select %<>% paste0('_STD')
            object <- read_maxquant_proteingroups(file, subgroups = select)
            expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
            expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
            expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
        # subtracted
            obj <- object
            obj %<>% subtract_baseline('subgroup', 'E00_STD')
            expect_true(sumexp_contains_fit(fit_limma(object, ~1), 'limma'))
    })

    test_that( " fit: fukuda20.proteingroups ", {
        file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics') 
        object <- read_maxquant_proteingroups(file)
        expect_true(sumexp_contains_fit(fit_wilcoxon(object), 'wilcoxon'))
        expect_true(sumexp_contains_fit(fit_lm(object),       'lm'))
        expect_true(sumexp_contains_fit(fit_limma(object),    'limma'))
    })

    test_that( " fit: atkin.somascan ", {
        # Original
            file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
            object <- read_somascan(file)
            expect_true( sumexp_contains_fit(    fit_limma( object, ~ Time + Diabetes, block = 'Subject' ), 'limma'   ) )
            expect_true( sumexp_contains_fit(       fit_lm( object, ~ Time + Diabetes                    ), 'lm'      ) )
            expect_true( sumexp_contains_fit(      fit_lme( object, ~ Time + Diabetes, block = 'Subject' ), 'lme'     ) )
            expect_true( sumexp_contains_fit(     fit_lmer( object, ~ Time + Diabetes, block = 'Subject' ), 'lmer'    ) )
            expect_true( sumexp_contains_fit( fit_wilcoxon( object, ~ Time,            block = 'Subject' ), 'wilcoxon') )
        # Subtracted
            object %<>% subtract_differences(block = 'Subject', subgroupvar ='Time')
            expect_true( sumexp_contains_fit(    fit_lm(object, ~ 0 + Time), 'lm'   ))
            expect_true( sumexp_contains_fit( fit_limma(object, ~ 0 + Time), 'limma'))
    })
    
    test_that( " fit: atkin.metabolon ", {
        # read
            file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
            object <- read_metabolon(file)
        # test
            expect_true(sumexp_contains_fit(                  fit_limma( object, ~ subgroup, block = 'Subject'),  'limma'   ))
            expect_true(sumexp_contains_fit(                     fit_lm( object, ~ subgroup, block = 'Subject'),  'lm'      ))
            expect_true(sumexp_contains_fit(                    fit_lme( object, ~ subgroup, block = 'Subject'),  'lme'     ))
            expect_true(sumexp_contains_fit( suppressWarnings( fit_lmer( object, ~ subgroup, block = 'Subject')), 'lmer'    ))
            expect_true(sumexp_contains_fit(               fit_wilcoxon( object, ~ subgroup, block = 'Subject'),  'wilcoxon'))
    })
    
    test_that( " fit: mcclain21 ", {
        file  <- download_mcclain21('counts')
        sfile <- download_mcclain21('samples')
        object <- .read_rnaseq_counts(file, sfile = sfile, by.y = 'rna_id')
        expect_s4_class(object, 'SummarizedExperiment')
        
    })

        
#============================================================================
#                                                                           #
#            context(" plot_contrastogram ")                                #
#                                                                           #
#============================================================================
    
if (requireNamespace('diagram', quietly = TRUE)){
    test_that(' plot_contrastogram: fukuda20.proteingroups', {
        # Read
            file <-  system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
            object <- read_maxquant_proteingroups(file)
        # Test
            expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
})}

    
if (requireNamespace('diagram', quietly = TRUE)){
    test_that(' plot_contrastogram: billing19.proteingroups', {
        # Read
            file <-  system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
            subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
            subgroups %<>% paste0('_STD')
            object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'limma')
        # Test
            expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup', curve = 0.8))
})}


# if (requireNamespace('diagram', quietly = TRUE)){
#     test_that(' plot_contrastogram: billing16.proteingroups', {  # ratios : self-contrasts
#         # Read
#             file <- download_data('billing16.proteingroups.txt')
#             invert <- c('EM_E', 'BM_E', 'BM_EM')
#             object <- read_maxquant_proteingroups(
#                         file, invert = invert, fit = 'limma', formula = ~ 0 + subgroup)
#         # Test
#             expect_no_error(plot_contrastogram(object, subgroupvar = 'subgroup'))
# })}

