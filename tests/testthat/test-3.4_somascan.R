metadata <- S4Vectors::metadata

#============================================================================
#                                                                           #
#            context(" read_somascan ")                                     #
#                                                                           #
#============================================================================

test_that( " read_somascan ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup' %in% svars(object))
    expect_true('SampleGroup' %in% svars(object))
})

test_that(  " read_somascan: groupvar = 'Subject' ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, groupvar = 'Subject')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup'    %in% svars(object))
    expect_true('Subject'    %in% svars(object))
})

test_that(  " read_somascan: pca = TRUE", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
    expect_true('sample_id~pca' %in% names(metadata(object)))
})

test_that(  " read_somascan: fit = 'limma', block = 'Subject' ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, fit = 'limma', block = 'Subject')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that(" read_somascan: fit = 'lme', block = 'Subject' ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, block = 'Subject', fit = 'lme')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lme'))))
})

test_that("read_somascan: fit = 'lmer', block = 'Subject' ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, block = 'Subject', fit = 'lmer')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lmer'))))
})

test_that("read_somascan: fit = 'wilcoxon', block = 'Subject' ", {
    file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
    object <- read_somascan(file, block = 'Subject', fit = 'wilcoxon')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
