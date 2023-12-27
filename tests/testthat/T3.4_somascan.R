context('read_somascan')       # late - early : 32 down, 75 up
metadata <- S4Vectors::metadata

test_that( "read_somascan(file) ", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup' %in% svars(object))
    expect_false('SampleGroup' %in% svars(object))
})

test_that(  "read_somascan(file, subgroupvar = 'Subject')", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, subgroupvar = 'Subject')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup'    %in% svars(object))
    expect_false('Subject'    %in% svars(object))
    expect_true('SampleGroup' %in% svars(object))
})

test_that(  "read_somascan(file, pca = TRUE)", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
    expect_true('sample_id~pca' %in% names(metadata(object)))
})

test_that(  "read_somascan(file, fit = 'limma', block = 'Subject', plot = TRUE)", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, fit = 'limma', block = 'Subject', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that("read_somascan(file, fit = 'lme', block = 'Subject', plot = TRUE)", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, block = 'Subject', fit = 'lme', plot = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lme'))))
})

test_that("read_somascan(file,  fit = 'lmer', block = 'Subject')", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, block = 'Subject', fit = 'lmer', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lmer'))))
})

test_that("read_somascan(file, fit = 'wilcoxon', block = 'Subject', plot = TRUE)", {
    file <- download_data('atkin.somascan.adat')
    object <- read_somascan(file, block = 'Subject', fit = 'wilcoxon', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
