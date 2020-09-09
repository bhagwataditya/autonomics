require(magrittr)
require(autonomics)
require(autonomics.find)
glutaminase <- autonomics.data::glutaminase
glutaminase %<>% autonomics.find::default_contrasts()

glutaminase %<>% autonomics.find:::add_limma2()
plot_contrastogram(glutaminase)
plot_contrast_features(glutaminase, contrast = contrastdefs(glutaminase)[1], nplot = 1, fvars = 'BIOCHEMICAL')
plot_contrast_features(glutaminase, contrast = contrastdefs(glutaminase)[1], nplot = 1, fvars = 'BIOCHEMICAL')




autonomics.find::plot_contrast_features(glutaminase)
autonomics.find::make_ref_contrasts(glutaminase)
autonomics.find::make_ref_contrasts_within_stratum(
    autonomics.import::split_values(levels(glutaminase$subgroup), keep = FALSE), sep = '_')

