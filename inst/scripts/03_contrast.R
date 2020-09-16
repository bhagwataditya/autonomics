require(magrittr)
require(autonomics)
require(autonomics.find)
glutaminase <- autonomics.data::glutaminase
glutaminase %<>% autonomics.find::default_contrasts()

glutaminase %<>% autonomics.find:::add_limma2()
plot_contrastogram(glutaminase)
plot_contrast_features(glutaminase, contrast = contrastdefs(glutaminase)[1], nplot = 1, fvars = 'BIOCHEMICAL')
plot_contrast_features(glutaminase, contrast = contrastdefs(glutaminase)[1], nplot = 1, fvars = 'BIOCHEMICAL')


