require(magrittr)
require(ggplot2)
require(autonomics.plot)
#object <- atkin.2014::soma %>% extract(, .$time %in% c('t0', 't2', 't3'))
hypo <- atkin.2014::metabolon %>% extract(, .$time %in% c('t0', 't2', 't3'))
glutaminase <- autonomics.data::glutaminase
#object %<>% extract(, !stringi::stri_startswith_fixed(.$subgroup, 'UT'))

blanken <- function(p, axis.text.x = element_text()){
    p + guides(color = 'none', fill = 'none') + theme(axis.text.x = axis.text.x)
}

hypoplots <- list(
    blanken(plot_pca_features(hypo, dim = 1, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),
    blanken(plot_pca_samples( hypo)),
    blanken(plot_pca_features(hypo, dim = 2, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),

    blanken(plot_pls_features(hypo, dim = 1, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),
    blanken(plot_pls_samples( hypo)),
    blanken(plot_pls_features(hypo, dim = 2, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()))

glutaminaseplots <- list(
    blanken(plot_pca_features(glutaminase, dim = 1, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),
    blanken(plot_pca_samples( glutaminase)),
    blanken(plot_pca_features(glutaminase, dim = 2, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),

    blanken(plot_pls_features(glutaminase, dim = 1, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank()),
    blanken(plot_pls_samples( glutaminase)),
    blanken(plot_pls_features(glutaminase, dim = 2, n = 1, color_var = 'subgroup', fvar = 'BIOCHEMICAL'), axis.text.x = element_blank())
)

multiplot(  plotlist=hypoplots,
            layout = matrix(c(1,4,
                              2,5,
                              3,6), nrow=3, byrow=TRUE))
multiplot(  plotlist=glutaminaseplots,
            layout = matrix(c(1,4,
                              2,5,
                              3,6), nrow=3, byrow=TRUE))

multiplot(  plotlist=plotlist,
            layout = matrix(c(1,2,3,
                              4,5,6,
                              7,8,9,
                              10,11,12), nrow=4, byrow=TRUE))

multiplot(  plotlist=plotlist,
            layout = matrix(c(1,1,1,4,4,4,0,7,7,7,10,10,10,
                              2,2,2,5,5,5,0,8,8,8,11,11,11,
                              3,3,3,6,6,6,0,9,9,9,12,12,12), nrow=3, byrow=TRUE))


plot_lda_features(object, dim = 1) + guides(color = 'none')
autonomics.plot::multiplot(p3, p4, layout = matrix(c(1,0,2), nrow=1))

p3 <- plot_pca_features(object)
p4 <- plot_lda_features(object)



plot_sample_distributions(object)
plot_sample_densities(object, color = subgroup, facet = subgroup)
