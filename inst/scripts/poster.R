file <- download_data('atkin18.metabolon.xlsx')
object <- read_metabolon(file)
object %<>% filter_samples(subgroup != 't3')
object %<>% filter_samples(subgroup != 't1')
object %<>% autonomics::keep_connected_blocks(block = 'SUB')
object$x <- paste0(substr(object$T2D, 1, 1), substr(object$subgroup, 2, 2))
object$x %<>% factor()
levels(object$x) <- c('C1', 'C2', 'D1', 'D2')
levels(object$subgroup) <- c('t1', 't2')

object %<>% merge_fdt(.fit_limma(object, ~ subgroup,       block = 'SUB', coefs = 't2'))
object %<>% merge_fdt(.fit_limma(object, ~ subgroup + T2D, block = 'SUB', coefs = 'T2DM'))
object %<>% merge_fdt(.fit_limma(object, ~ subgroup / T2D, block = 'SUB', coefs = c('t1:T2DM', 't2:T2DM')))
object %<>% merge_fdt(.fit_limma(object, ~ T2D / subgroup, block = 'SUB', coefs = c('Control:t2', 'T2DM:t2')))
fdt0 <- .fit_limma(object, ~ T2D * subgroup, block = 'SUB', coefs = 'T2DM:t2')
names(fdt0) %<>% stri_replace_first_fixed('T2DM:t2', 'T2DM*t2')
object %<>% merge_fdt(fdt0)

sub1 <- filter_features(object,    `p~t1:T2DM~limma` < 0.05 & `p~t2:T2DM~limma` > 0.05)
sub2 <- filter_features(object,    `p~t2:T2DM~limma` < 0.05 & `p~t1:T2DM~limma` > 0.05)
sub3 <- filter_features(object, `p~Control:t2~limma` < 0.05 & `p~T2DM:t2~limma` > 0.05)
sub4 <- filter_features(object, `p~Control:t2~limma` > 0.05 & `p~T2DM:t2~limma` < 0.05)

palette <- make_colors(c('T2DM', 'Control'))
    
plot_exprs(object, x = 'x', coefs = 't2',         block = 'SUB', color = 'T2D', colorpalette = palette) # t2   effect across controls/diabetics
plot_exprs(object, x = 'x', coefs = 'T2DM',       block = 'SUB', color = 'T2D', colorpalette = palette) # T2DM effect across timepoints
plot_exprs(sub1,   x = 'x', coefs = 't1:T2DM',    block = 'SUB', color = 'T2D', colorpalette = palette) # T2DM effect at t0 but not t2
plot_exprs(sub2,   x = 'x', coefs = 't2:T2DM',    block = 'SUB', color = 'T2D', colorpalette = palette) # T2DM effect at t2 but not t0
plot_exprs(sub3,   x = 'x', coefs = 'Control:t2', block = 'SUB', color = 'T2D', colorpalette = palette) # t2   effect in Controls but not Diabetics
plot_exprs(sub4,   x = 'x', coefs = 'T2DM:t2',    block = 'SUB', color = 'T2D', colorpalette = palette) # t2   effect in Diabetics but not Controls
plot_exprs(object, x = 'x', coefs = 'T2DM*t2',    block = 'SUB', color = 'T2D', colorpalette = palette) # t2   effect in Diabetics but not Controls

palette <- c(Control = '#00B8E5', T2DM = '#FF6C91')
plot_exprs(object['methionine', ], x = 'x', coefs = NULL, block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['X - 13848', ], x = 'x', coefs = NULL, block = 'SUB', color = 'T2D', colorpalette = palette) + 
    
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['gamma-glutamylalanine', ], x = 'x', coefs = NULL, block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['docosatrienoate (22:3n3)'],   x = 'x', coefs = NULL,    block = 'SUB', color = 'T2D', colorpalette = palette) +
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['2-hydroxypalmitate', ],   x = 'x', coefs = NULL, block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['alanine', ],   x = 'x', coefs = NULL,    block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')
plot_exprs(object['docosatrienoate (22:3n3)', ], x = 'x', coefs = NULL,    block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')


plot_exprs(object['glucose', ], x = 'x', coefs = NULL,    block = 'SUB', color = 'T2D', colorpalette = palette) + 
                                theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
                                ggtitle(NULL) + guides(color = 'none')

