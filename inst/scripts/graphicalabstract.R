require(devtools)
load_all()
file <- download_data('atkin.somascan.adat')
object <- read_somascan(file)
object %<>% fit_limma(block = 'Subject', codingfun = code_diff)
plot_exprs(object, coef = 't1-t0', block = 'Subject', n = 2)

require(ggplot2)

p <- ggplot() + theme_void() 
p <- p + ggforce::geom_circle(aes(x0 = 1, y0 = 1, r = 2.5))
p <- p + annotate('rect', xmin = 0,   ymin = 1.0, xmax = 1.0, ymax = 2.0, fill = 'white', colour = 'black')
p <- p + annotate('rect', xmin = 1,   ymin = 1,   xmax = 2,   ymax = 2,   fill = 'white', colour = 'black')
p <- p + annotate('rect', xmin = 0,   ymin = 0,   xmax = 1,   ymax = 1,   fill = 'white', colour = 'black')

p <- p + annotate('text', label = 'assays()', x = 0.5, y = 1.5, family = 'mono')
p <- p + annotate('text', label = 'fdt()',    x = 1.5, y = 1.5, family = 'mono')
p <- p + annotate('text', label = 'sdt()',    x = 0.5, y = 0.5, family = 'mono')

p <- p + annotate('text', x = c(0.2, 0.4, 1.2, 1.4), y = 2.2,        label = c('C0.t0', 'C0.t1', 'pathway', 'mass'), angle = 90, hjust = 0)
p <- p + annotate('text', x = -0.05, y = c(1.9, 1.6), label = c('methionine', 'palmitate'), hjust = 1)
p <- p + annotate('text', x = -0.05, y = c(0.9, 0.6, 0.3), label = c('subject', 'time', 't2d'), hjust = 1)

p <- p + annotate('segment', x = 0.1, y = 2.1, xend = 1.1, yend = 2.1) +
         annotate('segment', x = 0,   y = 2,   xend = 0.1, yend = 2.1) +
         annotate('segment', x = 1,   y = 2,   xend = 1.1, yend = 2.1) +
         annotate('segment', x = 1.1, y = 2.1, xend = 1.1, yend = 2.0)
p <- p + annotate('text',   x = 1, y = -0.3, label = 'object', hjust = 0.5, size = 4.5, family = 'mono')
#p <- p + annotate('text',   x = 1, y = -0.8, label = 'SummarizedExperiment', hjust = 0.5, size = 4.5)

dt <- data.table(
    label = c(
              'METABOLITE\nread_metabolon()',
              'PROTEIN : AFFINITY\nread_somascan()\nread_olink()', 
              'PROTEIN : DIA\nread_diann()\nread_fragpipe()', 
              'PROTEIN : DDA\nread_proteingroups()\nread_phosphosites()', 
              'RNA : ARRAY\nread_affymetrix()',
              'RNA : NGS\nread_counts()\nread_bams()'
              ),
    platform = c('MET', 'AFFPRO', 'DIAPRO', 'DDAPRO', 'ARRRNA', 'NGSRNA'),
    degrees  = (0:5)*180/5
)
dt[, platform := factor(platform, unique(platform))]
dt[, radians := 2*pi/360*degrees]
dt[, xend := 1 + 2.5*cos(radians)]
dt[, yend := 1 + 2.5*sin(radians)]
dt[, x    := 1 + 4.5*cos(radians)]
dt[, y    := 1 + 4.5*sin(radians)]
p <- p + geom_segment(data = dt, aes(x = x, xend = xend, y = y, yend = yend), arrow = arrow())
p <- p + geom_label(data = dt, aes(x = x, y = y, label = label, fill = platform))
p <- p + guides(fill = 'none') + xlim(1+c(-6, +6)) + ylim(1+c(-6, +6))
p <- p + scale_fill_manual(values = c(NGSRNA = "#F8766D", ARRRNA = "#F564E3", DDAPRO = "#619CFF", DIAPRO = "#00BFC4", AFFPRO = "#00BA38", MET = "yellow"  ))

p + annotate('text', x = -3, y = -2.0, hjust = 0, family = 'mono', label = "object %<>% fit_lm(       ~time)                                           # classic") + 
    annotate('text', x = -3, y = -2.3, hjust = 0, family = 'mono', label = "            fit_lme(      ~time, block = 'subject')                        # random effects") + 
    annotate('text', x = -3, y = -2.6, hjust = 0, family = 'mono', label = "            fit_lmer(     ~time, block = 'subject')                        # advanced random effects") + 
    annotate('text', x = -3, y = -2.9, hjust = 0, family = 'mono', label = "            fit_wilcoxon( ~time, block = 'subject')                        # unparametric") +
    
    annotate('text', x = -3, y = -3.5, hjust = 0, family = 'mono', label = "            fit_limma(    ~time, block = 'subject'                         # large-scale") + 
    annotate('text', x = -3, y = -3.8, hjust = 0, family = 'mono', label = "            fit_limma(    ~time, block = 'subject', codingfun = contr.treatment.explicit   # alternative codings") + 
    annotate('text', x = -3, y = -4.1, hjust = 0, family = 'mono', label = "            fit_limma(  ~0+time, block = 'subject', contrasts = 't1-t0')   # simplified contrasts")


fdt(object) %<>% extract(, 1:2)
object %<>% fit_lm(      codingfun = contr.treatment.explicit, coefs = 't3-t0')
object %<>% fit_lme(     codingfun = contr.treatment.explicit, coefs = 't3-t0', block = 'Subject')
object %<>% fit_lmer(    codingfun = contr.treatment.explicit, coefs = 't3-t0', block = 'Subject')
object %<>% fit_limma(   codingfun = contr.treatment.explicit, coefs = 't3-t0', block = 'Subject')
object %<>% fit_wilcoxon(block = 'Subject')



fdt(object) %<>% extract(, c(1:14, 17, 23, 20))

plot_contrast_venn(is_sig(object, fit = c('limma', 'lme',  'wilcoxon'), contrast = 't3-t0'))

plot(-log10(fdt(object)$`p~t3-t0~limma`), -log10(fdt(object)$`p~t3-t0~wilcoxon`))
plot(-log10(fdt(object)$`p~t3-t0~limma`), -log10(fdt(object)$`p~t3-t0~lme`))





