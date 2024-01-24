# Complex designs bring complex questions
require(devtools)
document()
load_all()
file <- download_data('atkin.metabolon.xlsx')
object <- read_metabolon(file)
object$time <- object$subgroup
object$subgroup <- paste0(object$time, '.', object$Diabetes)
object %<>% keep_connected_blocks(block = 'Subject')
object %<>% filter_samples(time != 't3')
# https://www.w3schools.com/colors/colors_picker.asp
#colorpal <- make_colors(c('t0', 't1', 't2', 't3'))
colorpal <- c(t0 = '#ff0000', t1 = '#339933', t2 = '#1D75E1')
#colorpal <- c(t0 = '#ff0000', t1 = '#339933', t2 = '#1D75E1', t3 = '#9900ff')
plotdt <- sumexp_to_longdt(object['glucose',], svars = c('subgroup', 'Subject', 'Diabetes', 'Time'))
plotdt[time=='t0', x := 0]
plotdt[time=='t1', x := 1]
plotdt[time=='t2', x := 2]
#plotdt[time=='t3', x := 3]

p <- ggplot2::ggplot(plotdt) + theme_bw() + 
     geom_point(aes(x = x, y = value, color = time, alpha = Diabetes, shape = Diabetes), size = 3) + 
     geom_line( aes(x = x, y = value, color = time, group = Subject, alpha = Diabetes), show.legend = FALSE) + 
     scale_alpha_manual(values = c(Control = 0.35, T2DM = 1)) + 
     scale_color_manual(values = colorpal) + 
     theme(panel.grid = element_blank())
p <- p + ggtitle('Atkin 2013 : Hypoglycemia in Type 2 Diabetes (NCT02205996)')
p <- p + xlab('time') + ylab('Glucose')
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + xlim(c(-1, 3.5))
p <- p + ylim(28, 31.5)
p <- p + annotate('text', x = -0.2, y = 29.7,  label = 'Diabetic',                                 hjust = 1)
p <- p + annotate('text', x = -0.2, y = 29.2,  label = 'Control',                                  hjust = 1, alpha = 0.5)
p <- p + annotate('text', x = -0.4, y = 30.20, label = 't0: baseline\nOvernight fast',             hjust = 0, color = '#ff0000')
p <- p + annotate('text', x =  0.7, y = 29.80, label = 't1: euglycaemia\nglc = 90 mg/dL for 1h',   hjust = 0, color = '#339933')
p <- p + annotate('text', x =  1.7, y = 29.40, label = 't2: hypoglycaemia\nglc = 50 mg/dL for 1h', hjust = 0, color = '#1D75E1')
#p <- p + annotate('text', x =  2.5,   y = 30.20, label = 't3: post-study\nOvernight fast',           hjust = 0, color = '#1D75E1')
p
#p<- p + scale_x_continuous(labels = c('', 't0', 't1', 't2', 't3'))

p <- p + annotate('text', x = 3.5,  y = 31.40,   label = "~ diabetes + time", hjust = 0)
#p <- p + annotate('text', x = 1.2,  y = 31.40,   label =    't[1]-t[0]',     hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 2.2,  y = 31.40,   label =    't[2]-t[0]',     hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.2,  y = 31.40,   label =    't[3]-t[0]',     hjust = 1, parse = TRUE) 
#p<- p + annotate('text', x = 3.8,  y = 31.40,   label =    't[3]-t[0]',     hjust = 1, parse = TRUE)

p <- p + annotate('text', x = 3.5,  y = 31.20,   label = "~ diabetes / time", hjust = 0)
p <- p + annotate('text', x = 1.2,  y = 31.20,   label =    'C:t[1]-t[0]',   hjust = 1, parse = TRUE) 
p <- p + annotate('text', x = 2.2,  y = 31.20,   label =    'C:t[2]-t[0]',   hjust = 1, parse = TRUE)
#p<- p + annotate('text', x = 3.8,  y = 31.20,   label =    'C:t[3]-t[0]',   hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.2,  y = 31.20,   label =    'C:t[3]-t[0]',   hjust = 1, parse = TRUE)

p <- p + annotate('text', x = 1.2,  y = 31.00,   label =    'D:t[1]-t[0]',   hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 2.2,  y = 31.00,   label =    'D:t[2]-t[0]',   hjust = 1, parse = TRUE)
#p<- p + annotate('text', x = 3.8,  y = 31.00,   label =    'D:t[3]-t[0]',   hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.2,  y = 31.00,   label =    'D:t[3]-t[0]',   hjust = 1, parse = TRUE)

p <- p + annotate('text', x = 3.5,  y = 30.80,   label = "~ diabetes * time", hjust = 0)
p <- p + annotate('text', x = 1.2,  y = 30.80,   label =    'D-C:t[1]-t[0]', hjust = 1, parse = TRUE)
#p<- p + annotate('text', x = 3.8,  y = 30.80,   label =    'D-C:t[3]-t[0]', hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 2.2,  y = 30.80,   label =    'D-C:t[2]-t[0]', hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.2,  y = 30.80,   label =    'D-C:t[3]-t[0]', hjust = 1, parse = TRUE)

p <- p + annotate('text', x = 3.5,  y = 28.20,   label = "~ time / diabetes", hjust = 0)
p <- p + annotate('text', x = 0.33,  y = 28.20,   label =    't[0]:D-C',      hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 1.33,  y = 28.20,   label =    't[1]:D-C',      hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 2.33,  y = 28.20,   label =    't[2]:D-C',      hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.33,  y = 28.20,   label =    't[3]:D-C',      hjust = 1, parse = TRUE)

p <- p + annotate('text', x = 3.5,  y = 28.00,   label = "~ time * diabetes", hjust = 0)
p <- p + annotate('text', x = 1.33,  y = 28.00,   label =    't[1]-t[0]:D-C', hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 2.33,  y = 28.00,   label =    't[2]-t[0]:D-C', hjust = 1, parse = TRUE)
p <- p + annotate('text', x = 3.33,  y = 28.00,   label =    't[3]-t[0]:D-C', hjust = 1, parse = TRUE)
p

# subgroups : common part
require(data.table)
require(ggplot2)
arrowlength <- 0.08
lw <- 1
bg  <- 'white'
fg1 <- '#ff0000'
fg2 <- '#d22d2d' # 'black' # 'white'
fg3 <- '#9f6060' # 'black' # 'white'
font <- 'sans'
fontsize <- 4.5
parse <- TRUE
x1 <- 2
x2 <- 4
x3 <- 6
x4 <- 8
x5 <- 10
x6 <- 12
y1 <- 1
y2 <- 2
y3 <- 3
y4 <- 4
y5 <- 5
y6 <- 6

palette <- c(t1=fg1, t2=fg2, t3=fg3)
dt <- data.table(xlab     = c('t[1]*.C', 't[1]*.D', 't[2]*.C', 't[2]*.D', 't[3]*.C', 't[3]*.D'),
                 x        = c( x1,        x2,        x3,        x4,        x5,        x6),
                 y        = c( y1,        y2,        y3,        y4,        y5,        y6), 
                 diabetes = c('C',       'D',       'C',        'D',      'C',       'D'), 
                 time     = c('t1',      't1',      't2',       't2',     't3',      't3'))
p <- ggplot(dt) + theme_bw()
p <- p + xlim(c( 1, x6 + 1.5))
p <- p + ylim(c(-1.5, y6 + 0.5))
p <- p + geom_text(aes(x = x, y = -1, label = xlab, colour = time), parse = parse, size = fontsize, family = font) + scale_color_manual(values = palette)
p <- p + geom_segment(aes(x = x-0.5, xend = x+0.5, y = y, yend = y, colour = time), linewidth = lw) + scale_color_manual(values = palette)
p <- p + theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_rect(color = fg, linewidth = lw))

# ~ time / diabetes
p1 <- p
p1 <- p1 + annotate('text', label = '~ time / diabetes', x = x1, y = y6, hjust = 0, colour = fg2, size = 5.5, family = font)
p1 <- p1 + annotate('segment', x = x1, xend = x1, y = 0,  yend = y1,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg1, linewidth = lw)
p1 <- p1 + annotate('segment', x = x2, xend = x2, y = y1, yend = y2,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg1, linewidth = lw)
p1 <- p1 + annotate('segment', x = x3, xend = x3, y = 0,  yend = y3,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg2, linewidth = lw)
p1 <- p1 + annotate('segment', x = x4, xend = x4, y = y3, yend = y4,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg2, linewidth = lw)
p1 <- p1 + annotate('segment', x = x5, xend = x5, y = 0,  yend = y5,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg3, linewidth = lw)
p1 <- p1 + annotate('segment', x = x6, xend = x6, y = y5, yend = y6,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg3, linewidth = lw)
p1 <- p1 + annotate('text',    x = x1+0.2,        y = y1-0.6,   label = 't[1]',   hjust = 0,   colour = fg1, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('text',    x = x2+0.2,        y = y2-0.6,   label = 't[1]:D', hjust = 0,   colour = fg1, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('text',    x = x3+0.2,        y = y3-0.6,   label = 't[2]',   hjust = 0,   colour = fg2, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('text',    x = x4+0.2,        y = y4-0.6,   label = 't[2]:D', hjust = 0,   colour = fg2, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('text',    x = x5+0.2,        y = y5-0.6,   label = 't[3]',   hjust = 0,   colour = fg3, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('text',    x = x6+0.2,        y = y6-0.6,   label = 't[3]:D', hjust = 0,   colour = fg3, parse = parse, size = fontsize, family = font)
p1 <- p1 + annotate('segment', x = x1-0.5, xend = x2+0.5, y = y1, yend = y1, linetype = 2, colour = fg1, linewidth = 0.5)
p1 <- p1 + annotate('segment', x = x3-0.5, xend = x4+0.5, y = y3, yend = y3, linetype = 2, colour = fg2, linewidth = 0.5)
p1 <- p1 + annotate('segment', x = x5-0.5, xend = x6+0.5, y = y5, yend = y5, linetype = 2, colour = fg3, linewidth = 0.5)
p1 <- p1 + guides(color = 'none')


# ~ diabetes / time
p2 <- p
p2 <- p2 + annotate('text', label = '~ diabetes / time', x = x1, y = y6, hjust = 0, colour = fg2, size = 5.5, family = font)
p2 <- p2 + annotate('segment', x = x1, xend = x1, y = 0,  yend = y1,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg1, linewidth = lw)
p2 <- p2 + annotate('segment', x = x2, xend = x2, y = 0,  yend = y2,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg1, linewidth = lw)
p2 <- p2 + annotate('segment', x = x3, xend = x3, y = y1, yend = y3,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg2, linewidth = lw)
p2 <- p2 + annotate('segment', x = x4, xend = x4, y = y2, yend = y4,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg2, linewidth = lw)
p2 <- p2 + annotate('segment', x = x5, xend = x5, y = y1, yend = y5,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg3, linewidth = lw)
p2 <- p2 + annotate('segment', x = x6, xend = x6, y = y2, yend = y6,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = fg3, linewidth = lw)
p2 <- p2 + annotate('segment', x = x1-0.5, xend = x6+0.5, y = y1, yend = y1, linetype = 2, colour = fg1, linewidth = 0.5)
p2 <- p2 + annotate('segment', x = x2-0.5, xend = x6+0.5, y = y2, yend = y2, linetype = 2, colour = fg1, linewidth = 0.5)
p2 <- p2 + annotate('text',    x = x1+0.2,         y = y1-0.6, label = 'C',      hjust = 0,   colour = fg1, parse = parse, size = fontsize, family = font)
p2 <- p2 + annotate('text',    x = x2+0.2,         y = y2-0.6, label = 'D',      hjust = 0,   colour = fg1, parse = parse, size = fontsize, family = font)
p2 <- p2 + annotate('text',    x = x3+0.2,         y = y3-0.6, label = 'C:t[1]', hjust = 0,   colour = fg2, parse = parse, size = fontsize, family = font)
p2 <- p2 + annotate('text',    x = x4+0.2,         y = y4-0.6, label = 'D:t[1]', hjust = 0,   colour = fg2, parse = parse, size = fontsize, family = font)
p2 <- p2 + annotate('text',    x = x5+0.2,         y = y5-0.6, label = 'C:t[2]', hjust = 0,   colour = fg3, parse = parse, size = fontsize, family = font)
p2 <- p2 + annotate('text',    x = x6+0.2,         y = y6-0.6, label = 'D:t[2]', hjust = 0,   colour = fg3, parse = parse, size = fontsize, family = font)
p2 <- p2 + guides(color = 'none')

gridExtra::grid.arrange(p1 + theme(axis.text = element_blank(), axis.ticks = element_blank()), 
                        p2 + theme(axis.text = element_blank(), axis.ticks = element_blank()), nrow = 2)


# Alternative order
dt <- data.table(xlab     = c('C*.t[1]', 'C*.t[2]', 'C*.t[3]', 'D*.t[1]', 'D*.t[2]', 'D*.t[3]'),
                 x        = c( x1,        x2,        x3,        x4,        x5,        x6),
                 y        = c( y1,        y2,        y3,        y4,        y5,        y6), 
                 diabetes = c('C',       'C',       'C',        'D',      'D',       'D'), 
                 time     = c('t1',      't2',      't3',       't1',     't2',      't3'))
p3 <- ggplot(dt) + theme_bw()
p3 <- p3 + xlim(c( 1, x6 + 1.5))
p3 <- p3 + ylim(c(-1.5, y6 + 0.5))
p3 <- p3 + geom_text(aes(x = x, y = -1, label = xlab, colour = diabetes), parse = parse, size = fontsize, family = font) + scale_color_manual(values = c(C = 'red', D = 'darkred'))
p3 <- p3 + geom_segment(aes(x = x-0.5, xend = x+0.5, y = y, yend = y, colour = diabetes), linewidth = lw)
p3 <- p3 + theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_rect(color = fg, linewidth = lw))
p3 <- p3 + annotate('text', label = '~ diabetes / time', x = 1, y = y6, hjust = 0, colour = fg2, size = 5.5, family = font)
p3 <- p3 + annotate('segment', x = x1, xend = x1, y =  0, yend = y1,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'red',     linewidth = lw)
p3 <- p3 + annotate('segment', x = x2, xend = x2, y = y1, yend = y2, arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'red',     linewidth = lw)
p3 <- p3 + annotate('segment', x = x3, xend = x3, y = y1, yend = y3,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'red',     linewidth = lw)
p3 <- p3 + annotate('segment', x = x4, xend = x4, y =  0, yend = y4,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'darkred', linewidth = lw)
p3 <- p3 + annotate('segment', x = x5, xend = x5, y = y4, yend = y5,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'darkred', linewidth = lw)
p3 <- p3 + annotate('segment', x = x6, xend = x6, y = y4, yend = y6,  arrow = arrow(length = unit(arrowlength, 'inches')), colour = 'darkred', linewidth = lw)
p3 <- p3 + annotate('segment', x = x1-0.5, xend = x3+0.5, y = y1, yend = y1, linetype = 2, colour = 'red', linewidth = 0.5)
p3 <- p3 + annotate('segment', x = x4-0.5, xend = x6+0.5, y = y4, yend = y4, linetype = 2, colour = 'darkred', linewidth = 0.5)
p3 <- p3 + annotate('text',    x = x1+0.2,         y = y1-0.5, label = 'C',      hjust = 0,   colour = 'red',     parse = parse, size = fontsize, family = font)
p3 <- p3 + annotate('text',    x = x2+0.2,         y = y2-0.5, label = 'C:t[2]', hjust = 0,   colour = 'red',     parse = parse, size = fontsize, family = font)
p3 <- p3 + annotate('text',    x = x3+0.2,         y = y3-0.5, label = 'C:t[3]', hjust = 0,   colour = 'red',     parse = parse, size = fontsize, family = font)
p3 <- p3 + annotate('text',    x = x4+0.2,         y = y1-0.5, label = 'D',      hjust = 0,   colour = 'darkred', parse = parse, size = fontsize, family = font)
p3 <- p3 + annotate('text',    x = x5+0.2,         y = y5-0.5, label = 'D:t[2]', hjust = 0,   colour = 'darkred', parse = parse, size = fontsize, family = font)
p3 <- p3 + annotate('text',    x = x6+0.2,         y = y6-0.5, label = 'D:t[3]', hjust = 0,   colour = 'darkred', parse = parse, size = fontsize, family = font)
p3 <- p3 + guides(color = 'none')
p3 + theme(axis.text = element_blank(), axis.ticks = element_blank())