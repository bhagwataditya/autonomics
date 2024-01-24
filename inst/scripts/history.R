devtools::document()
set.seed(42)

# normal distribution
dt <- data.table(y = rnorm(10000, 0))
p <- ggplot(dt) + 
     theme_minimal() + 
     theme(panel.grid.major = element_blank(), 
          panel.grid.minor  = element_blank(), 
          axis.title.x      = element_blank(), 
          axis.title.y      = element_blank(),
          axis.text.x       = element_blank(),
          axis.text.y       = element_blank(),
     ) + 
     geom_density(aes(y = y),   color = 'white') + 
     geom_vline(xintercept = 0, color = '#F8766D')
p <- p + annotate('point', y =  0,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = -2,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = +2,   x = 0, color ='#F8766D')
p <- p + annotate('point', y =  0.4, x = 0, color ='#F8766D')
p <- p + annotate('point', y = -0.4, x = 0, color ='#F8766D')
p <- p + annotate('point', y =  0.2, x = 0, color ='#F8766D')
p <- p + annotate('point', y = -0.2, x = 0, color ='#F8766D')
p <- p + annotate('point', y =  1,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = -1,   x = 0, color ='#F8766D')
p <- p + geom_density(aes(y = y),           color ='#F8766D')

# t distribution
p <- ggplot() + 
theme_minimal() + 
theme(panel.grid.major = element_blank(), 
  panel.grid.minor  = element_blank(), 
  axis.title.x      = element_blank(), 
  axis.title.y      = element_blank(),
  axis.text.x       = element_blank(),
  axis.text.y       = element_blank(),
)
p + geom_density(data = data.table(x = rnorm(10000)),   mapping = aes(y = x), color = '#F8766D') + 
    geom_density(data = data.table(x =    rt(10000,4)/sqrt(10)), mapping = aes(y = x), color = '#00BFC4') # 5 samples  # #00BFC4
p + geom_density(data = data.table(x = rnorm(10000)),   mapping = aes(y = x), color = '#F8766D') + 
    geom_density(data = data.table(x =    rt(10000,4)/sqrt(10)), mapping = aes(y = x), color = 'white') # 5 samples  # #00BFC4
p + geom_density(data = data.table(x = rnorm(10000)),   mapping = aes(y = x), color = 'white') + 
    geom_density(data = data.table(x =    rt(10000,4)/sqrt(10)), mapping = aes(y = x), color = '#00BFC4') # 5 samples  # #00BFC4


dt <- rbind(data.table(batch = '1', y = rnorm(10000, 13, sd = 0.5)), 
            data.table(batch = '2', y = rnorm(10000, 15, sd = 0.5)))
p <- ggplot(dt) + 
     theme_minimal() + 
     theme(panel.grid.major = element_blank(), 
          panel.grid.minor  = element_blank(), 
          axis.title.x      = element_blank(), 
          axis.title.y      = element_blank(),
          axis.text.x       = element_blank(),
          axis.text.y       = element_blank(),
     ) + 
     geom_density(aes(y = y, linetype = batch),   color = 'white') + 
     geom_vline(xintercept = 0, color = '#F8766D')
p <- p + annotate('point', y = 13,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = 12,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = 14,   x = 0, color ='#F8766D')
p <- p + annotate('point', y = 13.2, x = 0, color ='#F8766D')
p <- p + annotate('point', y = 12.8, x = 0, color ='#F8766D')
p <- p + annotate('point', y = 13.1, x = 0, color ='#F8766D')
p <- p + annotate('point', y = 13.5, x = 0, color ='#F8766D')
p <- p + annotate('point', y = 12.5, x = 0, color ='#F8766D')
p <- p + annotate('point', y = 15,   x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 14,   x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 16,   x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 15.2, x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 14.8, x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 15.1, x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 15.5, x = 0, color ='#F8766D', shape = 17)
p <- p + annotate('point', y = 14.5, x = 0, color ='#F8766D', shape = 17)
p <- p + geom_density(aes(y = y, linetype = batch),   color = '#F8766D')

# Two sample t test
#------------------

# Single curve
    set.seed(1)
    y0 <- seq(-4, +4, length.out = 100)
    plotdt <- rbind(data.table(y = y0 + 2, x = -dnorm(y0), type = 'sample'))
    p <- ggplot(plotdt, aes(x = x, y = y)) + theme_bw()  #+ facet_grid(rows = vars(type), scales = 'free_y')
    p <- p + theme(panel.grid.major.y = element_blank(), 
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank())
    p <- p + scale_y_continuous(breaks = c(0,2,5,8))
    p <- p + guides(color = 'none')
    p <- p + ylab(expression(Glc)) + xlab(expression(Frequency))
    p <- p + geom_path(aes(linetype = type), color = '#F8766D')

# Add triplets
    # p <- p + xlim(-0.4, 0.13)
    pointdt <- rbind(data.table(y = rnorm(3,2), x = rep(-0.07, 3), i = rep(1,3)), 
                     data.table(y = rnorm(3,2), x = rep(-0.04, 3), i = rep(2,3)), 
                     data.table(y = rnorm(3,2), x = rep(-0.01, 3), i = rep(3,3)), 
                     data.table(y = rnorm(3,2), x = rep( 0.02, 3), i = rep(4,3)), 
                     data.table(y = rnorm(3,2), x = rep( 0.05, 3), i = rep(5,3)), 
                     data.table(y = rnorm(3,2), x = rep( 0.08, 3), i = rep(6,3)))
    meandt <- pointdt[, .( y   = mean(y), 
                          ymin = mean(y) - sd(y)/2, 
                          ymax = mean(y) + sd(y)/2, 
                          x    = unique(x)),               by = c('i') ]
    (pp <- p                                                                      + geom_line(data = pointdt, color = '#F8766D', aes(group = i)) + geom_point(data = pointdt, color = '#F8766D') + guides(linetype = 'none'))
    (pp <- p + geom_point(data =  meandt, size = 3, shape = 3, color = '#F8766D') + geom_line(data = pointdt, color = '#F8766D', aes(group = i)) + geom_point(data = pointdt, color = '#F8766D') + guides(linetype = 'none'))
    t0 <- seq(-4, 4, length.out = 100)
    tdt <- data.table( y = 2 +t0, x = -dnorm(t0, sd = 1/sqrt(3)), type = 'samplemean(3)')
    pp  + geom_path(data = tdt, aes(x = x, y = y, linetype = type), color = '#F8766D')
    
# Add decets
    pointdt <- rbind(data.table(y = rnorm(10,2), x = rep(-0.07, 10), i = rep(1,10)), 
                     data.table(y = rnorm(10,2), x = rep(-0.04, 10), i = rep(2,10)), 
                     data.table(y = rnorm(10,2), x = rep(-0.01, 10), i = rep(3,10)), 
                     data.table(y = rnorm(10,2), x = rep( 0.02, 10), i = rep(4,10)), 
                     data.table(y = rnorm(10,2), x = rep( 0.05, 10), i = rep(5,10)), 
                     data.table(y = rnorm(10,2), x = rep( 0.08, 10), i = rep(6,10)))
    meandt <- pointdt[, .( y   = mean(y), 
                          ymin = mean(y) - sd(y)/2, 
                          ymax = mean(y) + sd(y)/2, 
                          x    = unique(x)),               by = c('i') ]
    (pp <- p                                                                      + geom_line(data = pointdt, color = '#F8766D', aes(group = i)) + geom_point(data = pointdt, color = '#F8766D') + guides(linetype = 'none'))
    (pp <- p + geom_point(data =  meandt, size = 3, shape = 3, color = '#F8766D') + geom_line(data = pointdt, color = '#F8766D', aes(group = i)) + geom_point(data = pointdt, color = '#F8766D') + guides(linetype = 'none'))
    t0 <- seq(-4, 4, length.out = 100)
    tdt <- data.table( y = 2 +t0, x = -dnorm(t0, sd = 1/sqrt(10)), type = 'samplemean(3)')
    pp  + geom_path(data = tdt, aes(x = x, y = y, linetype = type), color = '#F8766D')

# p value
    y0 <- seq(-4,+4, length.out = 100)
    x <- dt(y0, df = 3)
    ggplot(data.table(x = -x, y = y0), aes(x = x, y = y)) + theme_bw() + geom_path(color = '#F8766D') + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab(expression(density)) + ylab(expression(t))

# Two curves - very near
y0 <- seq(-4, +4, length.out = 100)
x1 <- dnorm(y0)
x2 <- dnorm(y0)
plotdt <- rbind(data.table(x = - x1, y = 2.0 + y0, group = 'group1'), 
                data.table(x = - x2, y = 2.2 + y0, group = 'group2'))
pointdt <- rbind(data.table(y = rnorm(3, 2.0), x = -0.05, group = 'group1'), 
                 data.table(y = rnorm(3, 2.2), x =  0,    group = 'group2'))
meandt <- pointdt[, .(x = unique(x), y = mean(y)), by = 'group']
ggplot(plotdt) + 
geom_path(aes(x = x, y = y, color = group)) + theme_bw() + 
theme(panel.grid = element_blank()) + 
guides(color = FALSE) + 
geom_point(data = pointdt, aes(x =x, y = y, color = group)) + 
geom_line(data = pointdt, aes(x = x, y = y, color = group, group = group)) + 
geom_point(data = meandt, aes(x = x, y = y, color = group), shape = 3)

        
# Two curves - nearer
    y0 <- seq(-4, +4, length.out = 100)
    x1 <- dnorm(y0)
    x2 <- dnorm(y0)
    plotdt <- rbind(data.table(x = - x1, y = 2 + y0, group = 'group1'), 
                    data.table(x = - x2, y = 5 + y0, group = 'group2'))
    pointdt <- rbind(data.table(y = rnorm(3, 2), x = 0, group = 'group1'), 
                     data.table(y = rnorm(3, 5), x = 0, group = 'group2'))
    meandt <- pointdt[, .(x = unique(x), y = mean(y)), by = 'group']
    ggplot(plotdt) + 
    geom_path(aes(x = x, y = y, color = group)) + theme_bw() + 
    theme(panel.grid = element_blank()) + 
    guides(color = FALSE) + 
    geom_point(data = pointdt, aes(x =x, y = y, color = group)) + 
    geom_line(data = pointdt, aes(x = x, y = y, color = group, group = group)) + 
    geom_point(data = meandt, aes(x = x, y = y, color = group), shape = 3)
    
        
# Two curves
    y0 <- seq(-4, +4, length.out = 100)
    plotdt <- rbind(
        data.table(y = y0 + 2, x = -dnorm(y0), group = 'Control',  type = 'sample',     block = 1), 
        data.table(y = y0 + 5, x = -dnorm(y0), group = 'Diabetic', type = 'sample',     block = 2))
    p <- ggplot(plotdt, aes(x = x, y = y, color = group, group = group)) + theme_bw()  #+ facet_grid(rows = vars(type), scales = 'free_y')
    p <- p + theme(panel.grid.major.y = element_blank(), 
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank())
    p <- p + scale_y_continuous(breaks = c(0,2,5,8))
    p <- p + guides(color = 'none')
    p <- p + ylab(expression(Glc)) + xlab(expression(Frequency))
    (p <- p + geom_path(linewidth = 0.8))
    (p <- p + xlim(-0.4, 0.13))
    # p <- p + annotate('text', y = 2, x = -0.25, label = 'Control',  color = '#F8766D', hjust = 0)
    # p <- p + annotate('text', y = 5, x = -0.25, label = 'Diabetic', color = '#00BFC4', hjust = 0)

# Samples
    pointdt <- rbind(data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(1,6)), 
                     data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(2,6)), 
                     data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(3,6)), 
                     data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(4,6)), 
                     data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(5,6)), 
                     data.table(y = c(rnorm(3,2), rnorm(3,5)), x = c(rep(0.08, 3), rep(0.08, 3)), group = c(rep('Control', 3), rep('Diabetic', 3)), i = rep(6,6)))
    meandt <- pointdt[, .( y   = mean(y), 
                          ymin = mean(y) - sd(y)/2, 
                          ymax = mean(y) + sd(y)/2, 
                          x    = 0.17),               by = c('group', 'i') ]
    j <- 1
    (pp <- p  + geom_point(   data = pointdt[i==j]) + geom_line(data = pointdt[i==j]))
    pp <- pp + xlim(c(-0.4, 0.50))
    (pp <- pp + geom_crossbar(data = meandt[ i==j], aes(y = y, ymin = ymin, ymax = ymax), width = 0.02))
    pp + annotate('segment', x = 0.35, xend = 0.35, y = meandt[i==j, min(ymax)]+0.2, yend = meandt[i==j, max(ymin)],     arrow = arrow(ends = 'both', length = unit(.2, 'cm'))) + 
         annotate('segment', x = 0.29, xend = 0.45, y = meandt[i==j, min(ymax)],     yend = meandt[i==j, min(ymax)]) + 
         annotate('segment', x = 0.35, xend = 0.35, y = meandt[i==j, min(ymin)+0.1], yend = meandt[i==j, min(ymax)-0.2], arrow = arrow(ends = 'both', length = unit(.2, 'cm'))) + 
         annotate('segment', x = 0.33, xend = 0.37, y = meandt[i==j, min(ymin)],     yend = meandt[i==j, min(ymin)]) + 
         annotate('point',   x = 0.33,              y = meandt[i==j, min(ymin)]-0.4) + 
         annotate('point',   x = 0.35,              y = meandt[i==j, min(ymin)]-0.4) + 
         annotate('point',   x = 0.37,              y = meandt[i==j, min(ymin)]-0.4) + 
         annotate('text',    x = 0.27, y = meandt[i==j, min(ymax)], label = 't = ')  + 
         annotate('text',    x = 0.37, y = meandt[i==j, mean(c(min(ymax), max(ymin)))], label = 'between', hjust = 0) + 
         annotate('text',    x = 0.38, y = meandt[i==j, mean(c(min(ymin), min(ymax)))], label = 'within',  hjust = 0) + 
         annotate('text',    x = 0.39, y = meandt[i==j, min(ymin)-0.3                ], label = 'sqrt(n)', hjust = 0, parse = TRUE)

    j <- 1; (p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]))
    j <- 1;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]) + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])
    j <- 2; (p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]))
    j <- 2;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]) + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])
    j <- 3;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j])
    j <- 3;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]) + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])
    j <- 4;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j])
    j <- 4;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]) + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])
    j <- 5;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j])
    j <- 5;  p + geom_point( data = pointdt[i==j]) + geom_line(data = pointdt[i==j]) + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])

# tvalues
    tdt <- pointdt[, 
              .(mean1 = mean(.SD[group == 'Control', x ]), 
                mean1 = mean(.SD[group == 'T2D',     x ]), 
                sd1   =   sd(.SD[group == 'Control', x ]), 
                sd2   =   sd(.SD[group == 'T2D',     x ]), 
                n1    = 3,
                n2    = 3), by = 'i']

    tdt[, t := (mean2 - mean1) / sqrt((sd1^2/n1 + sd2^2/n2)) ]

    tpoints <- data.table(x = seq(-4, +4, length.out = 100))
    tpoints[, y := dt(x, df = 4) ]
    tpoints[, group := 'T2D-Control']
    q <- ggplot(tpoints, aes(x = 2+x, y = y)) + theme_bw() + geom_line() + 
         xlab(expression(frac(Glc[T2D] - Glc[Control], frac(s, sqrt(n))))) + 
         ylab(expression(Density))

    gridExtra::grid.arrange(p, q, nrow = 2)
    j <- 1
    pp <- p + geom_point(   data = pointdt[i==j]) + geom_line(data = pointdt[i==j])
    pp <- pointdt[i==j & group == 'Control'][, pp + annotate('point', x = mean(x), y = mean(y), color = '#F8766D', shape = 108, size = 10) ]
    pp <- pointdt[i==j & group == 'T2D'    ][, pp + annotate('point', x = mean(x), y = mean(y), color = '#00BFC4', shape = 108, size = 10) ]
    pp + coord_flip()
    
    gridExtra::grid.arrange(pp, q, nrow = 2)
    pp <- pp + geom_crossbar(data = meandt[i==j], aes(x = x, xmin = xmin, xmax = xmax), width = 0.02, fill = 'white') + geom_point( data = pointdt[i==j])
    gridExtra::grid.arrange(pp, q, nrow = 2)
    pp + annotate(x = )

