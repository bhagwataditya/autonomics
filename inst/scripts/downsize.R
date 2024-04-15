#---------------------------------------------------------------------------
#
#   ATKIN HYPO
#   https://dom-pubs.pericles-prod.literatumonline.com/doi/10.1111/dom.13602
#
#----------------------------------------------------------------------------

mfile <- download_data('atkin.metabolon.xlsx')
sfile <- download_data('atkin.somascan.adat')
mobj <- read_metabolon(mfile)
sobj <- read_somascan(sfile)
sdt(mobj)
sdt(sobj)
mobj$Bmi <- NULL
sobj$Bmi <- NULL
mobj$Age <- NULL
sobj$Age <- NULL
sdt(mobj)
sdt(sobj)

PLOT_EXPRS <- function(obj)  plot_exprs(obj, block = 'Subject', coefs = NULL, shape = 'Diabetes', size = 'Diabetes') + 
                             scale_size_manual(values = c(Control = 1, T2DM = 2))

  mobj[ 'glucose',                                    ]  %>%  PLOT_EXPRS() # 1
  mobj[ 'cortisol',                                   ]  %>%  PLOT_EXPRS() # 2
  mobj[ 'docosatrienoate (22:3n3)',                   ]  %>%  PLOT_EXPRS() # 3
  mobj[ '10-nonadecenoate (19:1n9)',                  ]  %>%  PLOT_EXPRS() # 4 
# mobj[ 'X - 12450',                                  ]  %>%  PLOT_EXPRS() # 5
# mobj[ 'dihomo-linoleate (20:2n6)',                  ]  %>%  PLOT_EXPRS() # 6
# mobj[ 'linolenate [alpha or gamma; (18:3n3 or 6)]', ]  %>%  PLOT_EXPRS() # 7

  sobj['Activated Protein C',                         ]  %>%  PLOT_EXPRS() #  1
  sobj['BOC',                                         ]  %>%  PLOT_EXPRS() #  2
  sobj['CRDL1',                                       ]  %>%  PLOT_EXPRS() #  3
  sobj['Dynactin subunit 2',                          ]  %>%  PLOT_EXPRS() #  4
# sobj['Ephrin-A5',                                   ]  %>%  PLOT_EXPRS() #  5
# sobj['FABP',                                        ]  %>%  PLOT_EXPRS() #  6
# sobj['Insulin',                                     ]  %>%  PLOT_EXPRS() #  7
# sobj['Myoglobin',                                   ]  %>%  PLOT_EXPRS() #  8
# sobj['PH',                                          ]  %>%  PLOT_EXPRS() #  9
# sobj['PRL',                                         ]  %>%  PLOT_EXPRS() # 10
# sobj['RGMA',                                        ]  %>%  PLOT_EXPRS() # 11
# sobj['RSPO2',                                       ]  %>%  PLOT_EXPRS() # 12
# sobj['WIF-1',                                       ]  %>%  PLOT_EXPRS() # 13

  mfeatures <- c('glucose', 'cortisol', 'docosatrienoate (22:3n3)', '10-nonadecenoate (19:1n9)')
  sfeatures <- c('Activated Protein C', 'BOC', 'CRDL1', 'Dynactin subunit 2')
  
  mfeatures %<>% c(fnames(mobj)[order(rowVars(values(mobj), na.rm = TRUE))[1:16]])
  sfeatures %<>% c(fnames(sobj)[order(rowVars(values(sobj), na.rm = TRUE))[1:16]])

  mobj %<>% extract(mfeatures, )
  sobj %<>% extract(sfeatures, )
  
  biplot(pca(mobj))
  biplot(pca(sobj))
  
  plot_volcano(fit_limma(mobj), coefs = 't2')
  plot_volcano(fit_limma(sobj), coefs = 't2')

    
  read_metabolon('inst/extdata/atkin.metabolon.xlsx')
  read_somascan( 'inst/extdata/atkin.somascan.adat')

