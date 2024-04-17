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


#---------------------------------------------------------------------------
#
#   FUKUDA ZEBRAFISH DEVELOPMENT
#   https://www.embopress.org/doi/full/10.15252/embr.201949752
#
#----------------------------------------------------------------------------
  
# Read
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, contaminants = TRUE, reverse = TRUE)
    fdt(object)
    fdt(object)$protein    <- NULL # identical to feature_id
    fdt(object)$log2maxlfq <- NULL # not needed here
    fdt(object)$organism   <- NULL # DANRE everywhere
    fdt(object)$isoform    <- NULL # not needed now
    table(systematic_nas(object))
    object %<>% filter_exprs_replicated_in_some_subgroup()
    table(systematic_nas(object))
    
# Three types of detections. Impute systematic NAs.
    fdt(object)$detection <- 'none'
    fdt(object)$detection[ systematic_nas(object) ] <- 'systematic'
    fdt(object)$detection[     random_nas(object) ] <- 'random'
    fdt(object)$detection[         no_nas(object) ] <- 'full'
    table(fdt(object)$detection)
    object %<>% impute()
    fdt(object)
    
# LinMod
    object %<>% fit_limma(coefs = 'Adult')
    object %<>% extract( order(fdt(.)$`p~Adult~limma`) , )
    fdt(object)
    
# Lets compose a representative set                                                                                                                 dir    det    dif
#                                                                                                                                                   ---    ---    ---
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` < 0 ][ detection == 'systematic' ]  #    54  mcm6        F1R5P3_DANRE      -8e-9         dn     sys     1    (1)
                                                                                               #  5575  polr2d      Q6DRG4_DANRE      -5e-1                        0    (2)
                                                                                               #  6266  riox2       Q7T3G6_DANRE      -2e-1                        0    (3)
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` < 0 ][ detection == 'full'       ]  #  5045  hmbsb       Q503D2_DANRE      -4e-8                full    1    (4)
                                                                                               #  2769  synm        B7ZUY8_DANRE      -1                           0    (5)
                                                                                               #  4750  pmpcb       Q1L8E2_DANRE      -1                           0    (6)
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` < 0 ][ detection == 'random'     ]  #  6198  hbbe2       Q7T1B0_DANRE      -1e-7                rnd     1    (7)
                                                                                               #  6020  atg3        ATG3_DANRE        -1                           0    (8)
                                                                                               #   641  lztfl1      A0A0R4ISB5_DANRE  -1                           0    (9)
    
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` > 0 ][ detection == 'systematic' ]  #  3568  vtg4        F1Q7L0_DANRE      +2e-8         up     sys     1   (10)
                                                                                               #  5065  ovca2       OVCA2_DANRE       +2e-1                        0   (11)
                                                                                               #  3335  mief2       E7FFR5_DANRE      +1e-1                        0   (12)
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` > 0 ][ detection == 'full'       ]  #  5477  hbba2       Q6DGK4_DANRE      +3e-6                full    1   (13)
                                                                                               #  3416  lrrc57      LRC57_DANRE       +1                           0   (14)
                                                                                               #   471  ryr2b       A0A0R4IKT8_DANRE  +1                           0   (15)
    fdt(object)[ contaminant == '' ][ `effect~Adult~limma` > 0 ][ detection == 'random'     ]  #  5522  bfb         Q6DHC4_DANRE      +1e-5                rnd     1   (16)
                                                                                               #   599  utp20       A0A0R4IQW7_DANRE  +1                           0   (17)
                                                                                               #  2711  plbd1       B3DJ80_DANRE      +1                           0   (18)
    
    fdt(object)[ contaminant == '+' ]                                                          #  2936  ALB         CON__ALBU_BOVIN   -8e-1         dn     full    0   (19)
    fdt(object)[     reverse == '+' ]                                                          #  6754  REV__elmo3  REV__F1QSV8_DANRE -3e-1         dn     full    0   (20)

    features <- c( 'F1R5P3_DANRE',  'Q6DRG4_DANRE', 'Q7T3G6_DANRE',  'Q503D2_DANRE', 'B7ZUY8_DANRE',     'Q1L8E2_DANRE',  'Q7T1B0_DANRE',       'ATG3_DANRE', 'A0A0R4ISB5_DANRE', 
                   'F1Q7L0_DANRE',   'OVCA2_DANRE', 'E7FFR5_DANRE',  'Q6DGK4_DANRE',  'LRC57_DANRE', 'A0A0R4IKT8_DANRE',  'Q6DHC4_DANRE', 'A0A0R4IQW7_DANRE',     'B3DJ80_DANRE' , 
                   'CON__ALBU_BOVIN', 'REV__F1QSV8_DANRE' )
    all( features %in% fdt(object)$feature_id )
    fdt(object)[ feature_id %in% features ][ order(as.numeric(proId)) ]

  read_maxquant_proteingroups('inst/extdata/fukuda20.proteingroups.txt') # works
  read_maxquant_proteingroups(download_data('fukuda20.proteingroups.txt')) # works
  