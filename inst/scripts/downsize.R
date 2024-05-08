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
  
  
#---------------------------------------------------------------------------
#
#   BILLING STEM CELL COMPARISON
#   https://www.nature.com/articles/srep21507
#
#----------------------------------------------------------------------------

    file <- download_data('billing16.proteingroups.txt')
    object <- read_maxquant_proteingroups(file)
    sdt(object)
  

    
#---------------------------------------------------------------------------------
#
#   BILLING STEM CELL DIFFERENTIATION
#   https://www.sciencedirect.com/science/article/pii/S1535947620315243?via%3Dihub
#
#---------------------------------------------------------------------------------
    # Read
        rnafile <- download_data('billing19.rnacounts.txt')
        profile <- download_data('billing19.proteingroups.txt')
        fosfile <- download_data('billing19.phosphosites.txt')
        rnafile1 <- 'inst/extdata/billing19.rnacounts.txt'
        profile1 <- 'inst/extdata/billing19.proteingroups.txt'
        fosfile1 <- 'inst/extdata/billing19.phosphosites.txt'
        subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
        rna  <- read_rnaseq_counts(rnafile)
        rna1 <- read_rnaseq_counts(rnafile1)
        pro  <- read_maxquant_proteingroups(profile,  subgroups = subgroups, contaminant = TRUE, reverse = TRUE)
        pro1 <- read_maxquant_proteingroups(profile1, subgroups = subgroups)
        fos  <- read_maxquant_phosphosites( fosfile = fosfile,  profile = profile,  subgroups = subgroups, contaminant = TRUE, reverse = TRUE )
        fos1 <- read_maxquant_phosphosites( fosfile = fosfile1, profile = profile1, subgroups = subgroups)
        pro$subgroup  %<>% split_extract_fixed('_', 1)
        pro1$subgroup %<>% split_extract_fixed('_', 1)
        fos$subgroup  %<>% split_extract_fixed('_', 1)
        fos1$subgroup %<>% split_extract_fixed('_', 1)
        fdt(fos)$`Fasta headers` <- NULL
        fvars(rna) %<>% stri_replace_first_fixed('gene_name', 'gene')
    # Pca
        pcalist <- list( p1 = biplot(pca(rna ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('RNA'),   #  COL3A1 and COL11A1
                         p3 = biplot(pca(pro ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('PRO'),   #   HSPB6 and LCP1
                         p5 = biplot(pca(fos ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('FOS') )  #     VIM and NES
        gridExtra::grid.arrange(grobs = pcalist, nrow = 3)
    # Model
        rna  %>% impute() # no NA
        pro %<>% impute()
        fos %<>% impute()
        rna %<>% fit_limma(coefs = 'M00') # differentiation E00 -> M00
        pro %<>% fit_limma(coefs = 'M00')
        fos %<>% fit_limma(coefs = 'M00')
        rna %<>% extract(order(fdt(.)$`p~M00~limma`), )
        pro %<>% extract(order(fdt(.)$`p~M00~limma`), )
        fos %<>% extract(order(fdt(.)$`p~M00~limma`), )
    # Rm missing values
        fdt(rna)[is.na(`p~M00~limma`)]  #  0 % NA
        fdt(pro)[is.na(`p~M00~limma`)]  #  7 % NA :  586/8017
        fdt(fos)[is.na(`p~M00~limma`)]  # 27 % NA : 2314/8590
        pro %<>% filter_features(!is.na(`p~M00~limma`))
        fos %<>% filter_features(!is.na(`p~M00~limma`))
        fdt(rna) %<>% add_adjusted_pvalues('fdr')
        fdt(pro) %<>% add_adjusted_pvalues('fdr')
        fdt(fos) %<>% add_adjusted_pvalues('fdr')
        rna %<>% abstract_fit()
        pro %<>% abstract_fit()
        fos %<>% abstract_fit()
    # Drops features with no gene annotation
        fdt(rna)[gene=='' | is.na(gene) | gene == 'NA'] # 0
        fdt(pro)[gene=='' | is.na(gene) | gene == 'NA'] # 108 - incomplete fastahdrs - rm
        fdt(fos)[gene=='' | is.na(gene) | gene == 'NA'] #  71 - incomplete fastahdrs - rm
        pro %<>% filter_features( gene!='' & !is.na(gene) & gene != 'NA' )
        fos %<>% filter_features( gene!='' & !is.na(gene) & gene != 'NA' )
    # Fix Excel genes
        fdt(rna)[order(gene)]  # excel genes issue
        fdt(pro)[order(gene)]  # no issue
        fdt(fos)[order(gene)]  # no issue
        fdt(rna)$gene %<>% fix_xlgenes()
        fdt(pro)$gene %<>% fix_xlgenes()
        fdt(fos)$gene %<>% fix_xlgenes()
        fdt(rna)[order(gene)]  # excel genes issue
    # Contaminant
        procontaminants <- fdt(pro)[contaminant == '+' & gene != '', unique(gene) ]  # Lets take a contaminant in both PRO and FOS
        foscontaminants <- fdt(fos)[contaminant == '+' & gene != '', unique(gene) ]  # That makes things a bit easier
        commoncontaminants <- intersect( procontaminants, foscontaminants )          # KRT19 catches the eye. Well-known contaminant
        fdt(pro)[commoncontaminants, on = 'gene']   # PRO : 1 down - ok, gets filtered out anyways  - take
        fdt(rna)[gene == 'KRT19', on = 'gene']   # RNA : 1 flat - contributes to flat background - take
        fdt(fos)[gene == 'KRT19', on = 'gene']   # FOS : 4 down (imputed), 1 flat (detected) - take the detected
        contaminants <- union(procontaminants, foscontaminants)
        rna %<>% filter_features(! gene %in% c('', contaminants))
        pro %<>% filter_features(! gene %in% c('', contaminants))
        fos %<>% filter_features(! gene %in% c('', contaminants))
    # Reverse
        proreverse <- fdt(pro)[reverse == '+', unique(gene)]  # Lets find one in both PRO and FOS
        fosreverse <- fdt(fos)[reverse == '+', unique(gene)]  # REV__SYNE2 is the only such protein!
        intersect(proreverse, fosreverse)                     # PRO : 1 flat - take
        fdt(pro)[gene == 'REV__SYNE2']                        # FOS : 2 up (imp) - same protein - take first (12087)
        fdt(fos)[gene == 'REV__SYNE2']                        # Interestingly the fwd protein is also detected!
        fdt(rna)[gene == 'SYNE2']                             # RNA : 1 down - take
        fdt(pro)[gene == 'SYNE2']                             # PRO : 1 down - take
        fdt(fos)[gene == 'SYNE2']                             # FOS : 1 down (imputed)
        pro %<>% filter_features(reverse=='')                 #         will not contribute to modeling background because NA
        fos %<>% filter_features(reverse=='')                 #         still good to take because of the pairing with PRO
    # Multiprotein genes
        multiproteindt <- fdt(pro)[imputed==FALSE]
        multiproteindt %<>% extract(, .SD[.N>1], by = 'gene')                                                    # multiprotein gene
        multiproteindt %<>% extract(gene %in% unique(fdt(rna)[`p~M00~limma`>0.1]$gene))                          # RNA flat
        multiproteindt %<>% extract(, .SD[ sum(   `fdr~M00~limma`<0.05 & `effect~M00~limma`>0)>0], by = 'gene')  # PRO up (SYNE2 was already down)
        multiproteindt %<>% extract(, .SD[ sum(     `p~M00~limma`>0.5 )>0], by = 'gene')                         # PRO flat
        fdt(rna)[gene == 'ACOT9']  # RNA: 1 flat
        fdt(pro)[gene == 'ACOT9']  # PRO: 1 up, 1 flat, genenames intuitive
        fdt(fos)[gene == 'ACOT9']  # FOS: none
        multiproteingenes <- fdt(pro)[, .SD[.N>1], by = 'gene']$gene
        rna %<>% filter_features(!gene %in% multiproteingenes)
        pro %<>% filter_features(!gene %in% multiproteingenes)
        fos %<>% filter_features(!gene %in% multiproteingenes)
    # RNA down - PRO flat - FOS none
        unique(fdt(pro)[  `p~M00~limma`>0.5                         ]$gene)   %>% intersect(  
        unique(fdt(rna)[`fdr~M00~limma`<5e-8 & `effect~M00~limma`<0 ]$gene) ) %>% setdiff(
        unique(fdt(fos)$gene))
        fdt(rna)[gene == 'KIF5A']
        fdt(pro)[gene == 'KIF5A']
        fdt(fos)[gene == 'KIF5A']
    # RNA up - PRO flat - FOS none
        unique(fdt(pro)[  `p~M00~limma`>0.5                         ]$gene)   %>% intersect(  
        unique(fdt(rna)[`fdr~M00~limma`<5e-9 & `effect~M00~limma`>0 ]$gene) ) %>% setdiff(
        unique(fdt(fos)$gene))
        fdt(rna)[gene == 'CHD3']
        fdt(pro)[gene == 'CHD3']
        fdt(fos)[gene == 'CHD3']
    # RNA flat, PRO flat, FOS down+up (and some flat)
        candgenes <-             unique(fdt(rna)[  `p~M00~limma` > 0.5 ]$gene)                              # RNA flat
        candgenes %<>% intersect(unique(fdt(pro)[  `p~M00~limma` > 0.5 ]$gene))                             # PRO flat
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][`fdr~M00~limma` < 0.05 & `effect~M00~limma` < 0 ]$gene))   # FOS down
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][`fdr~M00~limma` < 0.05 & `effect~M00~limma` > 0 ]$gene))   # FOS up
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][  `p~M00~limma` > 0.05                          ]$gene))   # FOS flat
        fdt(rna)['TNKS1BP1', on = 'gene']
        fdt(pro)['TNKS1BP1', on = 'gene']
        fdt(fos)['TNKS1BP1', on = 'gene'][imputed==FALSE]
    # RNA flat, PRO flat, FOS down (and some flat)
        candgenes <-             unique(fdt(rna)[  `p~M00~limma` > 0.5 ]$gene)                              # RNA flat
        candgenes %<>% intersect(unique(fdt(pro)[  `p~M00~limma` > 0.5 ]$gene))                             # PRO flat
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][`fdr~M00~limma` < 5e-5 & `effect~M00~limma` < 0 ]$gene))   # FOS down
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][  `p~M00~limma` > 0.05                          ]$gene))   # FOS flat
        candgenes
        fdt(rna)[gene == candgenes[2]]
        fdt(rna)['KLC4', on = 'gene']  # flat
        fdt(pro)['KLC4', on = 'gene']  # flat
        fdt(fos)['KLC4', on = 'gene'][imputed==FALSE]  # one down, (one mildly up), one flat
    # RNA up, PRO flat, FOS up (and some flat)
        candgenes <-             unique(fdt(rna)[  `p~M00~limma` > 0.5 ]$gene)                              # RNA flat
        candgenes %<>% intersect(unique(fdt(pro)[  `p~M00~limma` > 0.5 ]$gene))                             # PRO flat
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][`fdr~M00~limma` < 5e-4 & `effect~M00~limma` > 0 ]$gene))   # FOS down
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][  `p~M00~limma` > 0.05                          ]$gene))   # FOS flat
        candgenes
        fdt(rna)['MFF', on = 'gene']                  # flat
        fdt(pro)['MFF', on = 'gene']                  # flat. but proteingroups file shows 3 isoforms
        fdt(fos)['MFF', on = 'gene'][imputed==FALSE]  # one up, one flat
        fdt(rna)['SLC4A7', on = 'gene']                  # flat
        fdt(pro)['SLC4A7', on = 'gene']                  # flat
        fdt(fos)['SLC4A7', on = 'gene'][imputed==FALSE]  # take 1 up and 2 flat
    # RNA flat, PRO impute down, FOS impute down
        candgenes <-      unique(fdt(rna)[ `p~M00~limma` > 0.5 ]$gene)                             # RNA flat
        candgenes %<>% intersect(fdt(pro)[imputed==TRUE]$gene)                                     # PRO imputed
        candgenes %<>% intersect(fdt(fos)[imputed==TRUE]$gene)                                     # FOS imputed
        candgenes %<>% intersect(fdt(pro)[ `fdr~M00~limma` < 5e-6 & `effect~M00~limma` < 0 ]$gene) # PRO down
        candgenes %<>% intersect(fdt(fos)[ `fdr~M00~limma` < 5e-6 & `effect~M00~limma` < 0 ]$gene) # FOS down
        candgenes
        fdt(rna)['MNAT1', on = 'gene']
        fdt(pro)['MNAT1', on = 'gene']
        fdt(fos)['MNAT1', on = 'gene']
    # RNA flat, PRO impute up, FOS impute up
        candgenes <-      unique(fdt(rna)[ `p~M00~limma` > 0.5 ]$gene)                             # RNA flat
        candgenes %<>% intersect(fdt(pro)[imputed==TRUE]$gene)                                     # PRO imputed
        candgenes %<>% intersect(fdt(fos)[imputed==TRUE]$gene)                                     # FOS imputed
        candgenes %<>% intersect(fdt(pro)[ `fdr~M00~limma` < 5e-3 & `effect~M00~limma` > 0 ]$gene) # PRO up
        candgenes %<>% intersect(fdt(fos)[ `fdr~M00~limma` < 5e-3 & `effect~M00~limma` > 0 ]$gene) # FOS up
        candgenes
        fdt(rna)['FYN', on = 'gene']
        fdt(pro)['FYN', on = 'gene']
        fdt(fos)['FYN', on = 'gene']
    # RNA flat, PRO flat, FOS flat, random NA
        candgenes <-             unique(fdt(rna)[ `p~M00~limma` > 0.2 ]$gene)   # RNA flat
        candgenes %<>% intersect(unique(fdt(pro)[ `p~M00~limma` > 0.2 ]$gene))  # PRO flat
        candgenes %<>% intersect(unique(fdt(fos)[ `p~M00~limma` > 0.2 ]$gene))  # FOS flat
        candgenes %<>% intersect(fdt(pro)$gene[random_nas(pro)])
        candgenes %<>% intersect(fdt(fos)$gene[random_nas(fos)])
        candgenes
        fdt(rna)['MFF', on = 'gene']
        fdt(pro)['MFF', on = 'gene'] # multiple isoforms
        fdt(fos)['MFF', on = 'gene']
        fdt(rna)['DDA1', on = 'gene']
        fdt(pro)['DDA1', on = 'gene']
        fdt(fos)['DDA1', on = 'gene']
        fdt(rna)['NEK1', on = 'gene']
        fdt(pro)['NEK1', on = 'gene']
        fdt(fos)['NEK1', on = 'gene']
        
    # Flatliner: RNA and FOS (not in PRO)
        candidates <- rnagenes %>% setdiff(progenes) %>% intersect(fosgenes)
        candidates %<>% intersect(fdt(rna)[`p~M00~limma` > 0.30, unique(gene)])  # flat rna - flat pro
        candidates %<>% intersect(fdt(fos)[`p~M00~limma` > 0.30, unique(gene)])  # flat rna - flat pro
        fdt(rna)['EAF1', on = 'gene']
        fdt(pro)['EAF1', on = 'gene']
        fdt(fos)['EAF1', on = 'gene']
    # Review
        rnafile <- download_data('billing19.rnacounts.txt')
        profile <- download_data('billing19.proteingroups.txt')
        fosfile <- download_data('billing19.phosphosites.txt')
        rnafile1 <- 'inst/extdata/billing19.rnacounts.txt'
        profile1 <- 'inst/extdata/billing19.proteingroups.txt'
        fosfile1 <- 'inst/extdata/billing19.phosphosites.txt'
        subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
        rna  <- read_rnaseq_counts(rnafile)
        rna1 <- read_rnaseq_counts(rnafile1)
        pro  <- read_maxquant_proteingroups(profile,  subgroups = subgroups, contaminant = TRUE, reverse = TRUE)
        pro1 <- read_maxquant_proteingroups(profile1, subgroups = subgroups)
        fos  <- read_maxquant_phosphosites( fosfile = fosfile,  profile = profile,  subgroups = subgroups, contaminant = TRUE, reverse = TRUE )
        fos1 <- read_maxquant_phosphosites( fosfile = fosfile1, profile = profile1, subgroups = subgroups)
        pro$subgroup  %<>% split_extract_fixed('_', 1)
        pro1$subgroup %<>% split_extract_fixed('_', 1)
        fos$subgroup  %<>% split_extract_fixed('_', 1)
        fos1$subgroup %<>% split_extract_fixed('_', 1)
        fdt(fos)$`Fasta headers` <- NULL
        fvars(rna) %<>% stri_replace_first_fixed('gene_name', 'gene')
        pcalist <- list( p1 = biplot(pca(rna ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('RNA'),  #  COL3A1
                         p2 = biplot(pca(rna1), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('rna'),  #  COL11A1
                         p3 = biplot(pca(pro ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('PRO'),  #  HSPB6
                         p4 = biplot(pca(pro1), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('pro'),  #  LCP1
                         p5 = biplot(pca(fos ), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('FOS'),  #  VIM
                         p6 = biplot(pca(fos1), nx = 1, ny = 1, feature_label = 'gene') + ggtitle('fos') ) #  NES
        gridExtra::grid.arrange(grobs = pcalist, layout_matrix = matrix(1:6, nrow = 3, byrow = TRUE))
        rna  %<>% fit_limma()
        pro  %<>% fit_limma()
        fos  %<>% fit_limma()
        rna1 %<>% fit_limma()
        pro1 %<>% fit_limma()
        fos1 %<>% fit_limma()
        gridExtra::grid.arrange( plot_volcano(rna,  label = 'gene_name'), 
                                 plot_volcano(rna1, label = 'gene'     ) )
        gridExtra::grid.arrange( plot_volcano(pro,  label = 'gene'     ), 
                                 plot_volcano(pro1, label = 'gene'     ) )
        gridExtra::grid.arrange( plot_volcano(fos,  label = 'gene'     ), 
                                 plot_volcano(fos1, label = 'gene'     ) )
        
        plotlist <- list( p1 = plot_volcano(rna,  label = 'gene_name' ),
                          p2 = plot_volcano(rna1, label = 'gene'),
                          p3 = plot_volcano(pro,  label = 'gene'),
                          p4 = plot_volcano(pro1, label = 'gene'),
                          p5 = plot_volcano(fos,  label = 'gene'),
                          p6 = plot_volcano(fos1, label = 'gene') )
        gridExtra::grid.arrange(grobs = plotlist, layout_matrix = matrix(1:6, nrow = 3, byrow = TRUE))


        
#---------------------------------------------------------------------------------
#
#   BILLING STEM CELL COMPARISON
#   https://www.nature.com/articles/srep21507
#
#---------------------------------------------------------------------------------
        # RNA, PRO, and SOMA exist
        # But the examples us only PRO
        # So for now, lets just downsize that.
        profile <-  download_data('billing16.proteingroups.txt')
        somafile <- download_data('billing16.somascan.adat')
        rna <- read_rnaseq_counts(rnafile)
        pro <- read_maxquant_proteingroups(profile)
        soma <- read_somascan(somafile)
        rna # already a reduced set file !
        pro
        soma
        
#---------------------------------------------------------------------------------
#
#   UNIPROT_HSA_20140515
#
#---------------------------------------------------------------------------------
        profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
        fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics')
        prodt <- fread(profile, select = c('id', 'Majority protein IDs', 'Gene names'))
        fosdt <- fread(fosfile, select = c('id', 'Protein group IDs', 'Proteins', 'Gene names'))
        
        uniprots <-         unique(unlist(stri_split_fixed(prodt$`Majority protein IDs`, ';')))
        uniprots %<>% union(unique(unlist(stri_split_fixed(fosdt$`Proteins`,             ';'))))
        uniprots %<>% extract(. != '')
        uniprots %<>% extract(!stri_detect_fixed(., 'REV__'))
        uniprots %<>% extract(!stri_detect_fixed(., 'CON__'))
      # uniprots %<>% split_extract_fixed('-', 1)
        uniprots %<>% unique()

        fastainfile <- download_data('uniprot_hsa_20140515.fasta')
        fastaobj <- Biostrings::readAAStringSet(fastainfile)
        fastaups <- names(fastaobj) %>% split_extract_fixed('|', 2)
        idx <- match(uniprots, fastaups)
        fastaobj %<>% extract(idx)
        fastaoutfile <- 'inst/extdata/uniprot_hsa_20140515.fasta'
        Biostrings::writeXStringSet(fastaobj, filepath = fastaoutfile)
        
        pro1 <- read_maxquant_proteingroups(profile)
        pro2 <- read_maxquant_proteingroups(profile, fastafile = fastafile)
        fdt(pro1)
        fdt(pro2)

        
#----------------
#
#   DIANN DATASET
#   https://www.nature.com/articles/s41467-023-40596-0
#
#----------------
        # PCA and PLS: major variation drivers
        file <- '../../../02_analysis/ag_serrano/ag_pogge-serrano-pkp2-001/atnomx_report_names.tsv'
        object <- read_diann_proteingroups(file)
        sdt(object)
        p1 <- biplot(pca(object), nx = 2, ny = 4, feature_label = 'gene')
        p2 <- biplot(pls(object), nx = 2, ny = 2, feature_label = 'gene') # maybe better these to raise less questions
        gridExtra::grid.arrange(p1, p2, nrow = 1)
        
        genes <- c(  
                    'PI16',        #  down in MA
                    'CPM',         #    up in MA
                    'RPL36AL',     #  down in KD (MA)
                    'ANKRD54',     #    up in KD (MA)
                    'TFB1M',       #  down in KD (MA + PA)
                    'NAV1'         #    up in KD (MA + PA)
                 )
        
        # Maybe lets look at the clusters
        object0 <- object
        fdt(object) %<>% extract(, 1:13)
        fdt(object)
        object %<>% fcluster()                      #    pamk : 2 major clusters
        fdt(object)
        table(fdt(object)$pamk)
        object %<>% fcluster(method = 'apcluster')  # apclust - if corrmat computation is bottleneck do onlc once
        

        # Flat backgrounders
        object
        object$sample_id %<>% factor( c(  sprintf('PA_ctrl_0%d', 1:4), 
                                          sprintf('PA_kd_0%d',   1:4), 
                                          sprintf('MA_ctrl_0%d', 1:4), 
                                          sprintf('MA_kd_0%d',   1:4)  ))
        plotdt <- sumexp_to_longdt(object, fvars = 'gene')
        ggplot(plotdt, aes(x = sample_id, y = value, group = feature_id, color = feature_id)) + 
        geom_point() + 
        geom_line() + 
        guides(color = 'none') + 
        theme_bw()
        
        subject <- filter_features(object, gene %in% genes)

        