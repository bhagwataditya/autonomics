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
        subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
        rna <- read_rnaseq_counts(rnafile)
        pro <- read_maxquant_proteingroups(profile, subgroups = subgroups, contaminant = TRUE, reverse = TRUE)
        fos <- read_maxquant_phosphosites( fosfile = fosfile, profile = profile, subgroups = subgroups, contaminant = TRUE, reverse = TRUE )
        pro$subgroup %<>% split_extract_fixed('_', 1)
        fos$subgroup %<>% split_extract_fixed('_', 1)
        fdt(fos)$`Fasta headers` <- NULL
        fvars(rna) %<>% stri_replace_first_fixed('gene_name', 'gene')
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
    # Rough numbers       
        fdt(rna)[ , sum(`p~M00~limma` > 0.05) / .N * 100 ]                            # 46 % flat
        fdt(rna)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` < 0 ) / .N * 100 ]  # 26 % down
        fdt(rna)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` > 0 ) / .N * 100 ]  # 28 % up
        
        fdt(pro)[ , sum(`p~M00~limma` > 0.05) / .N * 100 ]                            # 25 % flat
        fdt(pro)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` < 0 ) / .N * 100 ]  # 52 % down
        fdt(pro)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` > 0 ) / .N * 100 ]  # 23 % up

        fdt(fos)[ , sum(`p~M00~limma` > 0.05) / .N * 100 ]                            # 17 % flat
        fdt(fos)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` < 0 ) / .N * 100 ]  # 53 % down
        fdt(fos)[ , sum(`p~M00~limma` < 0.05 & `effect~M00~limma` > 0 ) / .N * 100 ]  # 30 % up
        
        length(rnagenes <- unique(fdt(rna)[, gene]))                                # 22 855
        length(progenes <- unique(fdt(pro)[contaminant=='' & reverse =='', gene]))  #  7 128
        length(fosgenes <- unique(fdt(fos)[contaminant=='' & reverse =='', gene]))  #  2 551
        length(rnagenes %>% intersect(progenes)) #  6903 : 30% rnagenes, 97% proteingenes
        length(rnagenes %>% setdiff(progenes))   # 15952 : 70% rnagenes
        length(progenes %>% setdiff(rnagenes))   #   225 :  3% progenes
        
        length(progenes %>% intersect(fosgenes)) # 2158 : 30% progenes  81% fosgenes
        length(progenes %>% setdiff(  fosgenes)) # 4970 : 70% progenes not fosquantified
        length(fosgenes %>% setdiff(  progenes)) #  393 : 15% fosgenes not proquantified
    # Summarize
        sumRNA <- fdt(rna)[, .(
                                    RNA.down = sum(  `M00~limma` == 'down' ), 
                                    RNA.flat = sum(`p~M00~limma` >   0.05  ),
                                    RNA.up   = sum(  `M00~limma` ==  'up'  )
                              ), by = 'gene']
        sumPRO <- fdt(pro)[, .(
                                    PRO.down = sum(  `M00~limma` == 'down' ), 
                                    PRO.flat = sum(`p~M00~limma` >   0.05  ),
                                    PRO.up   = sum(  `M00~limma` ==  'up'  ),
                                    PRO.dimp = sum(  `M00~limma` == 'down' &   imputed ), 
                                    PRO.fimp = sum(`p~M00~limma` >   0.05  &   imputed ),
                                    PRO.ump  = sum(  `M00~limma` ==  'up'  &   imputed )
                              ), by = 'gene']
        sumFOS <- fdt(fos)[, .(
                                    FOS.down = sum(  `M00~limma` == 'down' ), 
                                    FOS.flat = sum(`p~M00~limma` >   0.05  ),
                                    FOS.up   = sum(  `M00~limma` ==  'up'  ),
                                    FOS.dimp = sum(  `M00~limma` == 'down' &   imputed ), 
                                    FOS.fimp = sum(`p~M00~limma` >   0.05  &   imputed ),
                                    FOS.ump  = sum(  `M00~limma` ==  'up'  &   imputed )
                              ), by = 'gene']
        sumdt <- sumRNA
        sumdt %<>% merge.data.table(sumPRO, by = 'gene', all.x = TRUE, all.y = TRUE)
        sumdt %<>% merge.data.table(sumFOS, by = 'gene', all.x = TRUE, all.y = TRUE)
        summat <- dt2mat(sumdt)
        summat[1:5, 1:5]
        summat[is.na(summat)] <- 0
        sumdt <- mat2dt(summat, 'gene')
        
        contaminantgenes <- union( unique( fdt(pro)[contaminant == '+', gene] ) ,
                                   unique( fdt(fos)[contaminant == '+', gene] ) )
        sumdt[, contaminant := '']
        sumdt %<>% pull_columns(c('gene', 'contaminant'))
        sumdt[contaminantgenes, on = 'gene', contaminant := '+']
        sumdt[contaminant == '+']
        sumdt
        
        sumdt[, reverse := '']
        sumdt %<>% pull_columns(c('gene', 'contaminant', 'reverse'))
        sumdt[stri_detect_fixed(gene, 'REV__'), reverse := '+']
        sumdt[reverse == '+']
        sumdt
        
    # Contaminant
        procontaminants <- fdt(pro)[contaminant == '+' & gene != '', unique(gene) ]  # Lets take a contaminant in both PRO and FOS
        foscontaminants <- fdt(fos)[contaminant == '+' & gene != '', unique(gene) ]  # That makes things a bit easier
        commoncontaminants <- intersect( procontaminants, foscontaminants )          # KRT19 catches the eye. Well-known contaminant
        sumdt[commoncontaminants, on = 'gene'][, c(1, 4:6 )]                         # RNA : 1 flat - contributes to flat background - take
        sumdt[commoncontaminants, on = 'gene'][, c(1, 7:12)]                         # PRO : 1 down - ok, gets filtered out anyways  - take
        sumdt[commoncontaminants, on = 'gene'][, c(1,13:18)]                         # FOS : 4 down (imputed), 1 flat (detected) - take the detected
        fdt(rna)[ 'KRT19', on = 'gene' ]
        fdt(rna)[gene %in% 'KRT19'] # f  KRT19
        fdt(pro)[gene %in% 'KRT19'] # d  KRT19
        fdt(fos)[gene %in% 'KRT19'] # f  KRT19 3812
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
    # RNA flat - PRO down - FOS none
        unique(fdt(rna)[  `p~M00~limma`>0.5                          ]$gene)   %>% intersect(  
        unique(fdt(pro)[`fdr~M00~limma`<5e-12 & `effect~M00~limma`<0 ]$gene) ) %>% setdiff(
        unique(fdt(fos)$gene))
        fdt(rna)[gene == 'CA13']
        fdt(pro)[gene == 'CA13']
        fdt(fos)[gene == 'CA13']
    # RNA flat - PRO down - FOS none
        unique(fdt(rna)[  `p~M00~limma`>0.5                          ]$gene)   %>% intersect(  
        unique(fdt(pro)[`fdr~M00~limma`<5e-9  & `effect~M00~limma`>0 ]$gene) ) %>% setdiff(
        unique(fdt(fos)$gene))
        fdt(rna)[gene == 'OGFOD3']
        fdt(pro)[gene == 'OGFOD3']
        fdt(fos)[gene == 'OGFOD3']
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
    # RNA down, PRO flat, FOS flat
        candgenes <-             unique(fdt(rna)[ `fdr~M00~limma` < 5e-5 & `effect~M00~limma` < 0 ]$gene)
        candgenes %<>% intersect(unique(fdt(pro)[imputed==FALSE][   `p~M00~limma` > 0.1                           ]$gene))
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][   `p~M00~limma` > 0.1                           ]$gene))
        candgenes
        fdt(rna)['SGK223', on = 'gene']                  # down
        fdt(pro)['SGK223', on = 'gene']                  # flat. but proteingroups file shows 3 isoforms
        fdt(fos)['SGK223', on = 'gene'][imputed==FALSE]  # lets take the two flats
    # RNA up, PRO flat, FOS flat
        candgenes <-             unique(fdt(rna)[ `fdr~M00~limma` < 5e-6 & `effect~M00~limma` > 0 ]$gene)
        candgenes %<>% intersect(unique(fdt(pro)[imputed==FALSE][   `p~M00~limma` > 0.5                           ]$gene))
        candgenes %<>% intersect(unique(fdt(fos)[imputed==FALSE][   `p~M00~limma` > 0.5                           ]$gene))
        candgenes
        fdt(rna)[gene == 'TGOLN2']
        fdt(pro)[gene == 'TGOLN2']
        fdt(fos)[gene == 'TGOLN2']
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
        
    # RNA down, PRO down, FOS down
                       fdt(rna)[ `fdr~M00~limma` < 5e-7 & `effect~M00~limma` < 0, unique(gene)]   %>% 
            intersect( fdt(pro)[ `fdr~M00~limma` < 5e-7 & `effect~M00~limma` < 0, unique(gene)] ) %>% 
            intersect( fdt(fos)[   `p~M00~limma` > 0.50,                          unique(gene)] )
        fdt(rna)[gene == 'ADD2']
        fdt(pro)[gene == 'ADD2'][imputed==FALSE]
    # RNA up, PRO up, FOS up
                       fdt(rna)[ `fdr~M00~limma` < 5e-6 & `effect~M00~limma` > 0, unique(gene)]   %>% 
            intersect( fdt(pro)[ `fdr~M00~limma` < 5e-6 & `effect~M00~limma` > 0, unique(gene)] ) %>% 
            intersect( fdt(fos)[   `p~M00~limma` > 0.50,                          unique(gene)] )
        fdt(rna)[gene == 'MAP1A']
        fdt(pro)[gene == 'MAP1A'][imputed==FALSE]

    # Flatliner: RNA and FOS (not in PRO)
        candidates <- rnagenes %>% setdiff(progenes) %>% intersect(fosgenes)
        candidates %<>% intersect(fdt(rna)[`p~M00~limma` > 0.30, unique(gene)])  # flat rna - flat pro
        candidates %<>% intersect(fdt(fos)[`p~M00~limma` > 0.30, unique(gene)])  # flat rna - flat pro
        fdt(rna)['EAF1', on = 'gene']
        fdt(pro)['EAF1', on = 'gene']
        fdt(fos)['EAF1', on = 'gene']
    # Flatliner: RNA (not in PRO and FOS)
        candidates <- rnagenes %>% setdiff(progenes) %>% setdiff(fosgenes)
        candidates %<>% intersect(fdt(rna)[`p~M00~limma` > 0.50, unique(gene)])  # flat rna - flat pro
        candidates
        fdt(rna)['GSC2', on = 'gene']

    