---
title: importomics
subtitle: generifying and intuifying cross-platform omics analysis
---

![ ](https://bitbucket.org/graumannlabtools/importomics/downloads/read_prepro_analyze.png)


```
# Install
    install.packages('remotes')
    remotes::install_github('bhagwataditya/importomics',  
                             repos = BiocManager::repositories(), 
                             dependencies = TRUE, upgrade = FALSE)
                             
# Run
    file <- download_data('billing16.bam.zip')                     # RNAseq BAM
    read_bam(file, paired=TRUE, genome = 'hg38')
    
    file <- download_data('billing16.rnacounts.txt')               # RNAseq counts
    read_counts(file)
        
    file <- download_data('halama18.metabolon.xlsx')               # Metabolon xlsx
    read_metabolon(file)
    
    rm_subgroups <- c('BLANK_BM00', 'BM00_BM00', 'EM01_EM00', 'EM05_EM02', 'EM30_EM15')
    proteingroups <- download_data('billing19.proteingroups.txt')  # LCMSMS Proteingroups
    read_proteingroups(proteingroups, rm_subgroups = rm_subgroups)
    
    phosphosites <- download_data('billing19.phosphosites.txt')    # LCMSMS Phosphosites
    read_phosphosites(phosphosites, proteingroups, rm_subgroups = rm_subgroups)
    
    file <- download_data('atkin18.somascan.adat')                 # SOMAscan adat
    read_somascan(file)
```

