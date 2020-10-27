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
    file <- download_data('stemcells.bam.zip')                           # RNAseq BAM
    read_bam(file)
    
    file <- download_data('stemcells.rnacounts.txt')                     # RNAseq counts
    read_counts(file)
        
    file <- download_data('glutaminase.metabolon.xlsx')                  # Metabolon xlsx
    read_metabolon(file)
    
    proteingroups <- download_data('differentiation.proteinGroups.txt')  # LCMSMS Proteingroups
    read_proteingroups(proteingroups)
    
    phosphosites <- download_data('differentiation.phosphoSites.txt')    # LCMSMS Phosphosites
    read_phosphosites(phosphosites, proteingroupss)
    
    file <- download_data('hypoglycemia.somascan.adat')                  # SOMAscan adat
    read_somascan(file)
```

