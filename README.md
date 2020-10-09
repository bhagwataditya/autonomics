# Intro

Make omics data analysis flow :-).


# Install

    install.packages('remotes')
    remotes::install_github('bhagwataditya/importomics',       repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)

## Read and prepare

    # METABOLON
        require(importomics)
        file <- download_data("hypoglycemia.metabolon.xlsx")
        read_metabolon(file)
        
          require(magrittr)
          (object <- 'extdata/glutaminase/glutaminase.xlsx'    %>% 
                      system.file(package = 'autonomics.data') %>% 
                      autonomics::read_metabolon())
          object %<>% autonomics::prepare_metabolon()
    
    # SOMASCAN
          (object <- 'extdata/stemcomp/soma/stemcomp.adat'     %>% 
                      system.file(package = 'autonomics.data') %>% 
                      autonomics::read_somascan())
          object %<>% autonomics::prepare_somascan()
    
    # RNASEQ COUNTS
          object <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_counts(fid_var = 'gene_id')
          (object %<>% autonomics::prepare_rnaseq())

    # RNASEQ BAMFILES
    
         # Download example BAM files
           url <- "https://bitbucket.org/graumannlab/billing.stemcells/downloads/stemcomp.bamfiles.zip"
           dir.create('~/.autonomics', showWarnings = FALSE)
           destfile <- "~/.autonomics/stemcomp.bamfiles.zip"
           download.file(url, destfile = destfile, )
           utils::unzip(destfile, exdir = '~/.autonomics')
           unlink(destfile)

         # Download GTF file
           gtffile <- autonomics::download_gtf('Homo sapiens', 95)

         # Read BAM files into SummarizedExperiment
           object <- autonomics::read_bam(bamdir    = "~/.autonomics/stemcomp.bamfiles",
                                          gtffile   = gtffile,
                                          ispaired  = TRUE)
         # Prepare for analysis
           object %<>% autonomics::prepare_rnaseq()

    
    # PROTEINGROUPS
          object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_proteingroups()
          object %<>% autonomics::prepare_proteingroups()
    
    # EXIQON
          object <-  autonomics::read_exiqon(myfile)
          object %<>% autonomics::prepare_exiqon()

   
    # ANY OMICS DATASET
          object <- 'extdata/glutaminase/glutaminase.xlsx'   %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_omics(
                        sheet      = 2,
                        fid_rows   = 11:401,    fid_cols   = 5,
                        sid_rows   = 3,         sid_cols   = 15:86,
                        expr_rows  = 11:401,    expr_cols  = 15:86,
                        fvar_rows  = 10,        fvar_cols  = 1:14,
                        svar_rows  = 1:10,      svar_cols  = 14,
                        fdata_rows = 11:401,    fdata_cols = 1:14,
                        sdata_rows = 1:10,      sdata_cols = 15:86,
                        transpose  = FALSE)
                        
## Explore

    # Sample densities
        object <- autonomics.data::glutaminase
        object %>% autonomics::plot_sample_densities(color = subgroup, facet = subgroup)
        object %>% autonomics::plot_sample_boxplots( fill  = subgroup, coord_flip = FALSE)
        object %>% autonomics::plot_sample_violins(fill  = subgroup, coord_flip = FALSE)

    # Principal Component Analysis
        object <- autonomics.data::glutaminase
        object %>% autonomics::plot_pca_samples()
        object %>% autonomics::plot_pca_features()
        object %>% autonomics::plot_pca_samples_and_features(n=4)
        
    #  Linear Discriminant Analysis
        object %>% autonomics::plot_lda_samples()
        object %>% autonomics::plot_lda_features()
        object %>% autonomics::plot_lda_samples_and_features()
        
    # Partial Least Squares Analysis
        object %>% autonomics::plot_pls_samples()
        object %>% autonomics::plot_pls_features()
        object %>% autonomics::plot_pls_samples_and_features()
   
    # Spectral Map Analysis
        object %>% autonomics::plot_sma_samples()          
        object %>% autonomics::plot_sma_features() 
        object %>% autonomics::plot_sma_samples_and_features()


## Contrast

    # Add limma
       object <- autonomics.data::glutaminase
       table(object$subgroup)
       ctrdefs <- c(uM05.h10 = 'uM05_h10 - Veh_h10', 
                    uM10.h10 = 'uM10_h10 - Veh_h10')
       autonomics::contrastdefs(object) <- ctrdefs
       object %<>% autonomics::add_limma()
    
    # Plot contrast features
       object %>% autonomics::plot_contrast_features(contrast = ctrdefs[1], n=2)
       object %>% autonomics::plot_contrast_features(contrast = ctrdefs[2], n=2)
    
    # Write to file
       (file <- tempdir() %>% file.path('/glutaminase_results.txt'))
       object %>% autonomics::write_features(file)
    
    # Plot volcano
       object %>% autonomics::plot_volcano()
    
    # Plot venns
       object %>% autonomics::plot_contrast_venns(euler = TRUE)
    
