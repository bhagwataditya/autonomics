
#' Log2 transform exprs
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return Updated SummarizedExperiment
#' @export
log2transform <- function(object, verbose = FALSE){
    if (verbose) message('\t\tLog2 transform')
    exprs(object) %<>% log2()
    return(object)
}


#===========================================================
# EXIQON
#===========================================================

#' Prepare exiqon sumexp for analysis
#' @param filter_features              string: fvar condition on which to filter features
#' @param align_sample_means           logical
#' @param lod                          number: Ct value beyond which to consider signal as NA
#' @param qrilc_consistent_nondetects  logical
#' @param filter_conserved_in_human    logical
#' @param plot                         logical
#' @param verbose                      logical
#' @examples
#' if (require(subramanian.2016)){
#'    require(magrittr)
#'    object <- 'extdata/exiqon/subramanian.2016.exiqon.xlsx' %>%
#'               system.file(package = 'subramanian.2016')    %>%
#'               read_exiqon()
#'    object %>% prepare_exiqon(lod = 36)
#' }
#' @importFrom magrittr %>%
#' @export
prepare_exiqon <- function(
   object,
   filter_features             = '`#RefGenes`==0 & `#Spike`   ==0',
   align_sample_means          = TRUE,
   lod                         = 37,
   qrilc_consistent_nondetects = TRUE,
   filter_conserved_in_human   = FALSE,
   plot                        = TRUE,
   verbose                     = TRUE
){

   # Check
   #------
   number <- NULL

   # Plot function
   #--------------
   plotfun <- function(object, xlab, title){
               object %>%
               autonomics.plot::plot_overlayed_sample_distributions() +
               ggplot2::xlab(xlab) +
               ggplot2::ggtitle(title)}

   # sdata
   #------
   object$subgroup <- object$sample_id %>% autonomics.import::guess_subgroup_values(verbose = verbose)

   # exprs
   #------
     # Original
       if (plot) object %>% plotfun('Ct', 'Load Ct values') %>% print()
     # Filter features
       object %<>% autonomics.import::filter_features_(filter_features, verbose = TRUE)
     # Align sample means
       if (align_sample_means){
         autonomics.support::cmessage('\t\tAlign sample means')
         sample_means <- object %>% autonomics.import::exprs() %>% (function(y){y[y>32] <- NA; y}) %>% colMeans(na.rm = TRUE)
         sample_diffs <- sample_means - stats::median(sample_means)
         autonomics.import::exprs(object) %<>% sweep(2, sample_diffs)
         if (plot) object %>% plotfun('Ct', 'Align sample means') %>% print()
       }
     # NA beyond lod
       object %<>% autonomics.preprocess::na_exprs_weakly_gt(lod)
       if (plot) object %>% plotfun('Ct', 'NA Ct >= lod') %>% print()
     # Invert scale
       autonomics.support::cmessage('\t\tInvert exprs = %d - Ct', lod)
       autonomics.import::exprs(object) %<>% magrittr::subtract(lod, .)
       newmetric <- sprintf('%d - Ct', lod)
       if (plot) object %>% plotfun(newmetric, newmetric) %>% print()
     # Impute consistent NAs
       if (qrilc_consistent_nondetects){
         object %<>% autonomics.import::impute_consistent_nas(verbose = TRUE)
         if (plot) object %>% plotfun(newmetric, 'Impute consistent NA values') %>% print()
       }
     # Filter for mirs conserved in human
       if (filter_conserved_in_human){
         autonomics.import::fdata(object) %<>% (function(x){x$number <- x$feature_id %>% substr(9, nchar(.)); x})
         conserved_human_mirs <- autonomics.annotate::load_targetscan('H.sapiens') %>% magrittr::extract2('number')
         object %<>% autonomics.import::filter_features(number %in% conserved_human_mirs, verbose = TRUE)
       }

   # Return
   #-------
   object

}


#============================================================
# SOMASCAN
#============================================================

#' Prepare somascan
#'
#' @param object                 SummarizedExperiment
#' @param filter_sample_type     string vector: sample  types to be filtered for.     Probably a subset of c('Sample', 'QC', 'Buffer', 'Calibrator').
#' @param filter_feature_type    string vector: feature types to be filtered for.     Probably a subset of c('Protein', 'Hybridization Control Elution', 'Rat Protein').
#' @param filter_sample_quality  string vector: sample  qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param filter_feature_quality string vector: feature qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param rm_na_svars            logical  : whether to rm NA svars
#' @param rm_single_value_svars  logical  : whether to rm single value svars
#' @param infer_design           logical  : whether to infer design from sampleids
#' @param log2transform          logical  : whether to log2 transform
#' @return Summarizedexperiment
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- 'extdata/stemcomp/soma/stemcomp.adat' %>%
#'               system.file(package = 'autonomics.data') %>%
#'               read_somascan()
#'    object %>% prepare_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
#' @importFrom magrittr %>%
#' @export
prepare_somascan <- function(
   object,
   filter_sample_type      = 'Sample',
   filter_feature_type     = 'Protein',
   filter_sample_quality   = c('FLAG', 'PASS'),
   filter_feature_quality  = c('FLAG', 'PASS'),
   rm_na_svars             = TRUE,
   rm_single_value_svars   = TRUE,
   log2transform           = TRUE
){
   # Comply
      SampleType <- RowCheck <- Type <- ColCheck <- NULL

   # Assert
      assertive.types::assert_is_character(c(filter_sample_type, filter_feature_type, filter_sample_quality, filter_feature_quality))
      assertive.types::assert_is_logical(c(rm_na_svars, rm_single_value_svars, log2transform))

   # Filter
      if ('SampleType' %in% autonomics.import::svars(object)){ # sample type - older versions don't have it
         message('\t\t=========================================================================')
         autonomics.support::cmessage_df('\t\t%s', table(`Sample types` = autonomics.import::sdata(object)$SampleType))
         object %<>% autonomics.import::filter_samples(SampleType %in% !!rlang::enquo(filter_sample_type), verbose = TRUE)
      }
      if ('RowCheck'   %in% autonomics.import::svars(object)){ # sample quality
         message('\t\t=========================================================================')
         autonomics.support::cmessage_df('\t\t%s', table(`Sample qualities ("RowCheck")` = autonomics.import::sdata(object)$RowCheck))
         object %<>% autonomics.import::filter_samples(RowCheck %in% !!rlang::enquo(filter_sample_quality), verbose = TRUE)
      }
      if ('Type'       %in% autonomics.import::fvars(object)){ # feature type
         message('\t\t=========================================================================')
         autonomics.support::cmessage_df('\t\t%s', table(`Type` = autonomics.import::fdata(object)$Type))
         object %<>% autonomics.import::filter_features(Type %in% !!rlang::enquo(filter_feature_type), verbose = TRUE)
      }
      if ('ColCheck'   %in% autonomics.import::fvars(object)){ # feature quality
         message('\t\t=========================================================================')
         autonomics.support::cmessage_df('\t\t%s', table(`Feature qualities ("ColCheck")` = autonomics.import::fdata(object)$ColCheck))
         object %<>% autonomics.import::filter_features(ColCheck %in% !!rlang::enquo(filter_feature_quality), verbose = TRUE)
         message('\t\t=========================================================================')
      }

   # Select
      if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
      if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Log2 transform
      if (log2transform) object %<>% autonomics.import::log2transform(verbose = TRUE)

   # Return
   object
}


#======================
# METABOLON
#======================

#' Prepare metabolon sumexp for analysis
#' @param object                        SummarizedExperiment
#' @param log2transform                 logical: whether to log2 transform
#' @param qirlc_consistent_nondetects   logical(1)
#' @param add_kegg_pathways             logical(1): whether to add KEGG pathways to fdata
#' @param add_smiles                    logical(1): whether to add SMILES to fdata
#' @return SummarizedExperiment
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- 'extdata/glutaminase/glutaminase.xlsx' %>%
#'               system.file(package = 'autonomics.data') %>%
#'               read_metabolon()
#'    object %>% prepare_metabolon()
#' }
#' @importFrom magrittr %>%
#' @export
prepare_metabolon <- function(
   object,
   log2transform          = TRUE,
   impute_consistent_nas  = FALSE,
   add_kegg_pathways      = FALSE,
   add_smiles             = FALSE
){
   if (log2transform)         object %<>% autonomics.import::log2transform(verbose = TRUE)
   if (impute_consistent_nas) object %<>% autonomics.import::impute_consistent_nas()
   if (add_kegg_pathways)     object %<>% autonomics::kegg_entry_to_pathways(entry_var = 'KEGG',    pathway_var = 'KEGGPATHWAY')
   if (add_smiles)            object %<>% autonomics::pubchem_to_smiles(   pubchem_var = 'PUBCHEM', smiles_var  = 'SMILES')
   object
}


#===============================================
# LCMSMS PROTEINGROUPS and PHOSPHOSITES
#===============================================

#' Deconvolute proteingroups
#' @param object             SummerizedExperiment with proteinGroups data
#' @param fasta              path to fastafile (if character) or output of
#' load_uniprot_fasta_annotations (autonomics.annotate)
#' @param fastafields        character vector: fields to load from fastafile
#' @param drop_isoform_info  logical: whether to drop isoform info
#' @param verbose            logical: whether to report back or not
#' @return deconvoluted and annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    fastafile <- '../data/uniprot_hsa_20140515.fasta'
#'    if (file.exists(fastafile)){
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt'            %>%
#'                  system.file(package='autonomics.data')                  %>%
#'                  read_proteingroups()
#'       fdata(object) %>% head()
#'
#'       object %<>% magrittr::extract(1:100, )                             %>%
#'                   deconvolute_proteingroups(fastafile = fastafile)
#'       fdata(object) %>% head()
#'    }
#' }
#' @importFrom magrittr %>% %<>%
#' @importFrom data.table data.table := .N
#' @export
deconvolute_proteingroups <- function(
   object,
   fasta,
   fastafields = c('GENES', 'PROTEIN-NAMES', 'EXISTENCE', 'REVIEWED'),
   drop_isoform_info = FALSE,
   verbose = TRUE
){
   # Satisfy CHECK
   EXISTENCE <- GENES <- IS.FRAGMENT <- ISOFORM <- N <- NGENE <- NISOFORMS <- NPERACCESSION <- NULL
   `PROTEIN-NAMES` <- REVIEWED <- `Uniprot accessions` <- ngene <- nprotein <- nseq <- .SD <- NULL

   # Load fasta annotations
   if (verbose) autonomics.support::cmessage('\t\tLoad fasta file')
   fasta_annotations <-  if (inherits(fasta, "fasta_data.table")) {
      fasta
   } else {
      fasta %>% autonomics.annotate::load_uniprot_fasta_annotations(fastafields)
   }

   # Uncollapse
   fdata1 <- object %>%
      autonomics.import::fdata() %>%
      magrittr::extract(, c('feature_id', 'Uniprot accessions'), drop = FALSE) %>%
      tidyr::separate_rows(`Uniprot accessions`, sep = ';') %>%
      data.table::data.table()

   # Split into CANONICAL and isoform
   fdata1 %<>% magrittr::extract(, ISOFORM     := `Uniprot accessions`) %>%
      magrittr::extract(, `Uniprot accessions` := `Uniprot accessions` %>% stringi::stri_replace_first_regex('[-][0-9]+',     '')) # %>%
   #magrittr::extract(, ISOFORM := ISOFORM %>% sort() %>% unique() %>% paste0(collapse=';'), by = c('feature_id', 'Uniprot accessions')) %>%
   #unique()

   # Merge in uniprot fasta annotations
   nunmapped <- fdata1 %>%
      magrittr::extract(ISOFORM %in% setdiff(ISOFORM, fasta_annotations$UNIPROTKB)) %>%
      magrittr::extract(, .SD[1], by = 'feature_id') %>%
      nrow()
   if (verbose) autonomics.support::cmessage('\t\tDeconvolute %d/%d proteingroups with sequences from fastafile',
                                             nrow(object) - nunmapped, nrow(object))
   fdata1 %<>% merge(fasta_annotations, by.x = 'Uniprot accessions', by.y = 'UNIPROTKB', sort = FALSE)
   report_n <- function(dt, prefix='', suffix=''){
      n <- dt %>% magrittr::extract(, .SD[, list(nseq = .N,
                                                 ngene = length(unique(GENES)),
                                                 nprotein = length(unique(`PROTEIN-NAMES`)))],
                                    by = 'feature_id') %>%
         magrittr::extract(, list(ngroups = .N,
                                  nsinglegene    = sum(ngene==1),
                                  nsingleprotein = sum(nprotein==1),
                                  nsingleseq     = sum(nseq==1)))
      autonomics.support::cmessage('\t\t\t%s%d proteingroups -> %d singlegene -> %d singleprotein -> %d singleseq%s',
                                                stringi::stri_pad_right(prefix, width = 60), n$ngroups,          n$nsinglegene,   n$nsingleprotein,   n$nsingleseq, suffix)
   }
   if (verbose) message('')
   if (verbose) fdata1 %>% report_n(prefix = 'All proteingroups')

   # Prefer best existence
   fdata1 %<>% magrittr::extract(, .SD[EXISTENCE == min(EXISTENCE)], by = 'feature_id')
   if (verbose) fdata1 %>% report_n(prefix = 'Per group: drop inferior existences')

   # Drop trembl entries from swissprot groups
   fdata1 %<>% magrittr::extract(, .SD[REVIEWED == max(REVIEWED)], by = 'feature_id')
   swissprot <- fdata1[REVIEWED==1]
   trembl    <- fdata1[REVIEWED==0]
   if (verbose) fdata1 %>% report_n(prefix = 'Per group: drop trembl when swissprot available')

   # trembl groups
   #--------------
   if (nrow(trembl)>0){
      if (verbose) message('')
      if (verbose) trembl  %>% report_n(prefix = 'Trembl groups')

      # Drop fragments when full sequences available
      trembl  %>% magrittr::extract(, IS.FRAGMENT := 0)
      trembl  %>% magrittr::extract(, IS.FRAGMENT:= `PROTEIN-NAMES` %>% stringi::stri_detect_fixed('(Fragment)') %>% as.numeric())
      trembl %<>% magrittr::extract(, .SD[IS.FRAGMENT == min(IS.FRAGMENT)], by = 'feature_id')
      trembl[, IS.FRAGMENT:=NULL]
      if (verbose) trembl  %>% report_n(prefix = 'Per group: drop fragments when full available')

      # Use first sequence per gene
      trembl  %>% magrittr::extract(, N     := .N,                    by = 'feature_id')
      trembl  %>% magrittr::extract(, NGENE := length(unique(GENES)), by = 'feature_id')
      trembl %<>% magrittr::extract(, .SD[1],                         by = c('feature_id', 'GENES'))
      if (verbose) trembl  %>% report_n(prefix = 'Per group/gene: use first accession')
   }

   # swissprot groups
   #-----------------
   if (nrow(swissprot)>0){
      if (verbose) message('')
      if (verbose) swissprot  %>% report_n(prefix = 'Swissprot groups')

      # swissprot groups: collapse similar isoforms (shared accession)
      swissprot  %>% magrittr::extract(, ISOFORM := ISOFORM %>% paste0(collapse = ';'), by = c('feature_id', 'Uniprot accessions'))
      swissprot %<>% unique()
      if (verbose) swissprot  %>% report_n(prefix = 'Per group/accession: collapse spliceforms')

      # swissprot groups: collapse dissimilar isoforms: retain accession with maximum isoforms
      swissprot  %>% magrittr::extract(, NPERACCESSION := ISOFORM %>% stringi::stri_count_fixed(';'), by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, ISOFORM       := ISOFORM %>% paste0(collapse = ';'),         by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = c('feature_id', 'GENES'))
      swissprot %<>% magrittr::extract(, .SD[NPERACCESSION==max(NPERACCESSION)],                      by = c('feature_id', 'GENES'))
      if (verbose) swissprot  %>% report_n(prefix = 'Per group/gene: use spliceform with most accessions')
      swissprot %<>% magrittr::extract(, .SD[1],                                                      by = c('feature_id', 'GENES'))
      if (verbose) swissprot  %>% report_n(prefix = 'Per group/gene: use first spliceform')
   }

   # paralogs
   #---------
   # Collapse paralogs: choose gene with most isoforms
   if (verbose) message('')
   swissprot[, NPERACCESSION:=NULL]
   trembl[, N:=NULL]
   trembl[, NGENE:=NULL]
   fdata1 <- rbind(swissprot, trembl)
   monologs <- fdata1[, .SD[length(unique(GENES))==1], by = c('feature_id')]
   paralogs <- fdata1[, .SD[length(unique(GENES))>1], by = c('feature_id')]
   if (nrow(paralogs)>0){
      if (verbose) paralogs  %>% report_n(prefix = 'Paralog groups')
      paralogs  %>% magrittr::extract(,  NISOFORMS := 1+ISOFORM %>% stringi::stri_count_fixed(';'))
      paralogs  %>% magrittr::extract(,  GENES          := GENES           %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(,  ISOFORM        := ISOFORM         %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(, `PROTEIN-NAMES` := `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = 'feature_id')
      paralogs %<>% magrittr::extract(, .SD[NISOFORMS == max(NISOFORMS)], by = 'feature_id')
      if (verbose) paralogs  %>% report_n(prefix = 'Per group: use paralog with most spliceforms')
      paralogs %<>% magrittr::extract(, .SD[1], by = 'feature_id')
      if (verbose) paralogs  %>% report_n(prefix = 'Per group: use first paralog')
      paralogs  %>% magrittr::extract(, NISOFORMS := NULL)
   }
   fdata1 <- rbind(monologs, paralogs)

   if (verbose) message('')
   if (verbose) fdata1 %>% report_n(prefix = 'All groups (deconvoluted)')

   # Merge back
   nullify_fvars <- function(object, fvars){
      for (curfvar in fvars)   autonomics.import::fdata(object)[[curfvar]] <- NULL
      return(object)
   }
   object %<>% nullify_fvars(fvars = c('Uniprot accessions', 'Protein names', 'Gene names'))
   autonomics.import::fdata(object) %<>% merge(fdata1, by = 'feature_id', sort = FALSE, all.x = TRUE)

   # Rename (MaxQuant style)
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('GENES',         'Gene names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('PROTEIN-NAMES', 'Protein names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('ISOFORM',       'Isoforms')

   # Remove unimportant fvars
   autonomics.import::fdata(object)$REVIEWED  <- NULL
   autonomics.import::fdata(object)$EXISTENCE <- NULL
   if (drop_isoform_info) autonomics.import::fdata(object)$Isoforms <- NULL

   # Return
   object

}

filter_samples_new <- function(object, sample_filter) rlang::eval_tidy(rlang::enexpr(sample_filter), autonomics.import::sdata(object))
# do_filter_samples  <- function(file, sample_filter){
#                          object <- autonomics::read_proteingroups(file)
#                          object %>% filter_samples_new(object, !!rlang::enexpr(sample_filter))
#                       }
# file %>% do_filter_samples((subgroup %>% substr(nchar(.)-2, nchar(.))) == 'STD')

#' Prepare proteingroups
#' @param filter_reverse             logical:   filter out "reverse peptide" groups (used for peptide identification FDR computation)?
#' @param filter_contaminants        logical:   filter out contaminant groups?
#' @param filter_complete_nondetects logical:   filter out proteingroups with no quantification for any sample?
#' @param invert_subgroups           string vector: names of subgroups that require inversion (e.g. WT_KD -> KD_WT)
#' @param log2transform              logical:   log2 transform?
#' @param impute_consistent_nas      logical:   impute consistent NA values (i.e. those replicated in all samples of the same subgroup)?
#' @param deconvolution_fastafile    NULL or string: deconvolute proteingroups using this fastafile (optional)
#' @param verbose                    logical:   message progress?
#' @examples
#' require(magrittr)
#'
#' # Stem cell comparison: triple SILAC design
#'    if (require(autonomics.data)){
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  read_proteingroups()
#'       object %>% prepare_proteingroups(invert_subgroups = c('E_EM', 'E_BM', 'EM_BM'))
#'    }
#'
#' # Stem cell differentiation: internal standard design
#'    if (require(autonomics.data)){
#'       # Read
#'       object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  read_proteingroups()
#'
#'       # Check and remedy subgroup definitions
#'       object$subgroup %>% unique()
#'       object %<>% filter_samples(
#'                      stringi::stri_detect_fixed(subgroup, 'STD') &
#'                      stringi::stri_detect_fixed(subgroup, 'BLANK', negate = TRUE),
#'                      verbose = TRUE)
#'       object$subgroup %<>% stringi::stri_replace_first_fixed('_STD', '')
#'       object$subgroup %<>% factor(c('EM00', 'EM01', 'EM02', 'EM05', 'EM15', 'EM30', 'BM00'))
#'
#'       # Prepare
#'       object %<>% prepare_proteingroups()
#'    }
#'
#' # LFQ design
#'    if (require(graumann.lfq)){
#'       file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'       object <- file %>% autonomics::read_proteingroups()
#'       object %>% autonomics::prepare_proteingroups()
#'    }
#' @importFrom magrittr %>% %<>%
#' @export
prepare_proteingroups <- function(
   object,
   filter_reverse              = TRUE,
   filter_contaminants         = TRUE,
   filter_complete_nondetects  = TRUE,
   invert_subgroups            = character(0),
   log2transform               = TRUE,
   impute_consistent_nas       = stringi::stri_detect_fixed(autonomics.import::guess_maxquant_quantity(object), "Intensity", case_insensitive = TRUE),
   deconvolution_fastafile     = NULL,
   verbose                     = TRUE,
   plot                        = TRUE
){

   # Assert
   assertive.types::assert_is_logical(c(filter_reverse, filter_contaminants, filter_complete_nondetects, log2transform, impute_consistent_nas, verbose, plot))
   assertive.sets::assert_is_subset(invert_subgroups, autonomics.import::subgroup_levels(object))
   if (!is.null(deconvolution_fastafile)) assertive.files::assert_all_are_existing_files(deconvolution_fastafile)

   # Filter samples
   object %<>% autonomics.import::filter_samples_available_for_some_feature(verbose = verbose)

   # Filter features
   if (verbose) autonomics.support::cmessage('\tFilter features')
   if (filter_reverse)              object %<>% autonomics.import::filter_features_("Reverse != '+'", verbose = verbose)
   if (filter_contaminants) {       contaminant_var <- c('Potential contaminant', 'Contaminant') %>% intersect(autonomics.import::fvars(object))
                                    object %<>% autonomics.import::filter_features_(sprintf("`%s` != '+'", contaminant_var), verbose = verbose)}
   if (filter_complete_nondetects)  object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample(verbose = verbose)

   # Transform exprs
   if (verbose) autonomics.support::cmessage('\tTransform exprs')
   object %<>% autonomics.import::invert(subgroups = invert_subgroups)
   object %<>% autonomics.import::zero_to_na(verbose = verbose)
   object %<>% autonomics.import::nan_to_na(verbose = verbose)
   if (log2transform)         object %<>% autonomics.import::log2transform(verbose = verbose)
   if (impute_consistent_nas) object %<>% autonomics.import::impute_consistent_nas(verbose = verbose)
   if (plot)                  object %>%  autonomics.plot::plot_detects_per_subgroup() %>% print()

   # Process fdata
   if (verbose) autonomics.support::cmessage('\tImprove fdata')
   if (verbose) autonomics.support::cmessage("\t\tRename: 'Majority protein IDs' -> 'Uniprot accessions'")
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('Majority protein IDs', 'Uniprot accessions')
   if (!is.null(deconvolution_fastafile)) object %>% deconvolute_proteingroups(deconvolution_fastafile)

   # Return
   object
}

#' Prepare phosphosites
#' @param file                   character(1)
#' @param quantity               NULL or value in names(maxquant_patterns)
#' @param fvars                  character(n)
#' @param min_localization_prob  numeric(1)
#' @param proteingroups_file     character(1)
#' @param verbose                logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/phospho (STY)Sites.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_phosphosites()
#'    file %>% read_phosphosites() %>% prepare_phosphosites()
#' }
prepare_phosphosites <- function(
   object,
   min_localization_prob = 0.75,
   proteingroups_file    = paste0(dirname(file), '/proteinGroups.txt'),
   verbose               = FALSE
){

   # Filter phosphosites
   object %<>% autonomics.import::filter_features(!stringi::stri_detect_fixed(`Protein group IDs`, ';'), verbose = TRUE)
   assertive.numbers::assert_all_are_in_range(min_localization_prob, 0, 1)
   object %<>% autonomics.import::filter_features_(sprintf("`Localization prob` >= %s", as.character(min_localization_prob)),verbose = TRUE)

   # Return
   object
}

