#---------------------------------------------------------------------------
# 
#                       guess_compounddiscoverer_quantity
#
#---------------------------------------------------------------------------


#' compound discoverer quantity patterns
#' @examples
#' COMPOUNDDISCOVERER_PATTERNS
#' @export
COMPOUNDDISCOVERER_PATTERNS <- c(
    `area`           = '^Area: (.+)$',
    `normalizedarea` = '^Norm\\. Area: (.+)$')

#' Guess compound discoverer quantity from snames
#'
#' @param x character vector
#' @return  string: value from names(COMPOUNDDISCOVERER_PATTERNS)
#' @examples
#' \dontrun{
#' # file
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     guess_compounddiscoverer_quantity(file)
#' }
#'
#' # character vector
#'     x <- "Area: 20230908_F143_HILICNEG.raw (F11)"
#'     guess_compounddiscoverer_quantity(x)
#'
#'     x <- "Norm. Area: 20230908_F143_HILICNEG.raw (F11)"
#'     guess_compounddiscoverer_quantity(x)
#' @export
guess_compounddiscoverer_quantity <- function(x){
# Assert
    assert_is_character(x)
# read if x is filename
    if (is_scalar(x)){
        if (is_existing_file(x)){  
            x <- names(fread(x, header = TRUE, nrows = 1))
        }
    }
# guess from character vector
    for (quantity in names(COMPOUNDDISCOVERER_PATTERNS)){
        pattern <- COMPOUNDDISCOVERER_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#---------------------------------------------------------------------------
# 
#                       .read_compounddiscoverer
#
#---------------------------------------------------------------------------

#' Read compound discoverer files as-is
#' @param file               compoumd discoverer file
#' @param quantity           string
#' @param colname_formatter  function for column reformatting
#' @param verbose            TRUE / FALSE
#' @return data.table
#' @export
.read_compounddiscoverer <- function(
    file, quantity = guess_compounddiscoverer_quantity(file), 
    colname_formatter = NULL,
    modus_extractor = NULL,
    verbose = TRUE){
    . <- Name <- Value <- rowid <- Acquisition <- NULL
# Assert
    assert_compounddiscoverer_output(file)
    assert_is_subset(quantity, names(COMPOUNDDISCOVERER_PATTERNS))
    if (!is.null(colname_formatter)) assert_is_function(colname_formatter)
    assert_is_function(modus_extractor)
# Read
    if (verbose)  cmessage('%sRead%scompound discoverer output   %s', spaces(4), spaces(6), file)
    cddt <- fread(file, integer64 = 'numeric') %>% 
      un_int64()
    n0 <- nrow(cddt)
# Determine MS modus
    cddt %<>%
      dplyr::mutate(
        `MS Modus` = modus_extractor(names(.)) %>%
          unique() %>%
          assert_are_same_length(1))
# Determine ID qualifier
    cddt %<>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        `ID Quality` =  if (
          all(is.na(dplyr::c_across(dplyr::ends_with("Best Match")))))
          { NA } else {
            max(
              dplyr::c_across(dplyr::ends_with("Best Match")), na.rm = TRUE) },
        `ID Quality` =  dplyr::case_when(
          `ID Quality` > 80    ~ 2,
          `ID Quality` > 60    ~ 3,
          TRUE                 ~ 4)) %>%
      dplyr::ungroup()
# Extract require columns
    pattern <- COMPOUNDDISCOVERER_PATTERNS[[quantity]]
    anncols <- c(
      'Name', 'Compounds ID', 'MS Modus', 'Formula', 'Annot. DeltaMass [ppm]',
      'Calc. MW', 'm/z', 'RT [min]', 'ID Quality', 'mzCloud Best Match',
      'mzVault Best Match',
      Name)
    anncols %<>% intersect(names(cddt))
    valcols <- grep(pattern, names(cddt), value = TRUE)
    cddt %<>% extract(, c(anncols, valcols), with = FALSE)
# Format column names
    if (!is.null(colname_formatter))
      cddt %<>% magrittr::set_colnames(colname_formatter(names(.)))
# Uniquify Name; choose detection with largest area
    cddt %<>%
      dplyr::mutate(rowid = dplyr::row_number()) %>%
      dplyr::relocate(rowid) %>%
      tidyr::pivot_longer(
        cols      = dplyr::matches(COMPOUNDDISCOVERER_PATTERNS[quantity]),
        names_to  = "Acquisition",
        values_to = "Value") %>%
      dplyr::group_by(Name) %>%
      dplyr::arrange(dplyr::desc(Value), by_group = TRUE)
    if (verbose) {
       n <- cddt %>%
         dplyr::summarise(n = length(unique(rowid))) %>% 
         dplyr::filter(n > 1) %>%
         nrow()
       if (n > 0) cmessage('%sMerge %s multiple detects by max %s', spaces(18), n, quantity)
    } 
    cddt %<>%
      dplyr::filter(rowid == rowid[1]) %>%
        dplyr::select(-rowid) %>%
      tidyr::pivot_wider(
        names_from = c(Acquisition),
        values_from = Value) %>%
      dplyr::ungroup() %>%
      as.data.table(integer64 = "numeric")
# Return
    cddt[]
}

#---------------------------------------------------------------------------
#
#                 dequantify_compounddiscoverer
#
#---------------------------------------------------------------------------
#' dequantify_compounddiscoverer compound discoverer snames
#' 
#' Drop quantity.
#' 
#' `Norm. Area: 20230908_F143_HILICNEG.raw (F11)  ->  20230908_F143_HILICNEG.raw (F11)`
#'       `Area: 20230908_F143_HILICNEG.raw (F11)  ->  20230908_F143_HILICNEG.raw (F11)`
#' @param x        `character`
#' @param quantity `'area',              'normalizedarea'`
#' @param verbose  `TRUE` or `FALSE`
#' @return `character`
#' @examples
#' dequantify_compounddiscoverer("Norm. Area: 20230908_F143_HILICNEG.raw (F11)") # Norm. Area
#' dequantify_compounddiscoverer("Area: 20230908_F143_HILICNEG.raw (F11)")       # Area
#' @md
#' @export
dequantify_compounddiscoverer <- function(
    x, quantity = guess_compounddiscoverer_quantity(x), verbose  = FALSE
){
  pattern <- COMPOUNDDISCOVERER_PATTERNS[[quantity]]
  cleanx <- stri_replace_first_regex(x, pattern, '$1')
  if (verbose)  cmessage('%sStandardize snames: %s  ->  %s', spaces(14), x[1], cleanx[1])
  return(cleanx)
}

#' merge compound discoverer files
#' @param x        `list`
#' @param quantity `'area', 'normalizedarea'`
#' @param verbose  `TRUE` or `FALSE`
#' @return `data.table`
#' @export
merge_compounddiscoverer <- function(x, quantity = NULL, verbose = TRUE) {
    Modus <- Name <- Value <- Acquisition <- NULL
    assert_is_list(x)
    assert_has_names(x)
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(COMPOUNDDISCOVERER_PATTERNS))
    assert_all_have_setidentical_colnames(x)
# Introduce origin identifiers
    x <- lapply(
      names(x),
      function(y){
        x[[y]] %>%
          dplyr::mutate(
            Modus = y,
            ID    = paste(Modus, seq(dplyr::n()), sep = "_"))
      }) %>%
      magrittr::set_names(names(x))
# Ensure proper treatment of data from single source file
    if (length(x) == 1) return(utils::head(x, n = 1))
# Pivot longer and combine
    if (verbose)  cmessage('%sMerge%scompount discoverer output%s(choose origin with max %s)', spaces(4), spaces(5), spaces(3), quantity)
    x %<>%
      lapply(
        tidyr::pivot_longer,
        cols      = dplyr::matches(COMPOUNDDISCOVERER_PATTERNS[quantity]),
        names_to  = "Acquisition",
        values_to = "Value") %>%
      dplyr::bind_rows()
# Processing of merged data
    ## Group by `Name`
    x %<>% dplyr::group_by(Name)
    ## Add column indicating all sources
    x %<>% dplyr::mutate(Modi   = paste(unique(Modus), collapse = ";"))
    ## Identify source providing max. *Area and filter by it
    x %<>%
      dplyr::arrange(dplyr::desc(Value), .by_group = TRUE) %>%
      dplyr::filter(Modus == Modus[1])
# Pivot data wide again
    x %>% 
      tidyr::pivot_wider(
        names_from = c(Acquisition),
        values_from = Value) %>%
      dplyr::ungroup() %>%
      as.data.table(integer64 = "numeric")
}

#---------------------------------------------------------------------------
#
#                   read_compounddiscoverer
#
#---------------------------------------------------------------------------

#' Read maxquant proteingroups
#' @param dir                compound discoverer output directory
#' @param files              compound discoverer output files
#' @param colname_formatter  function for column reformatting on file read
#' @param quantity           'area', 'normalizedarea' or NULL
#' @param nonames            TRUE or FALSE: retain compunds without Names?
#' @param subgroups          NULL or string vector : subgroups to retain
#' @param impute             TRUE or FALSE: impute group-specific NA values?
#' @param plot               TRUE or FALSE: plot ?
#' @param label              fvar
#' @param pca                TRUE or FALSE: run pca ?
#' @param pls                TRUE or FALSE: run pls ?
#' @param fit                model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula            model formula
#' @param block              model blockvar: string or NULL
#' @param coefs              model coefficients    of interest: character vector or NULL
#' @param contrasts          coefficient contrasts of interest: character vector or NULL
#' @param palette            color palette : named character vector
#' @param verbose            TRUE or FALSE : message ?
#' @return SummarizedExperiment
#' @export
read_compounddiscoverer <- function(
    dir = getwd(),
    files                 = list.files(
                              path = dir, pattern = "(RP|HILIC).*\\.csv$",
                              full.names          = TRUE),
    colname_regex         = "^(.*)\\d{8,8}_+(.*)_+((HILIC|RP)(NEG|POS))\\.raw.*$",
    colname_formatter     = function(x)
                              stringi::stri_replace_first_regex(
                                x, colname_regex, "$1$2"),
    modus_extractor       = function(x)
                              stringi::stri_subset_regex(
                                x, colname_regex) %>%
                                stringi::stri_replace_first_regex(
                                  colname_regex, "$3"),
    quantity              = NULL,
    nonames               = FALSE,
    sample_filter_strings = c("blank", "QC", "RS"),
    subgroups             = NULL,
    impute                = FALSE,
    plot                  = FALSE,
    label                 = 'feature_id',
    pca                   = plot,
    pls                   = plot,
    fit                   = if (plot) 'limma' else NULL,
    formula               = ~ subgroup,
    block                 = NULL,
    coefs                 = NULL,
    contrasts             = NULL,
    palette               = NULL,
    verbose               = TRUE)
{
    Name <- NULL
# Select files
    assert_is_a_string(dir)
    assert_all_are_dirs(dir)
    assert_is_not_null(files)
    assert_all_are_greater_than(length(files), 0)
    files <- if (!has_names(files)) {
      magrittr::set_names(files, tools::file_path_sans_ext(basename(files)))
    } else { files }
# Assert
    for (fl in files) assert_compounddiscoverer_output(fl)
    if (!is.null(colname_regex)) assert_is_a_string(colname_regex)
    if (!is.null(colname_formatter)) assert_is_function(colname_formatter)
    assert_is_function(modus_extractor)
    quantity <- if (is.null(quantity)) {
      sapply(files, guess_compounddiscoverer_quantity) %>% unique()
    } else { quantity }
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(COMPOUNDDISCOVERER_PATTERNS))
    assert_is_a_bool(nonames)
    if (!is.null(sample_filter_strings))
      assert_is_character(sample_filter_strings)
    assert_is_a_bool(verbose)
# Read/Merge
    cddt <- lapply(
      files,
      .read_compounddiscoverer, colname_formatter = colname_formatter,
      modus_extractor = modus_extractor, quantity = quantity,
      verbose = verbose) %>%
      merge_compounddiscoverer(quantity = quantity, verbose = verbose)
# Currently intro of `feature_id` only, no `annotate_maxquant` equivalent (yet)
    cddt[, feature_id := `Compounds ID`]
# SumExp
    if (verbose)  cmessage('%sSumExp', spaces(4))
    pattern <- COMPOUNDDISCOVERER_PATTERNS[[quantity]]
    cdmat <- mqdt_to_mat(cddt, pattern, verbose = verbose)
    cddt %<>% extract(, names(cddt) %>% setdiff(colnames(cdmat)), with = FALSE)
    object <- list(cdmat) %>%
      magrittr::set_names(paste0('log2', quantity)) %>%
      SummarizedExperiment(rowData = cddt)
# dequantify_compounddiscoverer.
    object$cdcol <- colnames(object)
    colnames(object) %<>%
      dequantify_compounddiscoverer(quantity = quantity, verbose = verbose)
# process_compounddiscoverer
    object %<>% process_compounddiscoverer(
        nonames = nonames, sample_filter_strings = sample_filter_strings, 
        subgroups = subgroups, impute = impute, verbose = verbose)
    assays <- c(assayNames(object)[1])
    for (assay in assays)  object %<>% add_assay_means(assay)
# Process statistically - if requested
    object %<>% analyze(
        pca       = pca,          pls     = pls,
        fit       = fit,          formula = formula,
        block     = block,        coefs   = coefs,
        contrasts = contrasts,    plot    = plot,
        label     = label,        palette = palette,
        verbose   = verbose )
    object
}

#----------------------------------------------------------------------------
#
#                 process_compounddiscoverer
#
#----------------------------------------------------------------------------

process_compounddiscoverer <- function(
    object, nonames, sample_filter_strings, subgroups, impute, verbose
){
# Deal with nonames
    if (!nonames) object %<>% filter_features(
      !stringi::stri_isempty(Name), verbose = verbose)
# Deal with sample filter strings
    if (!is.null(sample_filter_strings)) {
      if (verbose)  cmessage(
        '%sFiltering against %s in sample names',
        spaces(14), paste0("'", sample_filter_strings, "'", collapse = ", "))
      for (str in sample_filter_strings)
        object %<>% filter_samples(
          stringi::stri_detect_fixed(cdcol, str, negate = TRUE),
          verbose = verbose)
    }
# Infer Subgroup
    object$sample_id <- colnames(object)
    object %<>% add_subgroup(verbose = verbose)
    cols <- c('sample_id', 'subgroup', 'replicate', 'cdcol') %>%
      intersect(svars(object))
    sdt(object) %<>% pull_columns(cols)
# Samples
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    if (!is.null(subgroups)) {
        assert_is_subset(subgroups, as.character(object$subgroup))
        object %<>% filter_samples(subgroup %in% subgroups, verbose = verbose)
    }
# Features
    analysis(object)$nfeatures <- c(nrow(object))
    object %<>% rm_missing_in_all_samples(verbose = verbose)
# Impute
    if (impute) object %<>% impute(plot = FALSE)
    object
}
