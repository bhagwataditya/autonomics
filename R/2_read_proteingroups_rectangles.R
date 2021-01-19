add_maxquant_sdata <- function(
    object, samplefile = default_samplefile(object), verbose
){
    snames(object) %<>% stri_replace_last_fixed('___1', '') # PHOSPHOSITES
    object %<>% standardize_maxquant_snames(verbose = verbose)
    object %<>% demultiplex(verbose = verbose)
    object %<>% merge_samplefile(samplefile = samplefile, verbose = verbose)
    object
}


.read_proteingroup_rectangles <- function(
    file, quantity = guess_maxquant_quantity(file),
    samplefile = default_samplefile(file),
    select_subgroups = NULL, invert_subgroups = character(0), verbose = TRUE
){
# Scan
    dt <- fread(file, integer64 = 'numeric', header = FALSE)
    columns <- as.character(dt[1])
    fvars <- intersect(columns, PROTEINGROUP_FVARS)
    fid_rows <- 2:nrow(dt);    fid_cols <- which(columns == 'id')
    sid_rows <- 1;             pattern  <- MAXQUANT_PATTERNS[[quantity]]
    sid_cols <- which(stri_detect_regex(columns, pattern))
    expr_rows  <- 2:nrow(dt);  expr_cols  <- sid_cols
    fvar_rows  <- 1;           fvar_cols  <- match(fvars, columns)
    fdata_rows <- 2:nrow(dt);  fdata_cols <- fvar_cols
# Read
    object <- read_omics(file,
                        fid_rows   = fid_rows,    fid_cols   = fid_cols,
                        sid_rows   = sid_rows,    sid_cols   = sid_cols,
                        expr_rows  = expr_rows,   expr_cols  = expr_cols,
                        fvar_rows  = fvar_rows,   fvar_cols  = fvar_cols,
                        fdata_rows = fdata_rows,  fdata_cols = fdata_cols,
                        transpose  = FALSE,       verbose    = verbose)
# Features: parse/filter
    if (verbose)  message('\tPrepare samples')
    if ('Gene names' %in% fvars(object))    fvars(object) %<>%
        stri_replace_first_fixed('Gene names', 'feature_name')

# Samples: parse/filter
    if (verbose)  message('\tPrepare samples')
    object %<>% add_maxquant_sdata(verbose=verbose, samplefile=samplefile)
    object %<>% add_subgroup()
    object %<>% filter_maxquant_samples(
                    select_subgroups = select_subgroups, verbose)
# Transform
    if (verbose)  message('\tLog2 transform')
    exprs(object) %<>% zero_to_na(verbose = verbose)
    exprs(object) %<>% nan_to_na(verbose = verbose)
    object %<>% log2transform(verbose=verbose)
    object %<>% invert(invert_subgroups)
# Return
    metadata(object)$quantity <- quantity
    metadata(object)$platform <- 'maxquant'
    object
}

.read_phosphosite_rectangles <- function(
    file, quantity = guess_maxquant_quantity(file),
    samplefile = default_samplefile(file),
    select_subgroups = subgroup_levels(object), invert_subgroups = character(0),
    verbose=TRUE
){
# Scan
    dt <- fread(file, integer64 = 'numeric', header = FALSE)
    columns <- as.character(dt[1])
    value_cols <- phospho_expr_columns(columns, quantity)
    fvars <- intersect(columns, PHOSPHOSITE_FVARS)
    fvar_cols  <- which(columns %in% fvars)
    fid_rows   <- 2:nrow(dt);     fid_cols   <- which(columns == 'id')
    sid_rows   <- 1;              sid_cols   <- value_cols
    expr_rows  <- 2:nrow(dt);     expr_cols  <- value_cols
    fvar_rows  <- 1;              fvar_cols  <- fvar_cols
    fdata_rows <- 2:nrow(dt);     fdata_cols <- fvar_cols
# Read
    object  <- read_omics(file,
                        fid_rows   = fid_rows,    fid_cols   = fid_cols,
                        sid_rows   = sid_rows,    sid_cols   = sid_cols,
                        expr_rows  = expr_rows,   expr_cols  = expr_cols,
                        fvar_rows  = fvar_rows,   fvar_cols  = fvar_cols,
                        fdata_rows = fdata_rows,  fdata_cols = fdata_cols,
                        transpose  = FALSE,       verbose    = verbose)
    snames(object)   %<>% stri_replace_all_fixed('___1', '')
    object$sample_id %<>% stri_replace_all_fixed('___1', '')
    metadata(object)$quantity <- quantity
    metadata(object)$platform <- 'maxquant'
# Samples: Parse/Filter
    if (verbose)  message('\tPrepare samples')
    object %<>% add_maxquant_sdata(verbose=verbose, samplefile=samplefile)
    object %<>% add_subgroup()
    object %<>% filter_maxquant_samples(
                    select_subgroups = select_subgroups, verbose)
# Log2 transform
    if (verbose)  message('\tLog2 transform')
    exprs(object) %<>% zero_to_na(verbose = verbose)
    exprs(object) %<>% nan_to_na(verbose = verbose)
    object %<>% log2transform(verbose=verbose)
    object %<>% invert(invert_subgroups)
# Return
    object
}


# proteinfile <- download_data('billing19.proteingroups.txt')
# select <- c('EM00', 'EM01', 'EM02', 'EM05', 'EM15', 'EM30', 'BM00')
# select %<>% paste0('_STD')
# object <- read_proteingroup_rectangles(proteinfile, select=select)
# object <- read_proteingroups(proteinfile, select=select)
read_proteingroup_rectangles <- function(
    proteinfile, quantity = guess_maxquant_quantity(proteinfile),
    samplefile = default_samplefile(proteinfile),
    select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, fastafile = NULL, invert_subgroups = character(0),
    impute = stri_detect_regex(quantity, "[Ii]ntensity"),
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object),
    verbose      = TRUE, plot = TRUE
){
# Assert
    assert_all_are_existing_files(proteinfile)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    if (!is.null(fastafile)) assert_all_are_existing_files(fastafile)
    assert_is_a_bool(verbose)
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
# Read
    object <- .read_proteingroup_rectangles(proteinfile, quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose = verbose)
# Prepare
    object %<>% filter_maxquant_features(reverse = reverse,
                    contaminants = contaminants, verbose = verbose)
    object %<>% rename_proteingroup_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=impute, verbose=verbose, plot=plot)
# Analyze
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Return
    if (plot)  plot_samples(object)
    object
}

# require(magrittr)
# phosphofile <- download_data("billing19.phosphosites.txt")
# proteinfile <- download_data("billing19.proteingroups.txt")
# select <- c('EM00', 'EM01', 'EM02', 'EM05', 'EM15', 'EM30', 'BM00')
# select %<>% paste0('_STD')
# object <- read_phosphosite_rectangles(phosphofile, proteinfile, select=select)
# object <- read_phosphosites(phosphofile, proteinfile, select=select)
read_phosphosite_rectangles <- function(
    phosphofile,
    proteinfile = paste0(dirname(phosphofile), '/proteinGroups.txt'),
    quantity = guess_maxquant_quantity(phosphofile),
    samplefile = default_samplefile(proteinfile),
    select_subgroups = NULL,
    contaminants = FALSE, reverse = FALSE, min_localization_prob = 0.75,
    fastafile = NULL, invert_subgroups = character(0),
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object),
    verbose      = TRUE, plot = TRUE
){
# Assert
    `Protein group IDs` <- `Localization prob` <- NULL
    assert_all_are_existing_files(c(phosphofile, proteinfile))
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
# Read
    proteingroups <- .read_proteingroup_rectangles(file=proteinfile, quantity = quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose=verbose)
    object  <- .read_phosphosite_rectangles(file = phosphofile, quantity = quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose = verbose)
    object %<>% filter_maxquant_features(
                    reverse = reverse, contaminants = contaminants,
                    min_localization_prob = min_localization_prob,
                    verbose = verbose)
    object %<>% add_occupancies(proteingroups, verbose)
# Prepare
    object %<>% rename_phospho_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=FALSE,verbose=verbose,plot=plot)
# Contrast
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Return
    if (plot)  plot_samples(object)
    object
}
