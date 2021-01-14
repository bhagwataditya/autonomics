#=========================================================
#
#         download_gtf (alternative: biomartr::getGTF)
#             make_gtf_url
#                 release_to_build
#
#=========================================================

release_to_build <- function(release, organism){
    if        (organism == 'Homo sapiens'){
        if (release >= 76)  'GRCh38'   else 'GRCh37'  }

    else if   (organism == 'Mus musculus'){
        if (release >= 68)  'GRCm38'   else 'NCBIM37'  }

    else if   (organism == 'Rattus norvegicus'){
        if (release >= 80)  'Rnor_6.0' else 'Rnor_5.0'  }
}


#' Make link to GTF file
#' @param organism 'Homo sapiens', 'Mus musculus', or 'Rattus norvegicus'
#' @param release   number
#' @examples
#' make_gtf_url(organism = 'Homo sapiens',            release = 95)
#' make_gtf_url(organism = 'Mus musculus',            release = 95)
#' make_gtf_url(organism = 'Sacharomyces cerevisiae', release = 100)
#' @noRd
make_gtf_url <- function(organism, release){
    sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz',
            release,
            stri_replace_first_fixed(tolower(organism), ' ', '_'),
            stri_replace_first_fixed(        organism,  ' ', '_'),
            release_to_build(release, organism),
            release)
}


#' Download GTF file
#'
#' Download GTF file with feature annotations
#' @param organism  'Homo sapiens', 'Mus musculus' or 'Rattus norvegicus'
#' @param release    GTF release (number)
#' @param gtffile    string: path to local GTF file
#' @return gtffile path
#' @examples
#' organism <- 'Homo sapiens'
#' # download_gtf(organism)
#' @export
download_gtf <- function(
    organism,
    release = 100,
    gtffile = sprintf("~/autonomicscache/gtf/%s",
        basename(make_gtf_url(organism, release) %>% substr(1, nchar(.)-3)))
){
    assert_is_subset(organism,
                    c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'))
    . <- NULL
    remote <- make_gtf_url(organism, release)

    if(file.exists(gtffile)){
        cmessage("GTF file already available at %s", gtffile)
    } else {
        cmessage("download   %s'" , remote)
        cmessage("to         %s", gtffile )
        dir.create(dirname(gtffile), showWarnings = FALSE, recursive = TRUE)
        gtffile %<>% paste0('.gz')
        tryCatch(expr = { download.file(
                                url = remote, destfile = gtffile, quiet = TRUE)
                            gunzip(gtffile,  remove = TRUE, overwrite = TRUE)},
                error = function(cond){
                    message('failed     repeat manually using browser\n');
                    message(cond)})
    }
    gtffile %<>% stri_replace_last_fixed('.gz', '')
    invisible(gtffile)
}


#==============================================================================
#
#                                      gtf
#               add_genenames: geneid -----> genename     (optional)
#
#==============================================================================


#' Add gene names
#'
#' Add gene names to SummarizedExperiment
#'
#' @param object     SummarizedExperiment
#' @param gtffile    string: path to gtffile
#' @param verbose    TRUE or FALSE
#' @noRd
add_genenames <- function(object, gtffile, verbose = TRUE){

    gene_name <- NULL

    if (is.null(gtffile)) return(object)

    if (verbose) message('\t\t\tRead ', gtffile)
    gtfdt <- rtracklayer::import(gtffile) %>%
            GenomicRanges::as.data.frame() %>%
            data.table()

    if (verbose) message('\t\t\tMap gene_id to gene_name')
    x <- fdata(object)$gene_name
    gtfdt %<>% extract(, .SD[1], by = 'gene_id')
    fdata(object)$feature_name <-  gtfdt[x, gene_name, on = "gene_id"]

    object
}


#==============================================================================
#
#                       count_reads
#                           entrezg_to_symbol
#
#==============================================================================


entrezg_to_symbol <- function(x, genome){

    assert_is_subset(genome, c('mm9', 'mm10', 'hg19', 'hg38'))

    orgdb <- if (genome %in% c('mm10', 'mm9')){
                if (!requireNamespace('org.Mm.eg.db', quietly = FALSE)){
                    stop("First: BiocManager::install('org.Mm.eg.db')")}
                orgdb <- org.Mm.eg.db::org.Mm.eg.db

            } else if (genome %in% c('hg19', 'hg38')){
                if (!requireNamespace('org.Hs.eg.db', quietly = FALSE)){
                    stop("First: BiocManager::install('org.Hs.eg.db')")}
                orgdb <- org.Hs.eg.db::org.Hs.eg.db
            }
    x %<>% as.character()
    suppressMessages(y <- AnnotationDbi::mapIds(orgdb, x, 'SYMBOL', 'ENTREZID'))
    y %<>% unname()
    y
}



count_reads <- function(files, paired, nthreads, genome){
# Common args
    . <- NULL
    args <- list(files = files, isPaired = paired, nthreads = nthreads)
# Inbuilt genome
    if (genome %in% c('mm10', 'mm9', 'hg38', 'hg19')){
        args %<>% c(list(annot.inbuilt = genome))
        fcounts <- do.call(featureCounts, args)
        fcounts$annotation$gene_name <-
            entrezg_to_symbol(fcounts$annotation$GeneID, genome)
# User GTF
    } else {
        assert_all_are_existing_files(genome)
        args %<>% c(list(annot.ext = genome, isGTFAnnotationFile = TRUE,
                        GTF.attrType.extra  = 'gene_name'))
        fcounts <- do.call(featureCounts, args)
    }
# Rename, Select, Return
    names(fcounts$annotation) %<>% gsub('GeneID',   'feature_id',   .)
    names(fcounts$annotation) %<>% gsub('gene_name', 'feature_name', .)
    fcounts$annotation %<>% extract(,
                intersect(names(.), c('feature_id','feature_name')),
                drop = FALSE)
    fcounts$annotation$feature_id %<>% as.character()
    rownames(fcounts$annotation) <- fcounts$annotation$feature_id
    fcounts
}


#=============================================================================
#
#               counts_to_cpm/cpm_to_counts
#                   scaledlibsizes
#
#=============================================================================


#' @rdname counts_to_cpm
#' @noRd
scaledlibsizes <- function(counts){
    colSums(counts) * edgeR::calcNormFactors(counts)
}


#' Convert between counts and cpm (counts per million scaled reads)
#' @param counts           count matrix
#' @param cpm              cpm matrix (counts per million scaled reads)
#' @param scaled_libsizes  numeric vector: scaled library sizes
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_counts(file)
#'
#' scaled_libsizes <- scaledlibsizes(exprs(object))
#' cpm <- counts_to_cpm(counts(object), scaled_libsizes)
#'
#' counts  <- cpm_to_counts(cpm, scaled_libsizes)
#' sum(counts(object) - counts)
#' @return cpm matrix
#' @noRd
counts_to_cpm <- function(counts, scaled_libsizes = scaledlibsizes(counts)){
    t(t(counts + 0.5)/(scaled_libsizes + 1) * 1e+06)
}


#' @noRd
cpm_to_counts <- function(cpm, scaled_libsizes){
    1e-06 * t(t(cpm) * (scaled_libsizes + 1)) - 0.5
}

#=============================================================================
#
#                  compute_precision_weights
#                      compute_precision_weights_once
#
#=============================================================================

compute_precision_weights_once <- function(
    object,
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    plot    = TRUE, ...
){
    counts   <- counts(object)
    scaled_libsizes <- scaledlibsizes(counts)

    design <- create_design(object, formula=!!enquo(formula))
    counts %>%
    voom(design = design, lib.size = scaled_libsizes, plot = plot, ...) %>%
    extract2('weights') %>%
    set_rownames(rownames(object)) %>%
    set_colnames(colnames(object))
}


#' Compute voom precision weights
#' @param object  SummarizedExperiment: exprs(.) returns log2cpm, counts(.)
#'                returns raw counts.
#' @param formula                formula to create design matrix (using svars)
#' @param plot    TRUE (default) or FALSE
#' @param ...     passed to limma::voom() -> limma::lmFit()
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_counts(file, plot = FALSE)
#' compute_precision_weights(object)[1:3, 1:3]
#'
#' object$subgroup <- 'billing16'
#' compute_precision_weights(object)[1:3, 1:3]
#'
#' object$subgroup <- object$sample_id
#' compute_precision_weights(object)[1:3, 1:3]
#'
#' file <- download_file('billing19.rnacounts.txt')
#' object <- read_counts(file)
#' counts(object) <- exprs(object)
#' compute_precision_weights(object)[1:3, 1:3]
#' @noRd
compute_precision_weights <- function(
    object,
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    plot    = TRUE
){
# Assert & message
    assert_is_not_null(counts(object))
    formula <- enquo(formula)
# Estimate precision weights
    has_block <- has_complete_block_values(object)
    weights   <- compute_precision_weights_once(
                    object, formula = !!formula, plot = !has_block & plot)
# Update precision weights using block correlation
    if (has_block){
        log2cpm <- log2(counts_to_cpm(counts(object)))
        correlation <-  duplicateCorrelation(
                            log2cpm, design = create_design(object, !!formula),
                            block = object$block, weights = weights)
        correlation %<>% extract2('consensus')
        weights <- compute_precision_weights_once(
                    object, formula = !!enquo(formula), block = object$block,
                    correlation = correlation, plot = plot)
    }
# Return
    dimnames(weights) <- dimnames(object)
    weights
}


explicitly_compute_precision_weights_once <- function(
    object, formula, plot = TRUE, ...
){
# Extract
    log2cpm  <- exprs(object)
    lib.size <- scaledlibsizes(counts(object))
    design   <- create_design(object, formula=!!enquo(formula))
# Assert
    n <- nrow(log2cpm)
    if (n < 2L) stop("Need at least two genes to fit a mean-variance trend")
# Fit linear model
    fit <- lmFit(log2cpm, design=design, ...)
# Predict
    if (is.null(fit$Amean)) fit$Amean <- rowMeans(log2cpm, na.rm = TRUE)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[seq_len(fit$rank)]
        fitted.log2cpm  <-  fit$coef[, j, drop = FALSE] %*%
                            t(fit$design[,j, drop = FALSE])
    } else {
        fitted.log2cpm <- fit$coef %*% t(fit$design)
    }
    fitted.log2count <-(2^fitted.log2cpm) %>% cpm_to_counts(lib.size) %>% log2()
# Fit mean-variance trend
    mean.log2count <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
    sdrt.log2count <- sqrt(fit$sigma)        # mean log2 count  &  sqrtsd(resid)
    all.identical <- matrixStats::rowVars(log2cpm)==0
    if (any(all.identical)) {
        mean.log2count <- mean.log2count[!all.identical]
        sdrt.log2count <- sdrt.log2count[!all.identical]
    }
    l <- lowess(mean.log2count, sdrt.log2count, f = 0.5)
    f <- approxfun(l, rule = 2)
# Compute precision weights
    w <- 1/f(fitted.log2count)^4     # f(.) = sqrt(sd(.)) --> f(.)^4 = var(.)
    dim(w) <- dim(fitted.log2count)
# Plot
    if (plot) {
        plot(mean.log2count, sdrt.log2count, xlab = "mean log2count",
            ylab = "sdrt log2count", pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }
# Return
    return(w)
}


#==============================================================================
#
#                   preprocess_counts
#                       filter_low_count_features
#
#==============================================================================

filter_low_count_features <- function(object,filter_count,verbose){
    idx <- filterByExpr(counts(object),
                        group     = object$subgroup,
                        lib.size  = scaledlibsizes(counts(object)),
                        min.count = filter_count)

    if (verbose) message('\t\tKeep ', sum(idx), '/', length(idx),
                ' features: count >= ',
                filter_count, ' in at least some samples')

    object[idx, ]
}


preprocess_counts <- function(
    object,
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    filter_count = 10,
    verbose      = TRUE,
    plot         = TRUE
){
# Filter lowly expressed features
    sdata(object)$libsize <- colSums(counts(object))
    object %<>% filter_low_count_features(
                    filter_count, verbose = verbose)
    sdata(object)$libsize.filtered <- colSums(counts(object))
    sdata(object)$libsize.scaled   <- scaledlibsizes(counts(object))
# Normalize: counts -> cpm
    if (verbose) message('\t\tTMM normalize and log2transform')
    exprs(object) <- counts(object) %>% counts_to_cpm() %>% log2()
# Compute and add precision weights
    assays(object)$weights <- compute_precision_weights(
        object, formula = !!enquo(formula), plot = plot)
# Return
    if (verbose){
        message('\t\tReturn object: counts(object)  = counts')
        message('\t\t               exprs(object)   = log2cpm')
        message('\t\t               weights(object) = voom weights (see plot)')
        message('\t\t               object$subgroup = subgroup values')
    }
    object
}


#=============================================================================
#
#                           read_bam
#                           read_counts
#
#=============================================================================


#' Read RNAseq SAM/BAM files into SummarizedExperiment
#'
#' @param bamdir        SAM/BAM dir path (string)
#' @param paired        whether reads are paired end (TRUE/FALSE)
#' @param genome        string: either "mm10", "hg38" etc. or a GTF file
#' @param nthreads      no of cores to be used by Rsubread::featureCounts()
#' @param filter_count  min feature count required by at least one sample
#' @param formula       formula to create design matrix (using svars)
#' @param contrastdefs  contrast definition vector/matrix/list
#' @param verbose       TRUE/FALSE
#' @param plot          TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # in-built genome
#' bamdir <- download_data("billing16.bam.zip")
#' # object <- read_bam(bamdir, paired=TRUE, genome="hg38")
#' # gtffile <- download_gtf("Homo sapiens")
#' # object <- read_bam(bamdir, paired=TRUE, genome=gtffile)
#' @author Aditya Bhagwat, Shahina Hayat
#' @export
read_bam <- function(bamdir, paired, genome, nthreads = detectCores(),
    filter_count = 10,
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object),
    verbose      = TRUE, plot = TRUE
){
# Assert
    assert_all_are_existing_files(bamdir)
    assert_is_a_bool(paired)
    assert_is_a_string(genome)
    assert_is_a_number(nthreads)
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
# Count reads
    files <- list.files(bamdir, pattern = ".sam$|.bam$", full.names = TRUE,
                        recursive = TRUE)
    fcounts <- count_reads(files, paired, nthreads=nthreads, genome=genome)
# Forge SummarizedExperiment
    filenames   <- basename(file_path_sans_ext(files))
    subdirnames <- basename(dirname(files))
    sample_names <- if (has_no_duplicates(filenames)){           filenames
                    } else if (has_no_duplicates(subdirnames)){  subdirnames
                    } else {           paste0(subdirnames, '_', filenames)}
    object <- SummarizedExperiment(assays = list(
                counts = fcounts$counts %>% set_colnames(sample_names)))
    metadata(object)$platform <- 'rnaseq'
    metadata(object)$file <- bamdir
# Add fdata
    message("\t\tAdd fdata")
    fcounts$annotation$feature_id %<>% as.character()
    rowData(object)  <- fcounts$annotation
# Add design. Preprocess
    object$sample_id <- sample_names
    object %<>% add_coldata(verbose = verbose)
    object %<>% preprocess_counts(formula = !!formula,
                    filter_count = filter_count, verbose = verbose, plot = plot)
# Contrast
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Plot
    if (plot)  plot_samples(object)
# Return
    object
}


#' Read RNAseq counts
#'
#' Read tsv file with rnaseq counts into SummarizedExperiment
#'
#' File format: header row
#'              feature annotations in first few columns
#'              feature counts      in next columns
#'
#' @param file      string: path to rnaseq counts file
#' @param fid_col   number of name of column with feature identifiers
#' @param fname_col string or number: feature name variable
#' @param filter_count number (default 10): filter out features
#' with less than 10 counts (in the smallest library) across samples. Filtering
#' performed with \code{\link[edgeR]{filterByExpr}}
#' @param colfile       coldata file (string) or NULL
#' @param by.x          merge var in file
#' @param by.y          merge var in colfile
#' @param formula       formula to create design matrix (using svars)
#' @param contrastdefs  contrastdef vector/matrix/list
#' @param verbose       TRUE/FALSE
#' @param plot          TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' read_counts(file)
#'
#' file <- download_data('billing19.rnacounts.txt')
#' read_counts(file)
#' @seealso merge_coldata, merge_rowdata
#' @export
read_counts <- function(file, fid_col = 1,
    fname_col = character(0), filter_count = 10,
    colfile = NULL, by.x = 'sample_id', by.y = 'sample_id',
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object), verbose = TRUE, plot = TRUE){
# Read
    assert_all_are_existing_files(file)
    formula      <- enexpr(formula)
    contrastdefs    <- enexpr(contrastdefs)
    dt <- fread(file, integer64='numeric')
    if (is.character(fid_col)){
        assert_is_subset(fid_col, names(dt))
        fid_col <- which(names(dt)==fid_col)
    }
    expr_cols  <- which(unname(vapply(dt, is.integer, logical(1))))
    fdata_cols <- c(fid_col, 1 + which(unname(!vapply(
                                        dt[, -fid_col, with = FALSE],
                                        is.integer,
                                        logical(1)))))
    object <- read_omics(file,
                        fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                        sid_rows   = 1,            sid_cols   = expr_cols,
                        expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                        fvar_rows  = 1,            fvar_cols  = fdata_cols,
                        fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                        transpose  = FALSE, verbose    = TRUE)
    if ('gene_name' %in% fvars(object)) fvars(object) %<>%
        stri_replace_first_fixed('gene_name', 'feature_name')
    metadata(object)$platform <- 'rnaseq'
    counts(object) <- exprs(object)
    assays(object)$exprs <- NULL
# Prepare
    object %<>% add_coldata(colfile = colfile, by.x = by.x, by.y = by.y)
    object %<>% preprocess_counts(formula = !!formula,
                    filter_count = filter_count, plot = plot, verbose = verbose)
    if (length(fname_col)>0){
        assert_is_subset(fname_col, fvars(object))
        fdata(object)$feature_name <- fdata(object)[[fname_col]]
        fdata(object) %<>% pull_columns(
            intersect(c('feature_id', 'feature_name'), names(.)))
    }
# Contrast
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Plot
    if (plot)  plot_samples(object)
# Return
    object
}


#=========================================================================

#' @title Get/Set counts
#' @description Get / Set counts matrix
#' @param object SummarizedExperiment
#' @param value count matrix (features x samples)
#' @return count matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_counts(file)
#' counts(object) <- exprs(object)
#' counts(object)[1:3, 1:3]
#' @rdname counts
#' @export
setGeneric('counts',   function(object)   standardGeneric("counts"))

#' @rdname counts
setMethod(
    "counts",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$counts)

#' @rdname counts
#' @export
setGeneric('counts<-',   function(object, value)   standardGeneric("counts<-"))

#' @rdname counts
setReplaceMethod(
    "counts",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$counts <- value
        object })

#' @rdname counts
setReplaceMethod(
    "counts",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$counts[] <- value
        object })




#=========================================================================
#' @title Get/Set weights
#' @description Get/Set weight matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @param ... addtional params
#' @return weight matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' weights(object)[1:3, 1:2]
#' weights(object) <- 1; weights(object)[1:3, 1:2]
#' @rdname weights
#' @export
setGeneric('weights', function(object)   standardGeneric("weights"))

#' @rdname weights
setMethod(
    "weights",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$weights)


#' @rdname weights
#' @export
setGeneric('weights<-', function(object, value) standardGeneric("weights<-"))

#' @rdname weights
setReplaceMethod(
    "weights",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$weights <- value
        object })

#' @rdname weights
setReplaceMethod(
    "weights",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        if (!'weights' %in% names(assays(object))){
            assays(object)$weights <- matrix(
            1, nrow=nrow(object), ncol=ncol(object), dimnames=dimnames(object))
        }
        assays(object)$weights[] <- value
        object })

