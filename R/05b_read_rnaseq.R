#=========================================================
# RNASEQ
#=========================================================
#download_gtf
release_to_build <- function(release, organism){
    if        (organism == 'Homo sapiens'){
        if (release >= 76)  'GRCh38'   else 'GRCh37'  }

    else if   (organism == 'Mus musculus'){
        if (release >= 68)  'GRCm38'   else 'NCBIM37'  }

    else if   (organism == 'Rattus norvegicus'){
        if (release >= 80)  'Rnor_6.0' else 'Rnor_5.0'  }
}

#-----------------------------------------------------------
# Following GTF functions are soft-deprecated.
# Better to outsource this functionality to biomartr::getGTF
#-----------------------------------------------------------

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
#' @param release    GTF release. By default release 95 selected
#' @param gtffile    string: path to local GTF file
#' @examples
#' \dontrun{ # requires internet and does not always work:
#'           # https://stackoverflow.com/questions/55532102
#'    download_gtf(organism = 'Homo sapiens')
#'    download_gtf(organism = 'Mus musculus')
#'    download_gtf(organism = 'Rattus norvegicus')
#' }
#' @noRd
download_gtf <- function(
    organism,
    release = 100,
    gtffile = sprintf("~/importomicscache/gtf/%s",
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


#' Read GTF file into data.table
#'
#' Read GTF file into a data.table. Filter for particular values of variable.
#'
#' @param gtffile    string: path to gtffile
#' @param var        string: variable on which to filter
#' @param values     filter variable for these values. If NULL, not filtering.
#' @param writefile  string: file to write gtf table to
#' @examples
#' \dontrun{ # requires internet connection
#'    require(magrittr)
#'    gtffile <- download_gtf(organism = 'Homo sapiens', release = 95)
#'    gtfdt <- read_gtf(gtffile, var = 'gene_id', values = 'ENSG00000198947')
#' }
#' @noRd
read_gtf <- function(
    gtffile,
    var       = 'gene_id',
    values    = NULL,
    writefile = NULL
){

    # Assert
    assert_all_are_existing_files(gtffile)
    assert_is_a_string(var)

    # Read
    dt <- rtracklayer::import(gtffile) %>%
          GenomicRanges::as.data.frame() %>% data.table()

    # Filter
    if (!is.null(values)){
        dt %>% setkeyv(var)
        dt %<>% extract(values)
    }

    # Write
    if (!is.null(writefile)){
        cmessage("\t\tWrite   %s", writefile)
        dir.create(dirname(writefile), recursive = TRUE, showWarnings = FALSE)
        fwrite(dt, writefile)
    }

    # Return
    dt

}


#' Read RNAseq BAM files into SummarizedExperiment
#'
#' @param bamdir       string: path to SAM or BAM file directory
#'                    (one SAM or BAM file per sample)
#' @param ispaired     TRUE or FALSE (default): paired end reads?
#' @param gtffile      NULL (use Rsubread's default human/mouse annotations) or
#'                     string (path to GTF file)
#' @param fvars        character vector: GTF variables to include in object.
#' @param nthreads     number of cores to be used by Rsubread::featureCounts()
#' @param filter_features_min_count  number
#' @param verbose      TRUE(default) / FALSE
#' @param ...          passed to Rsubread::featureCounts
#' @return SummarizedExperiment
#' @examples
#' bamdir <- download_data("stemcells.bam.zip")
#' # read_bam(bamdir, ispaired = TRUE)
#' @export
read_bam <- function(bamdir, ispaired = FALSE, gtffile = NULL,
    fvars = character(0), nthreads   = detectCores(),
    filter_features_min_count = 10, verbose = TRUE, ...
){
# Assert
    assert_all_are_existing_files(bamdir)
    assert_is_a_bool(ispaired)
    if (!is.null(gtffile))   assert_all_are_existing_files(gtffile)
    assert_is_a_number(nthreads)
# Count reads
    files <- list.files(
        bamdir, pattern = ".sam$|.bam$", full.names = TRUE, recursive = TRUE)
    fcounts <-  featureCounts(
                    files               =  files,
                    annot.ext           =  gtffile,
                    isGTFAnnotationFile = !is.null(gtffile),
                    GTF.attrType.extra  =  if(length(fvars)==0) NULL else fvars,
                    isPairedEnd         =  ispaired,
                    nthreads            =  nthreads,
                    ...)
# Forge SummarizedExperiment
    filenames   <- basename(file_path_sans_ext(files))
    subdirnames <- basename(dirname(files))
    sample_names <- if (has_no_duplicates(filenames)){           filenames
                    } else if (has_no_duplicates(subdirnames)){  subdirnames
                    } else {           paste0(subdirnames, '_', filenames)}
    object <- SummarizedExperiment(assays = list(
        counts = fcounts$counts %>% set_colnames(sample_names)))
# Add fdata
    message("\t\tAdd fdata")
    rowData(object) <- fcounts$annotation[ , c('GeneID', fvars), drop = FALSE]
    rownames(object) <- rowData(object)$GeneID
    fvars(object) %<>% stri_replace_first_fixed('GeneID', 'feature_id')
# Add design. Preprocess
    message("\t\tAdd sdata")
    colData(object) <- DataFrame(
                        sample_id = sample_names, row.names = sample_names)
    object$subgroup <- guess_subgroup_values(object, verbose = FALSE)
    object %<>% preprocess_counts(filter_features_min_count)
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
#' @param verbose   TRUE (default) or FALSE
#' @param filter_features_min_count number (default 10): filter out features
#' with less than 10 counts (in the smallest library) across samples. Filtering
#' performed with \code{\link[edgeR]{filterByExpr}}
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('differentiation.rnacounts.txt')
#' read_counts(file)
#' @seealso merge_sdata, merge_fdata
#' @export
read_counts <- function(
    file, fid_col = 1, fname_col = character(0),
    filter_features_min_count = 10, verbose = TRUE
){
# Read
    assert_all_are_existing_files(file)
    dt <- fread(file, integer64='numeric')
    if (is.character(fid_col)){
        assert_is_subset(fid_col, names(dt))
        fid_col <- which(names(dt)==fid_col)
    }
    expr_cols   <- which(unname(vapply(dt, is.integer, logical(1))))
    fdata_cols  <- c(fid_col,
                     1 + which(unname(!vapply(
                                        dt[, -fid_col, with = FALSE],
                                        is.integer,
                                        logical(1)))))
    message('\t\tRead ', "'", basename(file), "'", ' into SummarizedExperiment')
    object <- read_omics(
                file,
                fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                sid_rows   = 1,            sid_cols   = expr_cols,
                expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                fvar_rows  = 1,            fvar_cols  = fdata_cols,
                fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                transpose  = FALSE,
                verbose    = TRUE)
    counts(object) <- exprs(object)
    assays(object)$exprs <- NULL

# Add design. Preprocess
    sdata(object)$subgroup  <- guess_subgroup_values(object, verbose = TRUE)
    object %<>% preprocess_counts(filter_features_min_count)

# Align fdata
    if (length(fname_col)>0){
        assert_is_subset(fname_col, fvars(object))
        fdata(object)$feature_name <- fdata(object)[[fname_col]]
        fdata(object) %<>% pull_columns(c('feature_id', 'feature_name'))
    }

# Return
    object
}


preprocess_counts <- function(object, filter_features_min_count = 10){

    # Filter lowly expressed features
    sdata(object)$libsize <- colSums(counts(object))
    object %<>% filter_low_count_features(
                    filter_features_min_count, verbose = TRUE)
    sdata(object)$libsize.filtered <- colSums(counts(object))
    sdata(object)$libsize.scaled   <- scaledlibsizes(counts(object))

    # Normalize: counts -> cpm
    if (verbose) message('\t\tTMM normalize and log2transform')
    exprs(object) <- counts(object) %>% counts_to_cpm() %>% log2()

    # Compute and add precision weights
    assays(object)$weights <- compute_precision_weights(object, plot = plot)

    # Return
    if (verbose){
        message('\t\tReturn object: counts(object)  = counts')
        message('\t\t               exprs(object)   = log2cpm')
        message('\t\t               weights(object) = voom precision weights (see plot)')
        message('\t\t               object$subgroup = subgroup values')
    }
    object

}


filter_low_count_features <- function(
    object, filter_features_min_count, verbose
){
    idx <- filterByExpr(counts(object),
                        group     = object$subgroup,
                        lib.size  = scaledlibsizes(counts(object)),
                        min.count = filter_features_min_count)

    if (verbose) message('\t\tKeep ', sum(idx), '/', length(idx),
                ' features: count >= ',
                filter_features_min_count, ' in at least some samples')

    subset(object, idx)
}


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
#' file <- download_data('differentiation.rnacounts.txt')
#' object <- read_counts(file, 'gene_id')
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

compute_precision_weights_once <- function(
   object,
   design = create_design_matrix(object),
   plot   = TRUE,
   ...
){
    counts   <- counts(object)
    scaled_libsizes <- scaledlibsizes(counts)

    counts %>%
    voom(design = design, lib.size = scaled_libsizes, plot = plot, ...) %>%
    extract2('weights') %>%
    set_rownames(rownames(object)) %>%
    set_colnames(colnames(object))
}

#' Does object contain subgroup replicates?
#' @param object SummarizedExperiment
#' @return logical
#' @noRd
contains_replicates <- function(object){
    if (!'subgroup' %in% svars(object)) return(FALSE)
    any(duplicated(sdata(object)$subgroup))
}


#' Create design for voom transformation
#' @param object SummarizedExperiment
#' @param verbose logical(1)
#' @return NULL (if no replicates) or design matrix
#' @noRd
create_voom_design <- function(object, verbose = TRUE){

   # Replicates
   if (contains_replicates(object)) return(create_design_matrix(object))

   # No replicates
   message('\t\t\tsubgroup values not replicated: voom(design=NULL)')
   return(NULL)

}


#' Compute voom precision weights
#' @param object  SummarizedExperiment: exprs(.) returns log2cpm, counts(.)
#'                returns raw counts.
#' @param design  design matrix
#' @param plot    TRUE (default) or FALSE
#' @param ...     passed to limma::voom() -> limma::lmFit()
#' @examples
#' file <- download_file('differentiation.rnacounts.txt')
#' object <- read_counts(file)
#' counts(object) <- exprs(object)
#' compute_precision_weights(object)[1:3, 1:3]
#' @noRd
compute_precision_weights <- function(
    object, design = create_voom_design(object), plot = TRUE
){    # Assert & message
    assert_is_not_null(counts(object))

    # Estimate precision weights
    has_block <- has_complete_block_values(object)
    weights   <- compute_precision_weights_once(
                    object, design = design, plot = !has_block)

    # Update precision weights using block correlation
    if (has_block){
        log2cpm <- log2(counts_to_cpm(counts(object)))
        correlation <-  duplicateCorrelation(
                            log2cpm, design = design, block = object$block,
                            weights = weights)
        correlation %<>% extract2('consensus')
      weights <- compute_precision_weights_once(
                    object, design = design, block = object$block,
                    correlation = correlation, plot = TRUE)
   }

   # Return
   dimnames(weights) <- dimnames(object)
   weights
}

explicitly_compute_precision_weights_once <- function(
    object, design = create_voom_design(object), plot = TRUE, ...
){
    # Extract
    log2cpm  <- exprs(object)
    lib.size <- scaledlibsizes(counts(object))

    # Assert
    n <- nrow(log2cpm)
    if (n < 2L) stop("Need at least two genes to fit a mean-variance trend")

    # Fit linear model
    fit <- lmFit(log2cpm, design=design, ...)

    # Predict
    if (is.null(fit$Amean)) fit$Amean <- rowMeans(log2cpm, na.rm = TRUE)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
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


