#=========================================================================
#
#                           Getters/Setters
#
#=========================================================================

#' @title Get/Set counts
#' @description Get / Set counts matrix
#' @param object SummarizedExperiment
#' @param value count matrix (features x samples)
#' @return count matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file)
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


#' @title Get/Set cpm
#' @description Get / Set cpm matrix
#' @param object SummarizedExperiment
#' @param value cpm matrix (features x samples)
#' @return cpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file)
#' cpm(object) <- exprs(object)
#' cpm(object)[1:3, 1:3]
#' @rdname cpm
#' @export
setGeneric('cpm',   function(object)   standardGeneric("cpm"))

#' @rdname cpm
setMethod(
    "cpm",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$cpm)

#' @rdname cpm
#' @export
setGeneric('cpm<-',   function(object, value)   standardGeneric("cpm<-"))

#' @rdname cpm
setReplaceMethod(
    "cpm",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$cpm <- value
        object })

#' @rdname cpm
setReplaceMethod(
    "cpm",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$cpm[] <- value
        object })


#' @title Get/Set log2cpm
#' @description Get / Set cpm matrix
#' @param object SummarizedExperiment
#' @param value log2cpm matrix (features x samples)
#' @return log2cpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file)
#' log2cpm(object) <- exprs(object)
#' log2cpm(object)[1:3, 1:3]
#' @rdname log2cpm
#' @export
setGeneric('log2cpm',   function(object)   standardGeneric("log2cpm"))

#' @rdname log2cpm
setMethod(
    "log2cpm",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$log2cpm)

#' @rdname log2cpm
#' @export
setGeneric('log2cpm<-',   function(object, value)   standardGeneric("log2cpm<-"))

#' @rdname log2cpm
setReplaceMethod(
    "log2cpm",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$log2cpm <- value
        object })

#' @rdname log2cpm
setReplaceMethod(
    "log2cpm",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$log2cpm[] <- value
        object })



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

#' @rdname weights
setReplaceMethod(
    "weights",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){
        assays(object)$weights <- NULL
        object })

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

#==============================================================================
#
#                        .read_rnaseq_counts
#                        .read_rnaseq_bams
#
#==============================================================================

#' @rdname read_rnaseq_bams
#' @export
.read_rnaseq_bams <- function(
    bamdir, paired, genome, nthreads = detectCores(),
    samplefile = NULL, sampleidvar = 'sample_id', subgroupvar = 'subgroup',
    verbose = TRUE
){
# Assert
    assert_all_are_existing_files(bamdir)
    assert_is_a_bool(paired)
    assert_is_a_string(genome)
    assert_is_a_number(nthreads)
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
# Add sample/feature data
    message("\t\tAdd fdata")
    fcounts$annotation$feature_id %<>% as.character()
    rowData(object)  <- fcounts$annotation
    object$sample_id <- sample_names
    object %<>% merge_samplefile( samplefile = samplefile,
        sampleidvar = sampleidvar, subgroupvar = subgroupvar, verbose = verbose)
# Return
    object
}

#' @rdname read_rnaseq_counts
#' @export
.read_rnaseq_counts <- function(file, fid_col = 1,
    samplefile  = NULL, sampleidvar  = NULL, subgroupvar = NULL,
    featurefile = NULL, featureidvar = NULL, featurenamevar = NULL,
    verbose = TRUE
){
# scan
    assert_all_are_existing_files(file)
    if (verbose)  message('\tRead ', file)
    dt <- fread(file, integer64 = 'numeric')
    if (is.numeric(fid_col)) fid_col <- names(dt)[fid_col]
    idx <- vapply(dt, is.integer, logical(1))
    fdata1   <- dt[, !idx, with = FALSE]
    counts1  <- as.matrix(dt[,  idx, with = FALSE])
    rownames(counts1) <- fdata1[[fid_col]]
    object <- matrix2sumexp(counts1, featuredata = fdata1, featureidvar=fid_col)
    assayNames(object)[1] <- 'counts'
# sumexp
    object %<>% merge_samplefile(samplefile = samplefile,
                     sampleidvar = sampleidvar, subgroupvar = subgroupvar)
    object %<>% merge_featurefile(featurefile = featurefile,
                     featureidvar = featureidvar, featurenamevar=featurenamevar)
    metadata(object)$platform <- 'rnaseq'
    object$subgroup %<>% make.names()
    object$subgroup %<>% factor()
    object
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
#' object <- read_rnaseq_counts(file)
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
#                  compute_voom_weights
#                      compute_voom_weights
#
#=============================================================================


explicitly_compute_voom_weights <- function(
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
#                   preprocess_rnaseq_counts
#                       filter_low_count_features
#
#==============================================================================

#' Preprocess RNAseq counts
#' @param object     SummarizedExperiment
#' @param formula    formula
#' @param block      blocK var or NULL
#' @param min_count  min count
#' @param tmm        TRUE/FALSE
#' @param cpm        TRUE/FALSE
#' @param voomweight TRUE/FALSE
#' @param verbose    TRUE/FALSE
#' @param plot       TRUE/FALSE
#' @export
preprocess_rnaseq_counts <- function(object, formula, block = NULL,
    min_count = 10, tmm = TRUE, cpm  = TRUE, voomweight = TRUE,
    verbose = TRUE, plot = TRUE
){
# filter
    if (verbose) message('\t\tPreprocess')
    design <- create_design(object, formula)
    object$lib.size <- colSums(counts(object))
    idx <- filterByExpr(counts(object),
                        design    = design, # group = object$subgroup,
                        lib.size  = object$lib.size,
                        min.count = min_count)
    if (verbose) message('\t\tKeep ', sum(idx), '/', length(idx),
            ' features: count >= ', min_count, ' in at least some samples')
    object %<>% extract(idx, )
# tmm
    if (tmm){   if (verbose)  message('\t\t\ttmm lib.size')
                object$lib.size <- scaledlibsizes(counts(object)) }
# cpm
    if (cpm){   if (verbose)  message('\t\t\tcpm')
                cpm(object) <- counts_to_cpm(counts(object), object$lib.size)
                assays(object) %<>% extract(c('cpm', 'counts')) }
# log2
    if (verbose)  message('\t\t\tlog2')
    cpm(object) %<>% log2()
    assayNames(object) %<>% stri_replace_first_fixed('cpm', 'log2cpm')
# voom
    if (verbose)  message('\t\t\tvoom')
    object %<>% add_voom(design, verbose=verbose, plot = plot & !is.null(block))
    if (!is.null(block)){
        if (verbose)  message('\t\t\t\tvoom with block')
        object %<>% add_voom(design, block=block, verbose=verbose, plot=plot)
    }
# Return
    object
}

add_voom <- function(
    object, design, block = NULL, verbose = TRUE, plot = TRUE
){
    n0 <- ncol(object)
    object %<>% extract(, rownames(design))
    if (verbose & nrow(object)<n0)  message('\t\t\t\tRetain ',
            ncol(object), '/', n0, ' samples with subgroup definitions')

    object %<>% add_blockcor(block = block, design=design, verbose=verbose)
    blockvalues <- if (is.null(block)) NULL else sdata(object)[[block]]
    blockcor <- metadata(object)$blockcor
    weights <- voom(counts(object),
                            design      = design,
                            lib.size    = object$lib.size,
                            block       = blockvalues,
                            correlation = blockcor,
                            plot        = plot
                        )$weights
    rownames(weights) <- rownames(object)
    colnames(weights) <- colnames(object)
    weights(object) <- weights
    object
}



#=============================================================================
#
#                           read_rnaseq_bams
#                           read_rnaseq_counts
#
#=============================================================================


#' Read RNAseq SAM/BAM files into SummarizedExperiment
#' @param bamdir        SAM/BAM dir path (string)
#' @param paired        whether reads are paired end (TRUE/FALSE)
#' @param genome        string: either "mm10", "hg38" etc. or a GTF file
#' @param nthreads      no of cores to be used by Rsubread::featureCounts()
#' @param samplefile    sample file
#' @param sampleidvar   sampleid var
#' @param subgroupvar   subgroup var
#' @param formula       formula to create design matrix (using svars)
#' @param contrastdefs  contrast definition vector/matrix/list
#' @param filter_count  min feature count required by at least one sample
#' @param tmm           TRUE/FALSE: tmm-correct library sizes?
#' @param cpm           TRUE/FALSE: compute counts per million?
#' @param voomweight    TRUE/FALSE: compute voom precision weights?
#' @param verbose       TRUE/FALSE
#' @param plot          TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # in-built genome
#' bamdir <- download_data("billing16.bam.zip")
#' # object <- read_rnaseq_bams(bamdir, paired=TRUE, genome="hg38")
#' # object <- .read_rnaseq_bams(bamdir, paired=TRUE, genome="hg38")
#' # gtffile <- download_gtf("Homo sapiens")
#' # object <- read_rnaseq_bams(bamdir, paired=TRUE, genome=gtffile)
#' @author Aditya Bhagwat, Shahina Hayat
#' @export
read_rnaseq_bams <- function(
    bamdir, paired, genome, nthreads = detectCores(),
    samplefile = NULL, sampleidvar = 'sample_id',
    subgroupvar = 'subgroup',
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object), filter_count = 10,
    tmm = TRUE, cpm = TRUE, voomweight = TRUE, verbose = TRUE, plot=TRUE
){
# Read
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
    object <- .read_rnaseq_bams(bamdir      = bamdir,
                                paired      = paired,
                                genome      = genome,
                                nthreads    = nthreads,
                                samplefile  = samplefile,
                                sampleidvar = sampleidvar,
                                subgroupvar = subgroupvar)
# Preprocess/Analyze
    object %<>% preprocess_rnaseq_counts(formula  = !!formula,
                                min_count  = min_count,
                                tmm        = tmm,
                                cpm        = cpm,
                                voomweight = voomweight,
                                verbose    = verbose,
                                plot       = plot)
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Plot/Return
    if (plot)  plot_samples(object)
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
#' @param file          string: path to rnaseq counts file
#' @param fid_col       number of name of column with feature identifiers
#' @param fname_col     string or number: feature name variable
#' with less than 10 counts (in the smallest library) across samples. Filtering
#' performed with \code{\link[edgeR]{filterByExpr}}
#' @param samplefile    sample file
#' @param sampleidvar   sampleid var
#' @param subgroupvar   subgroup var
#' @param formula       formula to create design matrix (using svars)
#' @param contrastdefs  contrastdef vector/matrix/list
#' @param min_count     number (default 10): filter out features
#' @param tmm           TRUE/FALSE: tmm-correct library sizes?
#' @param cpm           TRUE/FALSE: compute counts per million?
#' @param voomweight    TRUE/FALSE: compute voom precision weights?
#' @param verbose       TRUE/FALSE
#' @param plot          TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # BILLING19
#'     # ~ 0 + subgroup
#'     file <- download_data('billing19.rnacounts.txt')
#'     object <- read_rnaseq_counts(file)
#'     extract_limma_summary(object)
#'
#'     # ~ 0 + subgroups | weights
#'     weights(object) <- NULL
#'     object %<>% add_limma(plot=FALSE)
#'     extract_limma_summary(object)
#'
#'# GSE161731
#'    # Download
#'    require(magrittr)
#'    require(GEOquery)
#'    basedir <- '~/autonomicscache/datasets'
#'    subdir  <- '~/autonomicscache/datasets/GSE161731'
#'    if (!dir.exists(subdir))  getGEOSuppFiles("GSE161731",baseDir = basedir)
#'    file       <- paste0(subdir,'/GSE161731_counts.csv.gz')
#'    samplefile <- paste0(subdir,'/GSE161731_counts_key.csv.gz')
#'
#'    # Read
#'    .read_rnaseq_counts(file)
#'    .read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id')
#'    .read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id',
#'                         subgroupvar = 'cohort')
#'
#'    # Read and analyze
#'    read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id')
#'    read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id',
#'                        subgroupvar = 'cohort')
#'
#'    # Read and analyze with block effect
#'    object <- .read_rnaseq_counts(file, samplefile = samplefile,
#'                          sampleidvar = 'rna_id', subgroupvar = 'cohort')
#'    object %<>% filter_samples(
#'          subject_id %in% names(which(table(subject_id)>1)), verbose=TRUE)
#'    object %<>% pca(plot=TRUE)
#'    object %<>% filter_samples(subgroup != 'healthy')
#'    object %<>% pca(plot=TRUE)
#'    tail(fnames(object)) == tail(fdata(object)$feature_id)  # NEEDS FIX!
#'    preprocess_rnaseq_counts(object, subgroupvar = 'cohort', block = 'subject_id')
#'
#'    read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id', subgroupvar = 'cohort', block = 'subject_id')
#'    read_rnaseq_counts(file, samplefile = samplefile, sampleidvar = 'rna_id', subgroupvar = 'gender')
#' @export
read_rnaseq_counts <- function(
    file, fid_col = 1,
    samplefile = NULL, sampleidvar = NULL, subgroupvar = NULL, block = NULL,
    featurefile = NULL, featureidvar = NULL, featurenamevar = NULL,
    formula = if (is.null(subgroupvar)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object), min_count = 10,
    tmm = TRUE, cpm = TRUE, voomweight = TRUE, verbose = TRUE, plot=TRUE
){
# Initialize
    contrastdefs <- enexpr(contrastdefs)
# Read
    object <- .read_rnaseq_counts(
                    file, fid_col = fid_col, samplefile = samplefile,
                    sampleidvar = sampleidvar, subgroupvar = subgroupvar)
    object %<>% preprocess_rnaseq_counts(formula = formula, block = block,
                    min_count = min_count, tmm = tmm, cpm = cpm,
                    voomweight = voomweight, plot = plot, verbose = verbose)
# Contrast
    object %<>% pca()
    object %<>% add_limma(formula = formula, block = block,
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Plot
    if (plot)  plot_samples(object)
# Return
    object
}


