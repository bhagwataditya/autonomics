#=========================================================================
#
#                           Getters/Setters
#
#=========================================================================


#=========

#' @title Get/Set counts
#' @description Get / Set counts matrix
#' @param object SummarizedExperiment
#' @param value count matrix (features x samples)
#' @return count matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' counts(object) <- exprs(object)
#' counts(object)[1:3, 1:3]
#' @rdname counts
#' @export
setGeneric('counts',   function(object)   standardGeneric("counts"))

#' @rdname counts
setMethod("counts", signature("SummarizedExperiment"),
function(object)   assays(object)$counts)

#' @rdname counts
#' @export
setGeneric('counts<-',   function(object, value)   standardGeneric("counts<-"))

#' @rdname counts
setReplaceMethod("counts", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$counts <- value
    object })

#' @rdname counts
setReplaceMethod("counts", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$counts[] <- value
    object })

#' @rdname counts
setReplaceMethod("counts", signature("SummarizedExperiment", "NULL"),
function(object, value){
    assays(object)$counts <- NULL
    object })


#===============

#' @title Get/Set log2counts
#' @description Get / Set log2counts matrix
#' @param object SummarizedExperiment
#' @param value log2count matrix (features x samples)
#' @return log2count matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2counts(object) <- exprs(object)
#' log2counts(object)[1:3, 1:3]
#' @rdname log2counts
#' @export
setGeneric('log2counts',   function(object)   standardGeneric("log2counts"))

#' @rdname log2counts
setMethod("log2counts", signature("SummarizedExperiment"),
function(object)   assays(object)$log2counts)

#' @rdname log2counts
#' @export
setGeneric('log2counts<-',
function(object, value) standardGeneric("log2counts<-"))

#' @rdname log2counts
setReplaceMethod("log2counts", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2counts <- value
    object })

#' @rdname log2counts
setReplaceMethod("log2counts", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2counts[] <- value
    object })


#===============

#' @title Get/Set log2countsratios
#' @description Get / Set log2countsratios matrix
#' @param object SummarizedExperiment
#' @param value log2countsratios matrix (features x samples)
#' @return log2countsratios matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2countsratios(object) <- exprs(object)
#' log2countsratios(object)[1:3, 1:3]
#' @rdname log2countsratios
#' @export
setGeneric('log2countsratios',   
function(object)   standardGeneric("log2countsratios"))

#' @rdname log2countsratios
setMethod("log2countsratios", signature("SummarizedExperiment"),
function(object)   assays(object)$log2countsratios)

#' @rdname log2countsratios
#' @export
setGeneric('log2countsratios<-',
function(object, value) standardGeneric("log2countsratios<-"))

#' @rdname log2countsratios
setReplaceMethod("log2countsratios", 
signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2countsratios <- value
    object })

#' @rdname log2countsratios
setReplaceMethod("log2countsratios", 
signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2countsratios[] <- value
    object })


#================

#' @title Get/Set cpm
#' @description Get / Set cpm matrix
#' @param object SummarizedExperiment
#' @param value cpm matrix (features x samples)
#' @return cpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' cpm(object) <- exprs(object)
#' cpm(object)[1:3, 1:3]
#' @rdname cpm
#' @export
setGeneric('cpm',   function(object)   standardGeneric("cpm"))

#' @rdname cpm
setMethod("cpm", signature("SummarizedExperiment"),
function(object)   assays(object)$cpm)

#' @rdname cpm
#' @export
setGeneric('cpm<-',   function(object, value)   standardGeneric("cpm<-"))

#' @rdname cpm
setReplaceMethod("cpm", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$cpm <- value
    object })

#' @rdname cpm
setReplaceMethod("cpm", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$cpm[] <- value
    object })


#================

#' @title Get/Set log2cpm
#' @description Get / Set log2cpm matrix
#' @param object SummarizedExperiment
#' @param value log2cpm matrix (features x samples)
#' @return log2cpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2cpm(object) <- exprs(object)
#' log2cpm(object)[1:3, 1:3]
#' @rdname log2cpm
#' @export
setGeneric('log2cpm',   function(object)   standardGeneric("log2cpm"))

#' @rdname log2cpm
setMethod("log2cpm", signature("SummarizedExperiment"),
function(object)   assays(object)$log2cpm)

#' @rdname log2cpm
#' @export
setGeneric('log2cpm<-',  function(object, value)  standardGeneric("log2cpm<-"))

#' @rdname log2cpm
setReplaceMethod("log2cpm", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2cpm <- value
    object })

#' @rdname log2cpm
setReplaceMethod("log2cpm", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2cpm[] <- value
    object })


#================

#' @title Get/Set log2cpmratios
#' @description Get / Set log2cpmratios matrix
#' @param object SummarizedExperiment
#' @param value log2cpmratios matrix (features x samples)
#' @return log2cpmratios matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2cpmratios(object) <- exprs(object)
#' log2cpmratios(object)[1:3, 1:3]
#' @rdname log2cpmratios
#' @export
setGeneric('log2cpmratios',   
function(object)   standardGeneric("log2cpmratios"))

#' @rdname log2cpmratios
setMethod("log2cpmratios", signature("SummarizedExperiment"),
function(object)   assays(object)$log2cpmratios)

#' @rdname log2cpmratios
#' @export
setGeneric('log2cpmratios<-',  
function(object, value)  standardGeneric("log2cpmratios<-"))

#' @rdname log2cpmratios
setReplaceMethod("log2cpmratios", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2cpmratios <- value
    object })

#' @rdname log2cpmratios
setReplaceMethod("log2cpmratios", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2cpmratios[] <- value
    object })


#================

#' @title Get/Set tpm
#' @description Get / Set tpm matrix
#' @param object SummarizedExperiment
#' @param value tpm matrix (features x samples)
#' @return tpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' tpm(object) <- exprs(object)
#' tpm(object)[1:3, 1:3]
#' @rdname tpm
#' @export
setGeneric('tpm',   function(object)   standardGeneric("tpm"))

#' @rdname tpm
setMethod("tpm", signature("SummarizedExperiment"),
function(object)   assays(object)$tpm)

#' @rdname tpm
#' @export
setGeneric('tpm<-',   function(object, value)   standardGeneric("tpm<-"))

#' @rdname tpm
setReplaceMethod("tpm", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$tpm <- value
    object })

#' @rdname tpm
setReplaceMethod("tpm", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$tpm[] <- value
    object })


#===============

#' @title Get/Set log2tpm
#' @description Get / Set log2tpm matrix
#' @param object SummarizedExperiment
#' @param value log2tpm matrix (features x samples)
#' @return log2tpm matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2tpm(object) <- exprs(object)
#' log2tpm(object)[1:3, 1:3]
#' @rdname log2tpm
#' @export
setGeneric('log2tpm',   function(object)   standardGeneric("log2tpm"))

#' @rdname log2tpm
setMethod("log2tpm", signature("SummarizedExperiment"),
function(object)   assays(object)$log2tpm)

#' @rdname log2tpm
#' @export
setGeneric('log2tpm<-',  function(object, value)  standardGeneric("log2tpm<-"))

#' @rdname log2tpm
setReplaceMethod("log2tpm", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2tpm <- value
    object })

#' @rdname log2tpm
setReplaceMethod("log2tpm", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2tpm[] <- value
    object })


#===============

#' @title Get/Set log2tpmratios
#' @description Get / Set log2tpmratios matrix
#' @param object SummarizedExperiment
#' @param value log2tpmratios matrix (features x samples)
#' @return log2tpmratios matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' log2tpmratios(object) <- exprs(object)
#' log2tpmratios(object)[1:3, 1:3]
#' @rdname log2tpmratios
#' @export
setGeneric('log2tpmratios',   
function(object)   standardGeneric("log2tpmratios"))

#' @rdname log2tpmratios
setMethod("log2tpmratios", signature("SummarizedExperiment"),
function(object)   assays(object)$log2tpmratios)

#' @rdname log2tpmratios
#' @export
setGeneric('log2tpmratios<-',  
function(object, value)  standardGeneric("log2tpmratios<-"))

#' @rdname log2tpmratios
setReplaceMethod("log2tpmratios", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2tpmratios <- value
    object })

#' @rdname log2tpmratios
setReplaceMethod("log2tpmratios", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2tpmratios[] <- value
    object })


#========

#' @title Get/Set weights
#' @description Get/Set weight matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @param ... addtional params
#' @return weight matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' weights(object)[1:3, 1:2]
#' weights(object) <- 1; weights(object)[1:3, 1:2]
#' @rdname weights
#' @export
setGeneric('weights', function(object)   standardGeneric("weights"))

#' @rdname weights
setMethod("weights", signature("SummarizedExperiment"),
function(object)   assays(object)$weights)


#' @rdname weights
#' @export
setGeneric('weights<-', function(object, value) standardGeneric("weights<-"))

#' @rdname weights
setReplaceMethod("weights", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$weights <- value
    object })

#' @rdname weights
setReplaceMethod("weights", signature("SummarizedExperiment", "numeric"),
function(object, value){
    if (!'weights' %in% names(assays(object))){
        assays(object)$weights <- matrix(
        1, nrow=nrow(object), ncol=ncol(object), dimnames=dimnames(object))
    }
    assays(object)$weights[] <- value
    object })

#' @rdname weights
setReplaceMethod("weights", signature("SummarizedExperiment", "NULL"),
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



#=============================================================================
#
#               counts2cpm/cpm2counts
#                   scaledlibsizes
#
#=============================================================================


#' Get tmm-scaled libsizes
#' @param counts  counts matri
#' @return scaled libsize vector
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, cpm=FALSE, log2=FALSE, plot=FALSE)
#' scaledlibsizes(counts(object))
#' @export
scaledlibsizes <- function(counts){
    colSums(counts) * edgeR::calcNormFactors(counts)
}


#' Convert between counts and cpm
#' @param x         count/cpm matrix
#' @param libsize  (scaled) libsize vector
#' @return cpm/tpm/count matrix
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, cpm=FALSE, log2=FALSE, plot=FALSE)
#' libsize <- scaledlibsizes(exprs(object))
#' tpm <- counts2tpm(counts(object), genesize = 1)
#' cpm <- counts2cpm(counts(object), libsize)
#' counts  <- cpm2counts(cpm, libsize)
#' sum(counts(object) - counts)
#' @export
counts2cpm <- function(x, libsize = scaledlibsizes(x)){
    t(t(x)/(libsize + 1)) * 1e6
}


#' @rdname counts2cpm
#' @export
cpm2counts <- function(x, libsize){
    1e-6 * t(t(x) * (libsize + 1))
}

#' counts to tpm
#' @param x count matrix
#' @param genesize  genesize vector (kilobase)
#' @return tpm matrix
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, cpm=FALSE, log2=FALSE, plot=FALSE)
#' counts2tpm(counts(object), genesize=1)[1:3, 1:3]
#' @export
counts2tpm <- function(x, genesize){
    x  %<>% '/'(genesize)
    libsize <- matrixStats::colSums2(x)
    t(t(x)/libsize) * 1e6
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
    libsize <- scaledlibsizes(counts(object))
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
    fitted.log2count <-(2^fitted.log2cpm) %>% cpm2counts(libsize) %>% log2()
# Fit mean-variance trend
    mean.log2count <- fit$Amean + mean(log2(libsize + 1)) - log2(1e+06)
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
#' @param object       SummarizedExperiment
#' @param formula      designmat formula
#' @param block        blocK svar
#' @param min_count    min count required in some samples
#' @param pseudocount  added pseudocount to avoid log(x)=-Inf
#' @param genesize     genesize fvar to compute tpm
#' @param cpm          whether to compute counts per million (scaled) reads
#' @param tmm          whether to tmm normalize
#' @param voom         whether to voom weight
#' @param log2         whether to log2
#' @param verbose      whether to msg
#' @param plot         whether to plot
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('billing19.rnacounts.txt')
#' object <- .read_rnaseq_counts(file)
#' object$subgroup
#' object %<>% preprocess_rnaseq_counts()
#' @export
preprocess_rnaseq_counts <- function(object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, 'limma'), block = NULL,
    min_count = 10, pseudocount = 0.5, genesize = NULL, cpm  = TRUE, tmm = cpm,
    voom = TRUE, log2 = TRUE, verbose = TRUE, plot = TRUE){
# Initialize
    if (is.null(subgroupvar))  subgroupvar <- default_subgroupvar(object)
    if (is.null(formula))      formula <- default_formula(
                                            object, subgroupvar, fit='limma')
# Filter
    if (verbose) message('\t\tPreprocess')
    design <- create_design(object, formula=formula, verbose = verbose)
    object$libsize <- matrixStats::colSums2(counts(object))
    idx <- filterByExpr(counts(object), design = design,#group=object$subgroup,
                lib.size  = object$libsize, min.count = min_count)
    if (verbose) message('\t\tKeep ', sum(idx), '/', length(idx),
            ' features: count >= ', min_count, ' in at least some samples')
    object %<>% extract(idx, )
# Add pseudocount
    if (pseudocount>0){ if (verbose)  message('\t\t\tpseudocount ', pseudocount)
                        counts(object) %<>% add(pseudocount) }
# Tpm/Cpm normalize
    if (!is.null(genesize)){
        assert_is_subset(genesize, fvars(object))
        if (verbose)  message('\t\t\ttpm')
        tpm(object) <- counts2tpm(counts(object), fdata(object)[[genesize]])}
    if (tmm){   if (verbose)  message('\t\t\tcpm:    tmm scale libsizes')
                object$libsize <- scaledlibsizes(counts(object)) }
    if (cpm){   if (verbose)  message('\t\t\t\tcpm')
                cpm(object) <- counts2cpm(counts(object), object$libsize)
                other <- setdiff(assayNames(object), 'cpm')
                assays(object) %<>% extract(c('cpm', other)) }
# Voom  weight (counts) & dupcor (log2(cpm))
    if (voom){
        if (verbose)  message('\t\t\tvoom:   voom')
        object %<>% add_voom(design, verbose=FALSE, plot=plot & is.null(block))
        if (!is.null(block)){
            object %<>%
                add_voom(design, block=block, verbose=verbose, plot=plot) }}
# Log2 transform
    if (log2){  if (verbose)  message('\t\t\tlog2')
                selectedassays <- c('counts','cpm','tpm')
                selectedassays %<>% intersect(assayNames(object))
                for (curassay in selectedassays){
                    i <- match(curassay, assayNames(object))
                    assays(object)[[i]] %<>% log2()
                    assayNames(object)[[i]] %<>% paste0('log2', .)}}
# Return
    object
}


add_voom <- function(
    object, design, block = NULL, verbose = TRUE, plot = TRUE
){
# Retain samples with subgroups
    n0 <- ncol(object)
    object %<>% extract(, rownames(design))
    if (verbose & nrow(object)<n0)  message('\t\t\t\tRetain ',
            ncol(object), '/', n0, ' samples with subgroup definitions')
# Prepare block
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  cmessage('\t\t\t\tdupcor `%s`', blockvar)
            metadata(object)$dupcor <- duplicateCorrelation(
                log2(cpm(object)), design=design, block=block
            )$consensus.correlation }}
# Run voom
    txt <- if (is.null(block)) '\t\t\t\tvoom' else paste0(
                                '\t\t\t\t`', blockvar, '`-blocked voom')
    if (verbose) message(txt)
    weights <- voom(counts(object),
                    design      = design,
                    lib.size    = object$libsize,
                    block       = block,
                    correlation = metadata(object)$dupcor,
                    plot        = plot)$weights
# Add/Return
    rownames(weights) <- rownames(object)
    colnames(weights) <- colnames(object)
    weights(object) <- weights
    object
}



#==============================================================================
#
#                        .read_rnaseq_counts
#                        .read_rnaseq_bams
#
#==============================================================================

#' @rdname read_rnaseq_counts
#' @export
.read_rnaseq_bams <- function(
    dir, paired, genome, nthreads = detectCores(),
    sfile = NULL, sfileby = NULL, subgroupvar = NULL,
    ffile = NULL, ffileby = NULL, fnamevar    = NULL, verbose = TRUE
){
# Assert
    assert_all_are_existing_files(dir)
    assert_is_a_bool(paired)
    assert_is_a_string(genome)
    assert_is_a_number(nthreads)
# Count reads
    files <- list.files(dir, pattern = ".sam$|.bam$", full.names = TRUE,
                        recursive = TRUE)
    fcounts <- count_reads(files, paired, nthreads=nthreads, genome=genome)
# Forge SummarizedExperiment
    filenames   <- basename(file_path_sans_ext(files))
    subdirnames <- basename(dirname(files))
    sample_names <- if (has_no_duplicates(filenames)){           filenames
                    } else if (has_no_duplicates(subdirnames)){  subdirnames
                    } else {           paste0(subdirnames, '_', filenames)}

    colnames(fcounts$counts) <- sample_names
    object <- matrix2sumexp(fcounts$counts,
                            fdt = fcounts$annotation, fdtby='feature_id')
    assayNames(object)[1] <- 'counts'
# Add sample/feature data
    object %<>% merge_sfile(sfile, by.x = 'sample_id',  by.y = sfileby,
                            subgroupvar = subgroupvar, verbose = verbose)
    object %<>% merge_ffile(ffile, by.x = 'feature_id', by.y = ffileby,
                            fnamevar = fnamevar, verbose = verbose)
    metadata(object)$platform <- 'rnaseq'
    metadata(object)$file <- dir
# Return
    object
}


#' @rdname read_rnaseq_counts
#' @export
.read_rnaseq_counts <- function(file, fid_col = 1,
    sfile = NULL, sfileby  = NULL, ffile = NULL, ffileby = NULL, verbose = TRUE
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
    object <- matrix2sumexp(
                counts1, fdt = fdata1, fdtby = fid_col, verbose = verbose)
    assayNames(object)[1] <- 'counts'
# sumexp
    object %<>% merge_sfile(sfile = sfile, by.x = 'sample_id',by.y = sfileby)
    object %<>% merge_ffile(ffile = ffile, by.x='feature_id', by.y = ffileby)
    metadata(object)$platform <- 'rnaseq'
    object$subgroup %<>% factor()
    levels(object$subgroup) %<>% make.names()
    object
}


#=============================================================================
#
#                           read_rnaseq_bams
#                           read_rnaseq_counts
#
#=============================================================================


#' @rdname read_rnaseq_counts
#' @export
read_rnaseq_bams <- function(
    dir, paired, genome, nthreads = detectCores(),
    sfile = NULL, sfileby = NULL, subgroupvar = NULL, block = NULL,
    ffile = NULL, ffileby = NULL, fnamevar = NULL,
    formula = NULL, min_count = 10, pseudocount = 0.5, genesize = NULL,
    cpm = TRUE, tmm = cpm, log2 = TRUE, pca = FALSE, fit = NULL, 
    voom = !is.null(fit), contrastdefs = NULL, verbose = TRUE, plot=TRUE
){
# Read
    object <- .read_rnaseq_bams(dir   = dir,
                                paired   = paired,
                                genome   = genome,
                                nthreads = nthreads,
                                sfile    = sfile,
                                sfileby  = sfileby,
                                ffile    = ffile,
                                ffileby  = ffileby,
                                fnamevar = fnamevar,
                                subgroupvar = subgroupvar)
# Preprocess
    object %<>% preprocess_rnaseq_counts(formula    = formula,
                                        block       = block,
                                        min_count   = min_count,
                                        pseudocount = pseudocount,
                                        genesize    = genesize,
                                        cpm         = cpm,
                                        tmm         = tmm,
                                        voom        = voom,
                                        log2        = log2,
                                        verbose     = verbose,
                                        plot        = plot)
# Analyze
    object %<>% analyze(pca=pca, fit=fit, formula = formula, block = block, 
                    contrastdefs = contrastdefs, verbose = verbose, plot = plot)
# Return
    object
}



#' Read rnaseq
#'
#' Read/analyze rnaseq counts / bamfiles
#'
#' @param dir   read_rnaseq_bams: bam/samfile dir
#' @param paired   read_rnaseq_bams: whether paired end reads
#' @param genome   read_rnaseq_bams: mm10"/"hg38"/etc. or GTF file
#' @param nthreads read_rnaseq_bams: nthreads used by Rsubread::featureCounts()
#' @param file     read_rnaseq_counts: count file
#' @param fid_col  featureid fvar
#' @param sfile    sample file
#' @param sfileby  sample file mergeby column
#' @param subgroupvar  subgroup svar
#' @param block    block svar
#' @param ffile    feature file
#' @param ffileby  feature file mergeby column
#' @param fnamevar featurename fvar
#' @param formula  designmat formula
#' @param contrastdefs contrastdef vector/matrix/list
#' @param min_count    min feature count required in some samples
#' @param pseudocount  added pseudocount to prevent -Inf log2 values
#' @param genesize     genesize fvar for tpm
#' @param tmm      whether to tmm-scale library sizes
#' @param cpm      whether to compute cpm
#' @param voom     whether to compute voom precision weights
#' @param log2     whether to log2 transform
#' @param pca      whether to pca
#' @param fit      fit model: NULL, 'limma', 'lm', 'lme', 'lmer', 'wilcoxon'
#' @param verbose  whether to message
#' @param plot     whether to plot
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, pca= TRUE, fit='limma')
#' @author Aditya Bhagwat, Shahina Hayat
#' @export
read_rnaseq_counts <- function(
    file, fid_col = 1,
    sfile = NULL, sfileby = NULL, subgroupvar = NULL, block = NULL,
    ffile = NULL, ffileby = NULL, fnamevar = NULL,
    formula = NULL, min_count = 10, pseudocount = 0.5, genesize = NULL,
    cpm = TRUE, tmm = cpm, log2 = TRUE, pca = FALSE, 
    fit = NULL, voom = !is.null(fit), contrastdefs = NULL, 
    verbose = TRUE, plot = TRUE
){
# Read
    object <- .read_rnaseq_counts(file,
                                fid_col     = fid_col,
                                sfile       = sfile,
                                sfileby     = sfileby,
                                subgroupvar = subgroupvar,
                                ffile       = ffile,
                                ffileby     = ffileby,
                                fnamevar    = fnamevar,
                                verbose     = verbose)
# Preprocess
    object %<>% preprocess_rnaseq_counts(formula    = formula,
                                        block       = block,
                                        min_count   = min_count,
                                        pseudocount = pseudocount,
                                        genesize    = genesize,
                                        cpm         = cpm,
                                        tmm         = tmm,
                                        voom        = voom,
                                        log2        = log2,
                                        verbose     = verbose,
                                        plot        = plot)
# Analyze
    object %<>% analyze(pca=pca, fit=fit, formula = formula, block = block, 
                        weightvar = if (voom) 'weights' else NULL,
                    contrastdefs = contrastdefs, verbose = verbose, plot=plot)
# Return
    object
}


#' Analyze
#' @param object       SummarizedExperiment
#' @param pca          whether to perform pca
#' @param fit          NULL, 'limma', 'lm', 'lme', 'lmer', or 'wilcoxon'
#' @param subgroupvar  subgroup svar
#' @param formula      model formula
#' @param block        block svar
#' @param weightvar    NULL or name of weight matrix in assays(object)
#' @param contrastdefs contrastdefs vector/matrix/list
#' @param verbose      whether to msg
#' @param plot         whether to plot
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' object %<>% analyze(pca=TRUE, fit='limma')
#' @export
analyze <- function(
    object,
    pca = FALSE,
    fit = NULL,
    subgroupvar = default_subgroupvar(object),
    formula = default_formula(object, subgroupvar, fit),
    block = NULL,
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    contrastdefs = contrast_coefs(object, formula),
    verbose = TRUE,
    plot = TRUE
){
    subgroup <- if (is.null(subgroupvar)) quo(NULL) else sym(subgroupvar)
    if (plot){
        grid.draw(
            grid.arrange(arrangeGrob(
                plot_sample_densities(
                    object[, seq_len(min(30, ncol(object)))], 
                    fill = !!subgroup), 
                plot_feature_densities(object[sample(ncol(object), 4)]), 
                ncol=2),
            plot_summarized_detections(
                object, group = !!subgroup, fill = !!subgroup),
            nrow=2))
    }
    if (pca)   object %<>% pca(verbose=verbose, plot=plot, color=!!subgroup)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object,subgroupvar,fit)
        if (is.null(contrastdefs)) contrastdefs<- contrast_coefs(object,formula)
        object %<>% fitfun( subgroupvar  = subgroupvar,
                            formula      = formula,
                            contrastdefs = contrastdefs,
                            block        = block,
                            weightvar    = weightvar,
                            verbose      = verbose,
                            plot         = plot) }
    object
}


