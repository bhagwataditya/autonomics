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
#' counts(object) <- values(object)
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
#' log2counts(object) <- values(object)
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
#' log2countsratios(object) <- values(object)
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
#' cpm(object) <- values(object)
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
#' log2cpm(object) <- values(object)
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
#' log2cpmratios(object) <- values(object)
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
#' tpm(object) <- values(object)
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
#' log2tpm(object) <- values(object)
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
#' log2tpmratios(object) <- values(object)
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
#' object <- read_maxquant_proteingroups(file, plot = FALSE)
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
    gtffile = sprintf("%s/gtf/%s",
        R_user_dir('autonomics', 'cache'),
        basename(make_gtf_url(organism, release) %>% substr(1, nchar(.)-3)))
){
    assert_is_subset(organism,
                    c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'))
    . <- NULL
    remote <- make_gtf_url(organism, release)

    if(file.exists(gtffile)){
        message("GTF file already available at ", gtffile)
    } else {
        message("download   " ,remote)
        message("to         ", gtffile )
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
    if (!requireNamespace('GenomicRanges', quietly = TRUE)){
        message("BiocManager::install('GenomicRanges'). Then re-run.")
        return(object) 
    }
    if (!requireNamespace('rtracklayer', quietly = TRUE)){
        message("BiocManager::install('rtracklayer'). Then re-run.")
        return(object) 
    }
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
# Assert
    if (!requireNamespace('Rsubread', quietly = TRUE)){
        stop("BiocManager::install('Rsubread'). Then re-run.") }
# Common args
    . <- NULL
    args <- list(files = files, isPaired = paired, nthreads = nthreads)
# Inbuilt genome
    if (genome %in% c('mm10', 'mm9', 'hg38', 'hg19')){
        args %<>% c(list(annot.inbuilt = genome))
        fcounts <- do.call(Rsubread::featureCounts, args)
        fcounts$annotation$gene_name <-
            entrezg_to_symbol(fcounts$annotation$GeneID, genome)
# User GTF
    } else {
        assert_all_are_existing_files(genome)
        args %<>% c(list(annot.ext = genome, isGTFAnnotationFile = TRUE,
                        GTF.attrType.extra  = 'gene_name'))
        fcounts <- do.call(Rsubread::featureCounts, args)
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
#' libsize <- scaledlibsizes(values(object))
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
    x  %<>% divide_by(genesize)
    libsize <- matrixStats::colSums2(x, na.rm = TRUE)
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
    log2cpm  <- values(object)
    libsize <- scaledlibsizes(counts(object))
    design   <- create_design(object, formula = !!enquo(formula))
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
#' @param subgroupvar  subgroup svar
#' @param formula      designmat formula
#' @param block        blocK svar
#' @param min_count    min count required in some samples
#' @param pseudo       added pseudocount to avoid log(x)=-Inf
#' @param tpm          TRUE/FALSE : tpm normalize?
#' @param genesizevar  string     : tpm normalization genesize fvar 
#' @param cpm          TRUE/FALSE : cpm normalize? (counts per million (scaled) reads)
#' @param voom         TRUE/FALSE : voom weight?
#' @param log2         TRUE/FALSE : log2 transform?
#' @param verbose      TRUE/FALSE : msg?
#' @param plot         TRUE/FALSE : plot?
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
    formula = default_formula(object, subgroupvar), block = NULL,
    min_count = 10, pseudo = 0.5, tpm  = FALSE, 
    cpm = TRUE, voom = TRUE, log2 = TRUE,
    verbose = TRUE, plot = TRUE
){
# Initialize
    . <- NULL
    if (is.null(subgroupvar))  subgroupvar <- default_subgroupvar(object)
    if (is.null(formula))      formula <- default_formula(object, subgroupvar)
# tpm
    if (verbose) message('\t\tPreprocess')
    if (tpm){   assert_is_subset('genesize', fvars(object))
                if (verbose)  message('\t\t\ttpm')
                tpm(object) <- counts2tpm(counts(object), fdata(object)$genesize) }
# Filter
    object$libsize <- matrixStats::colSums2(counts(object))
    if (min_count>0)  object %<>% filter_by_expr(
                          formula = formula, min_count = min_count, verbose = verbose)
# tmm/cpm/voom normalize
    if (pseudo>0){
        if (verbose)  message('\t\t\tcounts: add ', pseudo)
            counts(object) %<>% add(pseudo) 
        }
    if (cpm){
        if (verbose)  message('\t\t\tcpm:    tmm scale libsizes')
        object$libsize <- scaledlibsizes(counts(object))
        if (verbose)  message('\t\t\t\tcpm')
        cpm(object) <- counts2cpm(counts(object), object$libsize)
        other <- setdiff(assayNames(object), 'cpm')
        assays(object) %<>% extract(c('cpm', other)) 
    }
    if (voom){  
        object %<>% add_voom(
            formula, verbose = verbose, plot = plot & is.null(block))
        if (!is.null(block))  object %<>% add_voom(
            formula, block = block, verbose = verbose, plot = plot) 
    }
    if (pseudo>0){
        if (verbose)  message('\t\t\tcounts: rm ', pseudo)
        counts(object) %<>% subtract(pseudo) 
    }
# Log2 transform
    if (log2){    assays(object)$log2counts <- log2(pseudo + counts(object))
        if (tpm)  assays(object)$log2tpm    <- log2(pseudo + tpm(object))
        if (cpm)  assays(object)$log2cpm    <- log2(pseudo + cpm(object))
    }
# Rm pseudocounts
# Order assays
    ass <- c('log2cpm', 'log2tpm', 'log2counts', 'cpm', 'tpm', 'counts', 'weights')
    ass %<>% intersect(assayNames(object))
    assays(object) %<>% extract(ass)
# Return
    object
}


filter_by_expr <- function(object, formula, min_count, verbose){
    design <- create_design(object, formula = formula, verbose = FALSE)
    idx <- filterByExpr(counts(object), design = design,#group=object$subgroup,
                lib.size  = object$libsize, min.count = min_count)
    if (verbose) message('\t\t\tKeep ', sum(idx), '/', length(idx),
            ' features: count >= ', min_count, ' (', formula2str(formula), ')')
    object %<>% extract(idx, )
    object
}


add_voom <- function(
    object, formula, block = NULL, verbose = TRUE, plot = TRUE
){
# Retain samples with subgroups
    n0 <- ncol(object)
    design <- create_design(object, formula = formula, verbose = FALSE)
    object %<>% extract(, rownames(design))
    if (verbose & nrow(object)<n0)  message('\t\t\t\tRetain ',
            ncol(object), '/', n0, ' samples with subgroup definitions')
# Prepare block
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  message('\t\t\tdupcor `', blockvar, '`')
            metadata(object)$dupcor <- 0.2#duplicateCorrelation(
                #log2(cpm(object)), design=design, block=block
            #)$consensus.correlation 
            }}
# Run voom
    txt <- sprintf('\t\t\tvoom: %s', formula2str(formula))
    if (!is.null(block))  txt %<>% paste0(' | ', blockvar)
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

add_ensdb <- function(object, ensdb, verbose = TRUE){
# Assert
    if (is.null(ensdb)) return(object)
    if (!stri_startswith_fixed(fdt(object)$feature_id[1],'ENS'))  return(object)
# Genesize
    genesize <- ensembldb::lengthOf(
        ensdb, filter = ensembldb::GeneidFilter(fdt(object)$feature_id))
    genesize <- data.table(feature_id = names(genesize), genesize = genesize)
    object %<>% merge_fdata(genesize)
# Return
    object
}
    
    

#' @rdname read_rnaseq_counts
#' @export
.read_rnaseq_bams <- function(
    dir, paired, genome, nthreads = detectCores(),
    sfile = NULL, by.y = NULL, subgroupvar = 'subgroup',
    ensdb = NULL, verbose = TRUE
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
    object <- matrix2sumexp(fcounts$counts)
    object %<>% merge_fdata(data.table(fcounts$annotation), by.y = 'feature_id')
    assayNames(object)[1] <- 'counts'
# Add sample/feature data
    object %<>% merge_sfile(
        sfile, by.x = 'sample_id',  by.y = by.y, verbose = verbose)
    fdt(object)$feature_id %<>% split_extract_fixed('.', 1)  # drop ensemblid version
    object %<>% add_ensdb(ensdb)    
# Return
    object
}

is_numeric_character <- function(x)  all(!is.na(suppressWarnings(as.numeric(x))))

#' @rdname read_rnaseq_counts
#' @export
.read_rnaseq_counts <- function(file, fid_col = 1,
    sfile = NULL, by.y  = NULL, subgroupvar = 'subgroup', ensdb = NULL, 
    verbose = TRUE
){
# counts
    assert_all_are_existing_files(file)
    if (verbose)  message('\tRead ', file)
    dt <- fread(file, integer64 = 'numeric')
    if (is.numeric(fid_col)) fid_col <- names(dt)[fid_col]
    idx <- vapply(dt, is.integer,           logical(1)) |
           vapply(dt, is_numeric_character, logical(1))
    fdt0   <- dt[, !idx, with = FALSE]
    if (stri_startswith_fixed(fdt0[[fid_col]][1], 'ENS')){
        fdt0[[fid_col]] %<>% split_extract_fixed('.', 1)  # drop ensemblid version
    }
    assertive::assert_has_no_duplicates(fdt0[[fid_col]])
    counts1  <- as.matrix(dt[,  idx, with = FALSE])
    rownames(counts1) <- fdt0[[fid_col]]
    object <- matrix2sumexp(counts1, verbose = verbose)
    assayNames(object)[1] <- 'counts'
# fdata
    object %<>% merge_fdata(fdt0, by.y = fid_col, verbose = verbose)
    object %<>% add_ensdb(ensdb)
# sdata
    object %<>% merge_sfile(sfile = sfile, by.x = 'sample_id', by.y = by.y)
    object %<>% add_subgroup(subgroupvar, verbose = verbose)
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
    sfile = NULL, by.y = NULL, subgroupvar = 'subgroup', block = NULL,
    formula = NULL, min_count = 10, pseudo = 0.5, 
    ensdb = NULL, tpm = FALSE, cpm = TRUE, log2 = TRUE, 
    pca = TRUE, fit = 'limma', voom = cpm, contrasts = NULL, 
    verbose = TRUE, plot = pca & !is.null(fit), 
    feature_id = NULL, sample_id = NULL
){
# Read
    object <- .read_rnaseq_bams(
        dir         = dir,          paired      = paired,
        genome      = genome,       nthreads    = nthreads,
        sfile       = sfile,        by.y        = by.y,
        ensdb       = ensdb,        subgroupvar = subgroupvar)
# Preprocess
    object %<>% preprocess_rnaseq_counts(
        subgroupvar = subgroupvar,  formula     = formula,
        block       = block,        min_count   = min_count,
        pseudo      = pseudo,       tpm         = tpm,
        cpm         = cpm,          voom        = voom,
        log2        = log2,         verbose     = verbose,
        plot        = plot)
# Analyze
    object %<>% analyze(
        pca         = pca,          fit         = fit, 
        subgroupvar = subgroupvar,  formula     = formula, 
        block       = block,        contrasts   = contrasts, 
        verbose     = verbose,      plot        = plot, 
        feature_id  = feature_id,   sample_id   = sample_id)
# Return
    object
}


#' Read rnaseq counts/bams
#' @param dir           read_rnaseq_bams: bam/sam dir
#' @param paired        read_rnaseq_bams: TRUE/FALSE : paired end reads ?
#' @param genome        read_rnaseq_bams: 'mm10', 'hg38', etc. or GTF file
#' @param nthreads      read_rnaseq_bams: nthreads used by Rsubread::featureCounts()
#' @param file          count file
#' @param fid_col       featureid column (number or string)
#' @param sfile         sample file
#' @param by.y          sample file mergeby column
#' @param subgroupvar   subgroup svar
#' @param block         model blockvar: string or NULL
#' @param formula       model formula
#' @param min_count     min feature count required in some samples
#' @param pseudo        pseudocount added to prevent -Inf log2 values
#' @param tpm           TRUE/FALSE : add tpm to assays ( counts / libsize / genelength ) ?
#' @param ensdb         EnsDb with genesizes : e.g. AnnotationHub::AnnotationHub[['AH64923']]
#' @param cpm           TRUE/FALSE : add cpm to assays ( counts / effectivelibsize ) ?
#' @param log2          TRUE/FALSE : log2 transform ?
#' @param plot          TRUE/FALSE : plot?
#' @param pca           TRUE/FALSE : perform and plot pca?
#' @param fit           model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param voom          model weights to be computed?  TRUE/FALSE
#' @param coefs         model coefficients          of interest: string vector or NULL
#' @param contrasts     model coefficient contrasts of interest: string vector or NULL
#' @param verbose       TRUE/FALSE : message?
#' @param feature_id    string: feature for summary plot
#' @param sample_id     string: sample  for summary plot
#' @param palette       color palette : named string vector
#' @return SummarizedExperiment
#' @examples
#' # read_rnaseq_bams
#'     if (requireNamespace('Rsubread')){
#'        dir <- download_data('billing16.bam.zip')
#'         object <- read_rnaseq_bams(
#'            dir, paired = TRUE, genome = 'hg38', pca = TRUE, fit = 'limma', plot = TRUE)  
#'     }
#' # read_rnaseq_counts
#'     file <- download_data('billing19.rnacounts.txt')
#'     object <- read_rnaseq_counts(file, fit = 'limma', coefs = 'E05')
#'     object <- read_rnaseq_counts(file, fit = 'limma', coefs = 'E05', voom = FALSE)
#'     object <- read_rnaseq_counts(file, fit = 'limma', coefs = 'E05', voom = FALSE, cpm = FALSE)
#'     object <- read_rnaseq_counts(file, fit = 'limma', coefs = 'E05', voom = FALSE, cpm = FALSE, log2 = FALSE)
#'     
#' # read_rnaseq_counts(tpm = TRUE)
#'     \dontrun{
#'     ah <- AnnotationHub::AnnotationHub()
#'     ensdb <- ah[['AH64923']]
#'     object <- read_rnaseq_counts(file, fit = 'limma', coefs = 'E02', tpm = TRUE, ensdb = ensdb)
#'     }
#' @author Aditya Bhagwat, Shahina Hayat
#' @export
read_rnaseq_counts <- function(
    file, fid_col = 1, sfile = NULL, by.y = NULL, 
    subgroupvar = 'subgroup', block = NULL, formula = NULL, min_count = 10, 
    pseudo = 0.5, tpm = FALSE, ensdb = NULL, cpm = !tpm, log2 = TRUE,
    plot = FALSE, pca = plot, fit = if (plot) 'limma' else NULL, voom = cpm, 
    coefs = NULL, contrasts = NULL, verbose = TRUE, 
    feature_id = NULL, sample_id = NULL, palette = NULL
){
# Read
    object <- .read_rnaseq_counts(
        file,                       fid_col     = fid_col,
        sfile       = sfile,        by.y        = by.y,
        subgroupvar = subgroupvar,  ensdb       = ensdb,
        verbose     = verbose)
# Preprocess
    object %<>% preprocess_rnaseq_counts(
        subgroupvar = subgroupvar,  formula     = formula,
        block       = block,        min_count   = min_count,
        pseudo      = pseudo,       tpm         = tpm,
        cpm         = cpm,          voom        = voom,
        log2        = log2,         verbose     = verbose,
        plot        = FALSE)
# Analyze
    object %<>% analyze(
        pca          = pca,           fit       = fit, 
        subgroupvar  = subgroupvar,   formula   = formula, 
        block        = block,         weightvar = if (voom) 'weights' else NULL,
        coefs        = coefs,         contrasts = contrasts, 
        verbose      = verbose,       plot      = plot, 
        feature_id   = feature_id,    sample_id = sample_id, 
        palette      = palette)
# Return
    object
}


.read_salmon_sample <- function(sampledir){
    file <- sprintf('%s/quant.sf', sampledir)
    dt <- fread(file)[, .(Name, TPM)]
    setnames(dt, 'TPM', basename(sampledir))
    dt
}


#' Read salmon
#' @param dir   salmon results rootdir
#' @param ensdb EnsDb object
#' @param sfile samplefile
#' @param by    samplefile column to merge by
#' @return SummarizedExperiment
#' @examples 
#' # dir <- '../bh/salmon_quants'
#' # sfile <- '../bh/samplesheet.csv'
#' # by <- 'salmonDir'
#' # ah <- AnnotationHub::AnnotationHub()
#' # ensdb <- ah[['AH98078']]
#' # read_salmon(dir, sfile = sfile, by = 'salmonDir', ensdb = ensdb)
#' @export
read_salmon <- function(
    dir, sfile = NULL, by = NULL, ensdb = NULL
){
# Assert
    assertive.files::assert_all_are_dirs(dir)
    if (!is.null(ensdb))   assert_is_all_of(ensdb, 'EnsDb')
    if (!is.null(sfile))   assert_all_are_existing_files(sfile)
    if (!is.null(by)) assert_is_subset(by, names(fread(sfile)))
# Read
    samples <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
    object <- Map(.read_salmon_sample, samples)
    names(object) %<>% basename()
    object %<>% Reduce(function(x, y) merge(x, y, by = 'Name'), .)
    object %<>% dt2mat()
    object <- list(log2tpm = log2(0.00001 + object))  # assay
    object %<>% SummarizedExperiment()
    sdt(object)$sample_id <- snames(object)           # samples
    object %<>% merge_sfile(sfile, by.y = by)
    fdt(object)$feature_id <- fnames(object)          # features
    fdt(object)$enst <- fdt(object)$feature_id %>% split_extract_fixed('.', 1)
    fdt(object)$gene <- ensembldb::mapIds(ensdb, keys = fdt(object)$enst, keytype = 'TXID', column = 'GENENAME')
    object
}
