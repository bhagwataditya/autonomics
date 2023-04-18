
#' Download tcga example
#' @export
download_tcga_example <- function(){
# Assert
    # https://seandavi.github.io/post/2017-03-04-testing-the-genomicdatacommons-package/
    # https://seandavi.github.io/post/2018-03-25-extracting-clinical-information-using-the-genomicdatacommons-package/
    if (!requireNamespace('GenomicDataCommons', quietly = TRUE)){
        message("BiocManager::install('GenomicDataCommons'). Then re-run.")
    }
    if (!requireNamespace('AnnotationHub', quietly = TRUE)){
        message("BiocManager::install('AnnotationHub'). Then re-run.")
    }
    if (!requireNamespace('ensembldb', quietly = TRUE)){
        message("BiocManager::install('AnnotationHub'). Then re-run.")
    }
    file <- tools::R_user_dir('autonomics', 'cache')
    file %<>% file.path('datasets', 'tcga.rna.rds')
    if (file.exists(file)){ return(file)
    } else {
        message('Run the code in this function manually'); 
        return(invisible(NULL))
    }
    access <- ajcc_pathologic_stage <- age_at_index <- NULL
    analysis.workflow_type <- case_id <- cases.project.project_id  <- NULL
    cases.samples.sample_type <- days_to_death <- days_to_last_follow_up <- NULL
    event <- gender <- gene_id <- gene_name <- gene_type <- genesize <- NULL
    icd_10_code <- N <- sample_id <- sample_type <- timetoevent <- type <- NULL
    unstranded <- vital_status <- NULL
# sdt
    sdtTN  <-  GenomicDataCommons::files()
    sdtTN %<>% GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA')
    sdtTN %<>% GenomicDataCommons::filter(                     type == 'gene_expression')
    sdtTN %<>% GenomicDataCommons::filter(   analysis.workflow_type == 'STAR - Counts')
    sdtTN %<>% GenomicDataCommons::filter(                   access == 'open')
    sdtTN  %>% GenomicDataCommons::manifest() %>% nrow()  # 1231
    
    sdtN <- GenomicDataCommons::filter(sdtTN, cases.samples.sample_type == 'Solid Tissue Normal')
    sdtT <- GenomicDataCommons::filter(sdtTN, cases.samples.sample_type == 'Primary Tumor')
    GenomicDataCommons::manifest(sdtN) %>% nrow()  # 226
    GenomicDataCommons::manifest(sdtT) %>% nrow()  # 2222
    sdtN %<>% GenomicDataCommons::expand(c('cases', 'cases.samples'))
    sdtT %<>% GenomicDataCommons::expand(c('cases', 'cases.samples'))
    sdtN %<>% GenomicDataCommons::results_all()
    sdtT %<>% GenomicDataCommons::results_all()
    
    sdtN$sample_id    <- sdtN$cases %>% lapply(extract2, 'samples') %>% lapply(extract2, 1) %>% vapply(extract2, character(1), 'sample_id')   %>% unname()
    sdtT$sample_id    <- sdtT$cases %>% lapply(extract2, 'samples') %>% lapply(extract2, 1) %>% vapply(extract2, character(1), 'sample_id')   %>% unname()
    sdtN$sample_type  <- sdtN$cases %>% lapply(extract2, 'samples') %>% lapply(extract2, 1) %>% vapply(extract2, character(1), 'sample_type') %>% unname()
    sdtT$sample_type  <- sdtT$cases %>% lapply(extract2, 'samples') %>% lapply(extract2, 1) %>% vapply(extract2, character(1), 'sample_type') %>% unname()
    sdtN$primary_site <- sdtN$cases %>% vapply(extract2, character(1), 'primary_site') %>% unname()
    sdtT$primary_site <- sdtT$cases %>% vapply(extract2, character(1), 'primary_site') %>% unname()
    sdtN$disease_type <- sdtN$cases %>% vapply(extract2, character(1), 'disease_type') %>% unname()
    sdtT$disease_type <- sdtT$cases %>% vapply(extract2, character(1), 'disease_type') %>% unname()
    sdtN$case_id      <- sdtN$cases %>% vapply(extract2, character(1), 'case_id')      %>% unname()
    sdtT$case_id      <- sdtT$cases %>% vapply(extract2, character(1), 'case_id')      %>% unname()
    sdtN %<>% extract(c('case_id', 'sample_type', 'file_id'))
    sdtT %<>% extract(c('case_id', 'sample_type', 'file_id'))
    sdtN %<>% as.data.table()
    sdtT %<>% as.data.table()
    sdtN[, sample_type := 'N']
    sdtT[, sample_type := 'T']
    
    common <- intersect(sdtN$case_id, sdtT$case_id)
    sdtN %<>% extract(common, on = 'case_id')   # 113
    sdtT %<>% extract(common, on = 'case_id')   # 119
    sdtT %<>% extract(, .SD[1], by = 'case_id') # 113
    sampledt <- rbind(sdtN, sdtT)
    sampledt %<>% extract(order(case_id))
    sampledt[, sample_id := paste0(split_extract_fixed(case_id, '-', 1), '.', sample_type)]

# clindt
    clindt <- GenomicDataCommons::gdc_clinical(sampledt$case_id)
    clindt %<>% lapply(data.table)
    clindt$demographic %<>% extract(, .(case_id, gender, age_at_index, vital_status, days_to_death))
    clindt$diagnoses   %<>% extract(, .(case_id, icd_10_code, ajcc_pathologic_stage, days_to_last_follow_up))
    clindt <- merge(clindt$demographic, clindt$diagnoses, by = 'case_id')
    clindt[!is.na(days_to_death),          timetoevent := days_to_death]
    clindt[!is.na(days_to_last_follow_up), timetoevent := days_to_last_follow_up]
    clindt[, timetoevent := timetoevent / 365] # days -> years
    clindt[vital_status=='Alive', event := 0]
    clindt[vital_status=='Dead',  event := 1]
    clindt[, c('vital_status', 'days_to_death', 'days_to_last_follow_up') := NULL]

    clindt[,   case_id := as.character(case_id)]
    sampledt[, case_id := as.character(case_id)]
    sampledt %<>% merge(clindt, by = 'case_id') 
    sampledt %<>% pull_columns(c('sample_id', 'sample_type', 'case_id'))
    
# counts
    fnames <- lapply(sampledt$file_id, GenomicDataCommons::gdcdata, progress = FALSE)  # takes a long time :)
    fnames %<>% unlist()
    fnames <- data.table(file_id = names(fnames), file_path = unname(fnames))
    sampledt %<>% merge(fnames, by = 'file_id')
    dofread <- function(filename, sampleid){    dt <- fread(filename)
                                                dt[, sample_id := sampleid]
                                                dt    }
    cnts <- mapply(dofread, sampledt$file_path, sampledt$sample_id, SIMPLIFY = FALSE)
    cnts %<>% data.table::rbindlist()
    cnts %<>% extract(stri_detect_fixed(gene_id, 'ENSG'))
    cnts %<>% extract(gene_type == 'protein_coding')
    
    cnts %<>% extract(, .(sample_id, gene_id, gene_name, counts = unstranded))
    cnts %<>% unique() # stranded
    cnts[, N := .N, by = c('sample_id', 'gene_name')]
    cnts %<>% extract(N==1)
    fdt0 <- unique(cnts[, .(feature_id = gene_name, gene_id)])
    cnts %<>% data.table::dcast(gene_name ~ sample_id, value.var = 'counts')
    cnts %<>% dt2mat()
    rna <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = cnts))
    fdt(rna)$feature_id <- rownames(rna)
    sdt(rna)$sample_id  <- colnames(rna)
    rna %<>% merge_fdt(fdt0)
    sampledt[, c('file_path', 'file_id', 'case_id') := NULL]
    rna %<>% merge_sdt(sampledt)
    sdt(rna)$disease_entity <- 'breast'

# fdt    
    fdt(rna)$ensg <- fdt(rna)$gene_id %>% split_extract_fixed('.', 1)
    ah <- AnnotationHub::AnnotationHub()
    # AnnotationHub::query(ah, 'Homo sapiens', 'Ensembl', 'hg38')
    ensdb <- ah[['AH109336']]
    genesizedt <- ensembldb::lengthOf(ensdb, filter = ensembldb::GeneIdFilter(fdt(rna)$ensg))
    genesizedt <- data.table(ensg = names(genesizedt), genesize = genesizedt)
    rna %<>% merge_fdt(genesizedt, by.x = 'ensg', by.y = 'ensg')
    rna %<>% filter_features(!is.na(genesize))
    rna$case_id <- rna$sample_id
    rna$case_id %<>% split_extract_fixed('.', 1)
    rna %<>% preprocess_rnaseq_counts(formula = ~ sample_type, block = 'case_id', tpm = TRUE, cpm = TRUE, voom = TRUE)
    saveRDS(rna, file = file)
}




empty_survplot <- function(){
    # https://stackoverflow.com/questions/61907987/produce-empty-plot-with-ggsurvplot
    dt  <- data.table(exprlevel = c(rep("low", 10), rep("high", 10)), 
                      value   = c(rnorm(10,mean = 2), rnorm(10,mean = 3)))
    fit <- survival::survfit(survival::Surv(value) ~ exprlevel, data = dt)
    survminer::ggsurvplot(
        fit, 
        data             = dt, 
        surv.median.line = "none", 
        palette          = rep("white", 2), 
        legend           = "none")
}


.dichotomize_exprs <- function(subdt, percentile){
    value <- NULL
    if (all(is.na(subdt$value)))  return(cbind(subdt, exprlevel = 'no values available'))
    subdt %<>% extract(!is.na(value))
    subdt %<>% extract(order(value))
    n <- floor(0.01*percentile*nrow(subdt))
    lowervalue <- subdt$value[n]
    uppervalue <- rev(subdt$value)[n]
    if (length(lowervalue)==0 | lowervalue==uppervalue){
        subdt <- cbind(subdt[0], exprlevel = character(0))
    } else {
        lowergroup <- paste0(signif(lowervalue,1), '-')
        uppergroup <- paste0(signif(uppervalue,1), '+') 
        subdt <- rbind(cbind(subdt[value<=lowervalue], exprlevel = lowergroup),
                       cbind(subdt[value>=uppervalue], exprlevel = uppergroup))
        #subdt$exprlevel %<>% factor(c(lowergroup, uppergroup))
    }
    subdt
}

dichotomize_exprs <- function(dt, percentile){
    dt %<>% extract(, .dichotomize_exprs(.SD, percentile = percentile), by = 'feature_id')
    
}

.fit_survival <- function(subdt, samples = FALSE){
    timetoevent <- event <- exprlevel <- NULL
    diff <- survival::survdiff(survival::Surv(timetoevent, event) ~ exprlevel, data = subdt)
    coef <- suppressWarnings(coef(summary( survival::coxph(
                survival::Surv(subdt$timetoevent, subdt$event)~subdt$value)))[,'coef' ])
    exprlevels <- unique(subdt$exprlevel)
    exprlevels %<>% extract(order(as.numeric(substr(., 1, nchar(.)-1))))
    if (samples){
        lo <- unique(subdt[exprlevel == exprlevels[1]])$sample_id
        hi <- unique(subdt[exprlevel == exprlevels[2]])$sample_id
        lo %<>% as.character() %>% commonify_strings()
        hi %<>% as.character() %>% commonify_strings()
        return(data.table(
                   `effect~surv~LR` = sign(coef),
                        `p~surv~LR` = 1 - pchisq(diff$chisq, 1), 
                       `lo~surv~LR` = lo,
                       `hi~surv~LR` = hi ))
    } else {
        return(data.table(
                   `effect~surv~LR` = sign(coef),
                        `p~surv~LR` = 1 - pchisq(diff$chisq, 1)))
    }
}

#' @rdname dot-plot_survival
#' @export
fit_survival <- function(
    object, 
    assay      = assayNames(object)[1],
    percentile = 25, 
    samples    = if (ncol(object) < 50) TRUE else FALSE,
    verbose    = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
    assert_is_a_number(percentile)
    assert_all_are_in_left_open_range(percentile, 0, 50)
    event <- exprlevel <- timetoevent <- value <- NULL
# Fit
    if (verbose)  cmessage('\t\tsurvival ~ exprlevel')                                     # Filter across
    object %<>% filter_samples(!is.na(event) & !is.na(timetoevent))
    dt <- sumexp_to_longdt(object, assay = assay, svars = c('event', 'timetoevent'))       # Melt
    if (verbose)  message(
        sprintf("\t\t\texprlevel = 'Lo' (exprvalue <= %d%%)", percentile),                 # Dichotomize
        sprintf(            "  or  'Hi' (exprvalue >= %d%%)", 100 - percentile))
    dt %<>% dichotomize_exprs(percentile = percentile)                                     # Filter within 
    dt <- dt[, .SD[sum(event==1 & !is.na(value))>=3], by = c('feature_id', 'exprlevel')]   #    3 events     per feature/exprlevel
    dt <- dt[, .SD[    length(unique(exprlevel))==2], by = c('feature_id')             ]   #    2 exprlevels per feature
    if (verbose)  cmessage('\t\t\tp  =  survdiff(Surv(timetoevent, event) ~ exprlevel)')
    if (verbose)  cmessage('\t\t\teffect = coxph(Surv(timetoevent, event) ~ exprvalue)')
    dt %<>% extract(, .fit_survival(.SD, samples = samples), by = 'feature_id')            # Fit survival
# Return
    oldnames <- names(dt) %>% extract(stri_detect_regex(., '[~]LR$'))
    newnames <- paste0(oldnames, percentile)
    setnames(dt, oldnames, newnames) 
    for (col in newnames)  object[[col]] <- NULL
    object %<>% merge_fdt(dt)
    object
}



#' Plot survival 
#' @param object      SummarizedExperiment
#' @param assay       string
#' @param percentile  percentage (not greater than 50)
#' @param samples     TRUE or FALSE : record which samples in which stratum ?
#' @param verbose     TRUE or FALSE
#' @param title       string
#' @param subtitle    string
#' @param palette     color vector
#' @param n           number
#' @param ncol        number
#' @param nrow        number
#' @param file        filepath
#' @param width       number
#' @param height      number
#' @return ggsurvplot
#' @examples 
#' require(magrittr)
#' file <- download_tcga_example()
#' if (!is.null(file) & requireNamespace('survminer')){
#' # Read
#'     object <- readRDS(file)
#'     object %<>% extract(, .$sample_type == 'T')
#'     object %<>% extract(c('UGT3A2', 'NSUN3', 'XRCC4', 'WNT10A'), )
#' # Fit
#'     object %<>% fit_survival()
#'     object %<>% fit_survival(percentile = 50)
#'     fdt(object)
#' # Plot
#'     plot_survival(object)
#'     p1 <- .plot_survival(object[1, ])
#'     p2 <- .plot_survival(object[2, ])
#' }
#' @export
.plot_survival <- function(
    object,
    assay      = assayNames(object)[1],
    percentile = 25,
    title      = paste0(assay, ' ', percentile, '%'),
    subtitle   = NULL, #paste0(assay, ': ', percentile, '% split'),
    palette    = c("#009999", "#ff5050")
){
# Assert
    if (!requireNamespace('survminer', quietly = TRUE)){
        message("BiocManager::install('survminer'). Then re-run.")
        return(object) 
    }
    assert_is_valid_sumexp(object)
    if (nrow(object)==0)  return(empty_survplot())
    assert_is_subset(c('event', 'timetoevent'), svars(object))
    assert_is_identical_to_true(nrow(object)==1)
    feature <- unique(fdata(object)$feature_id)
    title %<>% paste(feature, ., sep = ' : ')
    assert_is_scalar(feature)
    assert_all_are_less_than_or_equal_to(percentile, 50)
    value <- exprlevel <- NULL
# Prepare
    subdt <- sumexp_to_longdt(
        object, 
        assay = assay, 
        svars = c('event', 'timetoevent'))
    subdt %<>% dichotomize_exprs(percentile = percentile)
# Plot
    fit <- survival::survfit(survival::Surv(timetoevent, event) ~ exprlevel, data = subdt)
    survminer::ggsurvplot(
        fit, data = subdt, conf.int = TRUE, palette = palette,
        risk.table = TRUE, risk.table.col = 'strata', risk.table.height = 0.25, 
        pval = TRUE, ggtheme = theme_bw(), title = title, subtitle = subtitle,
        legend.labs = unique(subdt$exprlevel), legend.title = assay)
}

#' survival percentiles
#' @param object SummarizedExperiment
#' @return numeric vector
#' @export
percentiles <- function(object){
    pvar(object, coefs = 'surv')  %>% 
    substr(nchar(.)-1, nchar(.))  %>% 
    as.numeric()
}


#' @rdname dot-plot_survival
#' @export
plot_survival <- function(
    object, 
    assay = assayNames(object)[1], 
    percentile = percentiles(object),
    title = paste0(assay, ' ', percentile, '%'),
    subtitle  = NULL,
    palette = c("#009999", "#ff5050"),
    n = 4,
    ncol = 4, 
    nrow = length(percentile), 
    file = NULL, 
    width  = 7*ncol, 
    height = 7*nrow
    
){
# Extract
    object %<>% order_on_p(fit = paste0('LR', percentile), coefs = 'surv')
    n %<>% min(nrow(object))
    object %<>% extract(1:n, )
# Plot
    if (!is.null(file))  pdf(file, width = width, height = height)
    npages <- ceiling(nrow(object)/ncol)
    for (i in 1:npages){
        idx1 <- (i-1)*ncol+1
        idxn <- max(i*ncol, nrow(object))
        idx <- idx1:idxn
        objlist <- object[idx, ]
        objlist %<>% split_features(by = 'feature_id')
        plots <- mapply(
            .plot_survival, 
            object     = rep(objlist, each = length(percentile)), 
            percentile = rep(percentile, times = length(objlist)),
            MoreArgs = list(assay = assay, palette = palette), SIMPLIFY = FALSE)
        survminer::arrange_ggsurvplots(plots, nrow = nrow, ncol = ncol)
    }
    if (!is.null(file))  dev.off()
}


