#-------------------------------------------------------------------------------
#
#               hdl proteome watch
#
#-------------------------------------------------------------------------------

#' hdl proteomewatch proteins
#' @return string vector: HDLProteomeWatch protein entries
#' @examples
#' hdlproteins()
#' @export
hdlproteins <- function(){
    url <- 'https://homepages.uc.edu/~davidswm/HDL%20Proteome%20Watch%202021.xlsx'
    dir <- file.path(tools::R_user_dir('autonomics', 'cache'), 'hdlproteomewatch')
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file <- file.path(dir, basename(url))
    if (!file.exists(file))  download.file(url, destfile = file, mode = 'wb')
    hdlproteins <- readxl::read_excel(file, sheet = 3, range = 'C9:E944', col_names = c('entry', 'mw', 'uniprot'))
    hdlproteins %<>% extract2('entry')
    hdlproteins
}

#' Tag hdlproteins
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @examples 
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file)
#' object %<>% tag_hdlproteins()
#' fdt(object)
#' @export
tag_hdlproteins <- function(object, verbose = TRUE){
    hdlproteins <- autonomics::hdlproteins()
    fdt0 <- fdt(object)[, c('feature_id', 'protein')]
    fdt0 %<>% uncollapse(protein, sep = ';')
    fdt0 <- fdt0[, hdl := protein %in% hdlproteins]
    fdt0 <- fdt0[, .(hdl = any(hdl)), by = feature_id]
    fdt0$hdl %<>% as.numeric()
    object %<>% merge_fdt(fdt0)
    if (verbose)  message(sprintf(
        '%d/%d HDL groups, %d/%d HDLProteomeWatch proteins', 
        sum(fdt(object)$hdl), 
        nrow(object), 
        length(intersect(fdt(object)[hdl==1, protein], hdlproteins)),
        length(hdlproteins)))
    object
}


#-------------------------------------------------------------------------------
#
#               open targets
#
#-------------------------------------------------------------------------------

#' opentargets dir
#' @export
OPENTARGETSDIR <- file.path(tools::R_user_dir('autonomics', 'cache'), 'opentargets', '22.04')

download_opentargets_targets <- function(){
    if (!dir.exists(file.path(OPENTARGETSDIR, 'targets'))){
        ftpdir <- 'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.04/output/etl/json/targets/'
        ftpfiles <- XML::getHTMLLinks(ftpdir) %>% setdiff(c('../', '_SUCCESS'))
        ftpfile <- ftpfiles[1]
        dir.create(file.path(OPENTARGETSDIR, 'targets'), showWarnings = FALSE, recursive = TRUE)
        for (i in seq_along(ftpfiles)){
            download.file(url = paste0(ftpdir, ftpfiles[i]),
                          destfile = file.path(OPENTARGETSDIR, 'targets', ftpfiles[i]))
        }
    }
    return(file.path(OPENTARGETSDIR, 'targets'))
}

null_to_str <- function(x){  x[is.null(x)] <- ''; x }

extract_proteinids <- function(x){
    x %<>% extract2('proteinIds')
    if (is.null(x)) return('')
    x %<>% data.table()
    x %<>% extract(source %in% c('uniprot_swissprot', 'uniprot_trembl'))
    x %<>% extract2('id')
    x %<>% paste0(collapse = ';')
    x
}

extract_functiondescriptions <- function(x){
    x %<>% extract2('functionDescriptions')
    x %<>% paste0(collapse = ';')
    x %<>% null_to_str()
    x
}

# file <- list.files(file.path(OPENTARGETSDIR, 'targets'))[1]
# read_opentargets_targets(file) 
.read_opentargets_targets <- function(file){
    lines <- readLines(file.path(OPENTARGETSDIR, 'targets', file))
    lines %<>% lapply(jsonlite::fromJSON)
    lines %<>% lapply( function(x){ data.table(
                                        ensembl    = x$id,
                                        genesymbol = x$approvedSymbol,
                                        genename   = x$approvedName,
                                        uniprot    = extract_proteinids(x),
                                       `function`  = extract_functiondescriptions(x))})
   lines %<>% rbindlist()
   lines
}


save_opentargets_targets <- function(){
    files <- dir(file.path(OPENTARGETSDIR, 'targets'))
    dt <- lapply(files, .read_opentargets_targets)
    dt %<>% rbindlist()
    fwrite(dt, file.path(OPENTARGETSDIR, 'targets.tsv'), sep = '\t')
}

#' Add opentargets annotations
#' @param object SummarizedExperiment
#' @return  SummarizedExperiment
#' @examples 
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file)
#' object %<>% add_opentargets_by_uniprot()
#' @export
add_opentargets_by_uniprot <- function(object){
    
    file <- file.path(OPENTARGETSDIR, 'targets.tsv')
    targetsdt <- fread(file, select = c('uniprot', 'genesymbol', 'genename', 'function'))
    targetsdt %<>% unique()
    targetsdt %<>% separate_rows(uniprot)
    targetsdt %<>% data.table()
    targetsdt %<>% extract(uniprot != '')
    targetsdt %<>% unique()
    targetsdt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = 'uniprot')

    idvar <- if ('fosId' %in% fvars(object)) 'fosId' else 'proId'
    fdt0 <- fdt(object)[, c(idvar, 'canonical'), with = FALSE]
    fdt0 %<>% separate_rows(canonical)
    fdt0 %<>% data.table()
    fdt0 %<>% merge(targetsdt, by.x = 'canonical', by.y = 'uniprot', all.x = TRUE, sort = FALSE)
    fdt0[, canonical := NULL]
    fdt0 %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idvar)
    object %<>% merge_fdata(fdt0, by.x = idvar, by.y = idvar)
    object    
}


