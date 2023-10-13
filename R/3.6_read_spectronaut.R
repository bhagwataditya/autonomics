# file <- '../spectronaut/20231006_095632_hsfa-1-10_Report.tsv'
.read_spectronaut_sites <- function(
    file, modification = 'Phospho (STY)', pg.pep = 1, verbose = TRUE
){
# Assert
    assert_all_are_existing_files(file)
    assert_is_subset(c('R.FileName', 'PG.ProteinGroups', 'PTM.CollapseKey', 'PG.Quantity', 
                       'PTM.Quantity'), cols(file))
    assert_is_subset(modification, unique(fread(file, select = 'PTM.ModificationTitle')[[1]]))
    assert_is_fraction(pg.pep)
    assert_is_a_bool(verbose)
# Read
    scols <- c('R.FileName', 'R.Condition', 'R.Replicate')
    fcols <- c('PG.Genes', 'PG.ProteinGroups', 'PTM.CollapseKey', 'PTM.ModificationTitle')
    vcols <- c('PG.Quantity', 'PTM.Quantity', 'PG.PEP')
    allcols <- c(scols, fcols, vcols)
    allcols %<>% intersect(cols(file))
    dt <- fread(file, select = c(scols, fcols, vcols))

# Filter
    # on modification
        n0 <- dt[, length(unique(PTM.CollapseKey))]
        dt %<>% extract(PTM.ModificationTitle == modification)
        n1 <- dt[, length(unique(PTM.CollapseKey))]
        if (verbose & n1<n0)  message('\t\tRetain ', n1, '/', n0, ' sites with a ', modification, ' modification')
    # on pep
        dt %<>% extract(PG.PEP < pg.pep)
        n2 <- dt[, length(unique(PTM.CollapseKey))]
        if (verbose & n2<n1)  message('\t\tRetain ', n2, '/', n1, ' sites: pg.pep < ', pg.pep)
# Rename
    setnames(dt, 'R.FileName',            'sample_id')
    setnames(dt, 'R.Condition',           'subgroup')
    setnames(dt, 'R.Replicate',           'replicate')
    setnames(dt, 'PG.Genes',              'gene')
    setnames(dt, 'PG.ProteinGroups',      'uniprot')
    setnames(dt, 'PTM.CollapseKey',       'site')
    setnames(dt, 'PTM.ModificationTitle', 'modification')
    setnames(dt, 'PG.Quantity',           'pg.intensity')
    setnames(dt, 'PTM.Quantity',          'site.intensity')
    setnames(dt, 'PG.PEP',                'pep')

# Return
    dt[]
}

read_spectronaut_sites <- function(
    file, modification = 'Phospho (STY)', pg.pep = 1, verbose = TRUE
){
# read
    dt <- .read_spectronaut_sites(
            file = file, modification = modification, pg.pep = pg.pep, verbose = verbose)
    log2sites    <- dcast(dt, site ~ sample_id, value.var = 'site.intensity')
    log2proteins <- dcast(dt, site ~ sample_id, value.var = 'pg.intensity')
    log2sites    %<>% dt2mat() %>% log2()
    log2proteins %<>% dt2mat() %>% log2()
# assays
    object <- list(log2sites    = log2sites, 
                   log2proteins = log2proteins, 
                   log2diffs    = log2sites - log2proteins)
    object %<>% SummarizedExperiment::SummarizedExperiment()
    sdt(object)$sample_id  <- snames(object)
    fdt(object)$feature_id <- fnames(object)
# sdt/fdt
    scols <- c('sample_id', 'subgroup', 'replicate')
    fcols <- c('gene', 'uniprot', 'site', 'modification', 'pep')
    scols %<>% intersect(names(dt))
    fcols %<>% intersect(names(dt))
    sdt0 <- unique(dt[, scols, with = FALSE])
    fdt0 <- unique(dt[, fcols, with = FALSE])
    setnames(fdt0, 'site', 'feature_id')
    object %<>% merge_sdt(sdt0)
    object %<>% merge_fdt(fdt0)
    object 
}


