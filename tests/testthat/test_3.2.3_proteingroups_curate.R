context('curate_uniprots')

profile   <- download_data('billing19.proteingroups.txt')
fosfile   <- download_data('billing19.phosphosites.txt')
fastafile <- download_data('uniprot_hsa_20140515.fasta')
pro0 <- read_pro(profile, curate = FALSE)
fos0 <- read_fos(fosfile, curate = FALSE)
pro1 <- curate_uniprots(pro0)
fos1 <- curate_uniprots(fos0)
pro2 <- curate_uniprots(pro0, fastafile = fastafile)
fos2 <- curate_uniprots(fos0, fastafile = fastafile)

test_equivalence <- function(dt0, dt1){
    # nrows
        dt0 <- copy(dt0)
    # proId/fosIds
        idcol <- if ('fosId' %in% names(dt0))  'fosId' else 'proId' 
        expect_identical(dt0[[idcol]], dt1[[idcol]])
    # snames identical
        quantity <- guess_maxquant_quantity(names(dt0))
        pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
        expect_identical(grep(pattern, names(dt0), value = TRUE),
                         grep(pattern, names(dt1), value = TRUE))
    # values identical
        expect_identical(dt0[, .SD, .SDcols = patterns(pattern)],
                         dt1[, .SD, .SDcols = patterns(pattern)])
    # Uniprots identical
        expect_identical(dt0[, Uniprot], dt1[, Uniprot])
    #`Potential contaminant` identical
        expect_identical(dt0[, `Potential contaminant`], 
                         dt1[, `Potential contaminant`])
    # Reverse identical
        expect_identical(dt0[, Reverse], dt1[, Reverse])
    # Curated is subset of Uniprot
        dt01 <- merge(dt0[, c(idcol, 'Uniprot', 'Reverse', 'Potential contaminant'), with = FALSE], 
                      dt1[, c(idcol, 'Curated'), with = FALSE], by = idcol, sort = FALSE)
        dt01 %<>% uncollapse(Curated, sep = ';')
        dt01 %<>% uncollapse(Uniprot, sep = ';')
        dt01 %<>% extract(Reverse=='' & `Potential contaminant`=='')
        expect_true(all(dt01[, is_subset(Curated, Uniprot), by = idcol][[2]]))
}

test_equivalence(pro0, pro1)
test_equivalence(pro0, pro2)
test_equivalence(fos0, fos1)
test_equivalence(fos0, fos2)

