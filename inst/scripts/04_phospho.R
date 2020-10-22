phospho_dt

# 8967
phospho_dt[, length(unique(phospho_id))]

# 1281: NA in all samples
phospho_dt[ , .SD[all(is.na(protein_value))], by = c('phospho_id')][, length(unique(phospho_id))]
phospho_dt[ , .SD[all(is.na(protein_value))], by = c('phospho_id')][phospho_id==unique(phospho_id)[1]]
phospho_dt[ , .SD[all(is.na(protein_value))], by = c('phospho_id')][phospho_id==unique(phospho_id)[2]]
phospho_dt[ , .SD[all(is.na(protein_value))], by = c('phospho_id')][phospho_id==unique(phospho_id)[4]]


# 0: NA in all samples of some subgroup
phospho_dt[, all.na  := all(is.na(protein_value)) , by = c('phospho_id', 'subgroup')]
phospho_dt[, some.na := !all(is.na(protein_value)) & any(is.na(protein_value)) , by = c('phospho_id', 'subgroup')]

phospho_dt[, .SD[ all(all.na==TRUE)                     ], by = 'phospho_id'][phospho_id==unique(phospho_id)[1]] #  all subgroups have all missing
phospho_dt[, .SD[!all(all.na==TRUE) & all(some.na==TRUE)], by = 'phospho_id'][phospho_id==unique(phospho_id)[1]] #  all subgroups have some missing
phospho_dt[, .SD[!all(all.na==TRUE) & any(all.na ==TRUE)], by = 'phospho_id'][phospho_id==unique(phospho_id)[1]] # some subgroups have all missing: none


phospho_dt[, .SD[!all(all_na_subgroup==TRUE) & any(all_na==TRUE)], by = 'phospho_id']

# NA in some samples of all subgroups
phospho_dt

phospho_dt[, .SD[!all(na_subgroup==TRUE) & any(na_subgroup==TRUE)], by = 'phospho_id']

phospho_dt[ , .SD[, all(is.na(protein_value)), by = 'subgroup'], by = c('phospho_id')][, length(unique(phospho_id))]


# 925: some proteinvalues NA
phospho_dt[ , .SD[!all(is.na(protein_value)) & any(is.na(protein_value))], by = c('phospho_id')][, length(unique(phospho_id))]
phospho_dt[ , .SD[!all(is.na(protein_value)) & any(is.na(protein_value))], by = c('phospho_id')][phospho_id==unique(phospho_id)[1]]
