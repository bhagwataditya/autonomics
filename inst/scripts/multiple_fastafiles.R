# Go to R_user_dir('autonomics', 'cache')
# Create dir `uniprot_sprot-only2014_01`
# Download https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2014_01/knowledgebase/uniprot_sprot-only2014_01.tar.gz
# Untar etc
dir <- file.path(R_user_dir('autonomics', 'cache'))
dir %<>% file.path('uniprot/uniprot_sprot-only2014_01/')
file <- file.path(dir, 'uniprot_sprot-only2014_01.fasta')
file.exists(file)
fasta <- Biostrings::readAAStringSet(file)
organism <- names(fasta) %>% split_extract_fixed(' ', 1) %>% split_extract_fixed('_', 2)
human <- fasta %>% extract(organism == 'HUMAN')
mouse <- fasta %>% extract(organism == 'MOUSE')
Biostrings::writeXStringSet(human, file.path(dir, 'uniprot_sprot_human.fasta'))
Biostrings::writeXStringSet(mouse, file.path(dir, 'uniprot_sprot_mouse.fasta'))