# Go to R_user_dir('autonomics', 'cache')
# Create dir `uniprot_sprot-only2014_01`
# Download https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2014_01/knowledgebase/uniprot_sprot-only2014_01.tar.gz
# Untar etc
dir <- file.path(R_user_dir('autonomics', 'cache'))
dir %<>% file.path('uniprot/uniprot_sprot-only2014_01/')
file <- file.path(dir, 'uniprot_sprot-only2014_01.fasta')
file.exists(file)
fasta <- seqinr::read.fasta(file)
human <- fasta %>% extract(split_extract_fixed(names(.), '_', 2) == 'HUMAN')
mouse <- fasta %>% extract(split_extract_fixed(names(.), '_', 2) == 'MOUSE')
seqinr::write.fasta(human, names = names(human), file.out = file.path(dir, 'uniprot_sprot_human.fasta'))
seqinr::write.fasta(mouse, names = names(mouse), file.out = file.path(dir, 'uniprot_sprot_mouse.fasta'))