# https://f1000research.com/articles/5-1438
require(Rsubread)
require(magrittr)

# Read quality
QS <- qualityScores("~/importomicscache/stemcells.fastq/SRR1909613_1.fastq")
boxplot(QS, ylab="Quality score", xlab="Base position", main="SRR1909613_1.fastq", cex=0.25, col="orange")

# Index
url <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz'
download.file(url, '~/importomicscache/genomes/GRCh38.primary_assembly.genome.fa.gz')
setwd('~/importomicscache/genomes/GRCh38')
buildindex(
    basename = "GRCh38",
    reference = '~/importomicscache/genomes/GRCh38.primary_assembly.genome.fa.gz')

# Align
dir('~/importomicscache/stemcells.fastq')
reads1 <- list.files('~/importomicscache/datasets/stemcells/stemcells.fastq', pattern = '_1', full.names = TRUE)
reads2 <- list.files('~/importomicscache/datasets/stemcells/stemcells.fastq', pattern = '_2', full.names = TRUE)
sample_design <- c(
    SRR1909613 = 'BM.R1', SRR1909637 = 'BM.R2', SRR1909638 = 'BM.R3', SRR1909639 = 'BM.R4',
    SRR1926136 = 'EM.R1', SRR1926132 = 'EM.R2', SRR1926133 = 'EM.R3',
    SRR1926134 = 'E.R1',  SRR1926135 = 'E.R2',  SRR1867792 = 'E.R3')
bams <- paste0(unname(sample_design[sub('_1.fastq', '', basename(reads1))]), '.bam')
bams %<>% paste0('~/importomicscache/datasets/stemcells/stemcells.bam/', .)
setwd('~/importomicscache/genomes/GRCh38')
#align(index = "GRCh38", readfile1 = reads1, readfile2 = reads2, output_file = bams, nthreads = 10)
    #    BM.R1.bam    BM.R2    BM.R3    BM.R4    EM.R1    EM.R2    EM.R3   E.R1    E.R2    E.R3
    #                  10       18        13         18     17        9     3      9              => approx 2 hours in total

# Read bam into analysis-ready SummarizedExperiment
object <- read_bam('~/importomicscache/datasets/stemcells/stemcells.bam/', ispaired = TRUE)

# Alternative: specify manually
# Count
#fc <- featureCounts(bams, annot.inbuilt = 'hg38') # 5 mins

object <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = fc$counts))
fdata(object) <- data.frame(feature_id = fnames(object), row.names = fnames(object), stringsAsFactors = FALSE)
fdata(object)$feature_name <- AnnotationDbi::mapIds(
                            org.Hs.eg.db::org.Hs.eg.db, fdata(object)$feature_id, keytype = 'ENTREZID', column = 'SYMBOL')
snames(object) %<>% substr(1, nchar(.)-4)
sdata(object) <- data.frame(sample_id = snames(object), row.names = snames(object), stringsAsFactors = FALSE)
object$subgroup <-  guess_subgroup_values(object$sample_id)
object

# Filter
object %<>% subset(!is.na(feature_name), ) # 28 395 -> 27 885

