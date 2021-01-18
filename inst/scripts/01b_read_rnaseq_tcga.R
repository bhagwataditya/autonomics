# Download tcga rnacounts
# https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html
require(TCGAbiolinks)
GDCprojects <- getGDCprojects()
data.table(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-LIHC") # liver hepatocelluar carcinoma
query_TCGA <- GDCquery(project  = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
          experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Primary solid Tumor", "Solid Tissue Normal"))
GDCdownload(    query = query_TCGA, directory='~/autonomicscache/datasets/tcga')
tcgaobject <- GDCprepare(query_TCGA, directory='~/autonomicscache/datasets/tcga')
idx <- vapply(colData(tcgaobject), is.list, logical(1))
colData(tcgaobject) %<>% extract(, !idx)
matchedpatients <- as.data.table(colData(tcgaobject))[, .N, by = 'patient'][N==2]$patient
tcgaobject %<>% filter_samples(patient %in% matchedpatients)

counts1  <- data.table(assays(tcgaobject)$`HTSeq - Counts`, keep.rownames = TRUE)
coldata1 <- as.data.table(colData(tcgaobject))
setnames(counts1, 'rn', 'gene')
fwrite(counts1,  "~/autonomicscache/datasets/tcga/TCGA-LIHC/lihc.counts.tsv",   sep='\t')
fwrite(coldata1, "~/autonomicscache/datasets/tcga/TCGA-LIHC/lihc.coldata1.tsv", sep='\t')

# Test read_rnaseq_counts
file         <- "~/autonomicscache/datasets/tcga/TCGA-LIHC/lihc.counts.tsv"
samplefile <- "~/autonomicscache/datasets/tcga/TCGA-LIHC/lihc.coldata1.tsv"
.read_rnaseq_counts(file = file,  samplefile = samplefile,
                    sampleidvar = 'barcode', subgroupvar = 'definition')

read_rnaseq_counts(file = file, samplefile = samplefile,
                    sampleidvar = 'barcode', subgroupvar = 'definition')
