# What does PLS do to truly uncorrelated samples
#-----------------------------------------------
## All correlated
devtools::load_all()
set.seed(42)
controlmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10,  2), geneC = rnorm(10,3), geneD = rnorm(10,4))
patientmat <- rbind(geneA = rnorm(10, 4), geneB = rnorm(10, -1), geneC = rnorm(10,6), geneD = rnorm(10,1))
colnames(controlmat) <- sprintf('C%02d', 1:10)
colnames(patientmat) <- sprintf('P%02d', 1:10)
omat <- cbind(controlmat, patientmat)
obj <- SummarizedExperiment(list(exprs = omat))
obj$sample_id <- snames(obj)
fdt(obj)$feature_id <- fnames(obj)
obj$subgroup <- obj$sample_id %>% substr(1,1)

S4Vectors::metadata(pca(obj, ndim = 4))
S4Vectors::metadata(pls(obj, ndim = 4))

opca <- biplot(pca(obj, ndim = 4), x = 'pca1', y = 'pca2')
opls <- biplot(pls(obj, ndim = 4), x = 'pls1', y = 'pls2')

## None correlated
set.seed(42)
controlmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10, 2), geneC = rnorm(10, 3), geneD = rnorm(10, 4))
patientmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10, 2), geneC = rnorm(10, 3), geneD = rnorm(10, 4))
colnames(controlmat) <- sprintf('C%02d', 1:10)
colnames(patientmat) <- sprintf('P%02d', 1:10)
umat <- cbind(controlmat, patientmat)
ubj <- SummarizedExperiment(list(exprs = umat))
ubj$sample_id <- snames(ubj)
fdt(ubj)$feature_id <- fnames(ubj)
ubj$subgroup <- ubj$sample_id %>% substr(1,1)

ubj$sample_id <- snames(ubj)
fdt(ubj)$feature_id <- fnames(ubj)
ubj$subgroup <- ubj$sample_id %>% substr(1,1)

S4Vectors::metadata(pca(ubj, ndim = 4))
S4Vectors::metadata(pls(ubj, ndim = 4))

upca <- biplot(pca(ubj, ndim = 4), x = 'pca1', y = 'pca2')
upls <- biplot(pls(ubj, ndim = 4), x = 'pls1', y = 'pls2')

# Plot
gridExtra::grid.arrange(
    opca + ggtitle( 'Corr PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    opls + ggtitle(  'Cor PLS') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    upca + ggtitle('UnCor PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    upls + ggtitle('Uncor PLs') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    nrow = 2
)

# Effect of centering operations on PCA
#--------------------------------------
set.seed(42)
controlmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10, 2), geneC = rnorm(10, 3), geneD = rnorm(10,4))
patientmat <- rbind(geneA = rnorm(10, 4), geneB = rnorm(10, 5), geneC = rnorm(10, 6), geneD = rnorm(10,7))
colnames(controlmat) <- sprintf('C%02d', 1:10)
colnames(patientmat) <- sprintf('P%02d', 1:10)
omat <- cbind(controlmat, patientmat)
obj <- SummarizedExperiment(list(exprs = omat))
obj$sample_id <- snames(obj)
fdt(obj)$feature_id <- fnames(obj)
obj$subgroup <- obj$sample_id %>% substr(1,1)

## None correlated
set.seed(42)
controlmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10, 2), geneC = rnorm(10, 3), geneD = rnorm(10, 4))
patientmat <- rbind(geneA = rnorm(10, 1), geneB = rnorm(10, 2), geneC = rnorm(10, 3), geneD = rnorm(10, 4))
colnames(controlmat) <- sprintf('C%02d', 1:10)
colnames(patientmat) <- sprintf('P%02d', 1:10)
umat <- cbind(controlmat, patientmat)
ubj <- SummarizedExperiment(list(exprs = umat))
ubj$sample_id <- snames(ubj)
fdt(ubj)$feature_id <- fnames(ubj)
ubj$subgroup <- ubj$sample_id %>% substr(1,1)

ubj$sample_id <- snames(ubj)
fdt(ubj)$feature_id <- fnames(ubj)
ubj$subgroup <- ubj$sample_id %>% substr(1,1)

# Plot
gridExtra::grid.arrange(
    biplot(pca(obj)) + ggtitle( 'Corr PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    biplot(pca(ubj)) + ggtitle('UnCor PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    nrow = 1
)

gridExtra::grid.arrange(
    biplot(pca(obj, center_samples = FALSE)) + ggtitle( 'Corr PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    biplot(pca(ubj, center_samples = FALSE)) + ggtitle('UnCor PCA') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    nrow = 1
)

# Atkin data
file <- download_data('atkin18.metabolon.xlsx')
object <- read_metabolon(file)
gridExtra::grid.arrange(
    biplot(pca(object                        )) + ggtitle( 'center_samples = TRUE')  + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'), 
    biplot(pca(object, center_samples = FALSE)) + ggtitle( 'center_samples = FALSE') + theme(plot.title = element_text(hjust = 0.5)) + guides(color = 'none'),
    nrow = 1
)


