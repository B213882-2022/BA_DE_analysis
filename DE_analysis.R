library(affy)
library(limma)
library(EnhancedVolcano)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggfortify)

# prepare
setwd('/home/riley/online_test_2')
save.image('result.RData')
load('result.RData')

# get normalised data (only use single-end sequencing results here)
# normalised method: rpkm
# total: 121 Biliary Atresia samples (BA); 7 normal control (NC)
filenames <- list.files('/home/riley/online_test_2/GSE122340_RAW')
filename <- paste0('GSE122340_RAW/', filenames[1])
df <- read.csv(filename,sep='\t')
colnames(df) <- c('Symbol',filenames[1])
for(name in filenames[2:length(filenames)]){
  filename <- paste0('GSE122340_RAW/', name)
  df_new <- read.csv(filename, sep='\t')
  df <- cbind(df,df_new$Value)
}
filenames <- gsub('_.*gz','',filenames)
colnames(df) <- c('Symbol',filenames)
normalised_df <- df[,2:ncol(df)]
rownames(normalised_df) <- make.unique(df$Symbol)  # prevent duplicated gene names
rm(df_new)
rm(filename)
rm(name)
rm(df)

# annotation file
coldata <- read.csv('coldata.csv', row.names = 1)
match(colnames(normalised_df), rownames(coldata))  #check order of sample names

# gene filtering (same filtering standard as the orginal paper's)
rowSums(normalised_df<0.5)>ncol(normalised_df)/2
filted_df <- normalised_df[!rowSums(normalised_df<0.5)>ncol(normalised_df)/2,]

# check normalisation
png('boxplot.png',1000,700)
boxplot(log2(filted_df+0.1))
dev.off()

# rpkm is not suitable for DE analysis, but I got no choice. 
# This log transformation was based on the advice on https://support.bioconductor.org/p/56275/ 
# suggested by Gordon Smyth (one of the authors who developed limma and edgeR packages)
filted_df <- log2(filted_df+0.1)  

# sample similarities
# hierarchical clustering
values <- filted_df
colnames(values) <- make.unique(coldata$condition)
values[1:5,1:5]
hc <- hclust(as.dist(1-cor(values, method="pearson")), method="average")
png('hclust.png',1700,1000)
plot(hc)
dev.off()
# PCA
pca <- prcomp(t(values), scale=T)
pca_sum <- summary(pca)
# 2D plot
cluster <- coldata$condition
autoplot(pca,data=data.frame(t(values),cluster),colour='cluster', 
         main='2D PCA') + 
  geom_text(label=colnames(values),nudge_x=0.005, nudge_y=0.005, size=2.5)
ggsave('PCA_2D.jpg',device='jpg', width = 8, height = 6)

# DE analysis with limma
# design matrix
design <- model.matrix(~0+factor(c(rep(1,121),rep(2,7))))
colnames(design) <- c("BA","NC")
design

# contrast matrix
contrastmatrix <- makeContrasts(BA-NC,levels=design)
contrastmatrix

# calcualte DE genes with eBayes function from Limma package
fit <- lmFit(filted_df, design)
fit2 <- contrasts.fit(fit, contrastmatrix)
fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
DE <- topTable(fit2,coef=1, number=Inf, p.value = 0.05, adjust.method = "BH")
nrow(DE)
head(DE)
DE_all <- topTable(fit2,coef=1, number=Inf, adjust.method = "BH")
nrow(DE_all)

# volvano plot to show some of the DE genes
volcano_plot <- EnhancedVolcano(DE_all, lab = rownames(DE_all), x = 'logFC', y = 'adj.P.Val',
                                pCutoff = 0.01, FCcutoff = log2(1.5),drawConnectors = TRUE,
                                widthConnectors = 0.75,title = "Volcano plot (thres: log2(1.5), p.adjust=0.01)",
                                subtitle = "by function 'eBayes'",
                                max.overlaps = 2)
ggsave('volcano.pdf',volcano_plot, device='pdf', width = 15, height = 15)

# functions to enable 'geom_split_violin'
# ref: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


# create violin plot to show DE analysis results
selected_genes <- rownames(DE[DE$logFC>0,][1:5,])  # positive fold-change
selected_genes
test_df <- t(filted_df[selected_genes,])
test_df <- data.frame(test_df,
                      condition=coldata,
                      sample=rownames(test_df))
head(test_df)
test_df <- melt(test_df, id=c('condition', 'sample'))
colnames(test_df) <- c('condition','sample','gene','expression')
head(test_df)
violin_plot <- ggplot(data=test_df, aes(x=gene, y=expression, fill=condition)) + 
                  geom_split_violin() +
                  ggtitle('Top 5 DE genes') +
                  theme(plot.title = element_text(hjust = 0.5))
ggsave('violin.pdf',violin_plot, device='pdf', width = 10, height = 7)

nega_selected_genes <- rownames(DE[order(DE$logFC,decreasing=FALSE),][1:5,])
nega_selected_genes
nega_test_df <- t(filted_df[nega_selected_genes,])
nega_test_df <- data.frame(nega_test_df,
                      condition=coldata,
                      sample=rownames(nega_test_df))
head(nega_test_df)
nega_test_df <- melt(nega_test_df, id=c('condition', 'sample'))
colnames(nega_test_df) <- c('condition','sample','gene','expression')
head(nega_test_df)
nega_violin_plot <- ggplot(data=nega_test_df, aes(x=gene, y=expression, fill=condition)) + 
                      geom_split_violin() +
                      ggtitle('Top 5 DE genes (negative)') +
                      theme(plot.title = element_text(hjust = 0.5))
ggsave('nega_violin.pdf',nega_violin_plot, device='pdf', width = 10, height = 7)



