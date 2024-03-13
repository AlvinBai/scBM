### loading packages -------------------------------------------------
rm(list = ls())
options(stringasfac)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(homologene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(ggrepel)

### functions --------------------------------------
DoubletEst <- function (object, pN = 0.25, ndims = 30) {
  for (pkg in c("Seurat", "DoubletFinder")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste(pkg, " package needed for this function to work. Please install it.", 
                 sep = ""), call. = FALSE)
    }
  }
  library(Seurat)
  library(DoubletFinder)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object, features = rownames(object))
  object <- RunPCA(object)
  object <- FindNeighbors(object, dims = 1:ndims)
  object <- FindClusters(object)
  bcmvn <- find.pK(summarizeSweep(paramSweep_v3(object, PCs = 1:30), 
                                  GT = FALSE))
  homotypic.prop <- modelHomotypic(object@meta.data$seurat_clusters)
  nExp_poi <- round((0.076 * ncol(object)/10000) * length(colnames(object)))
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  object <- doubletFinder_v3(object, PCs = 1:10, pN = pN, 
                             pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
  doublet <- colnames(object)[object[[paste0("DF.classifications_", 
                                             as.character(pN), "_", as.character(mpK), "_", as.character(nExp_poi))]] == "Doublet"]
  return(doublet)
}
NoiseGeneDetec <- function(obs){
  cycle_genes <- intersect(c("UBE2C","HMGB2", "HMGN2", "TUBA1B", "MKI67",
                             "CCNB1", "TUBB", "TOP2A", "TUBB4B") %>% human2mouse() %>% pull(mouseGene), rownames(bm.int))
  cycle_score <-  as.numeric(colMeans(as.matrix(bm.int@assays$RNA@data)[cycle_genes,]))
  
  library(WGCNA)
  gene_data <- data.frame(cycling = cycle_score, 
                          t(as.matrix(bm.int@assays$RNA@data)),check.names = F)
  cor_data <- cor(gene_data, method = 'pearson')[ ,c('cycling'), drop = F]
  genes <- rownames(cor_data)[cor_data[,'cycling'] > 0.3]
  
  return(genes)
} ###  find genes related to genes offered in dataset
compareGO <- function(markers, orgdb = org.Mm.eg.db, ont = 'BP', pvaluecutoff = 0.05, width.y = 50, labs.y = 'Description', 
                      labs.x = 'Cluster', show.cat = 10) {
  library(clusterProfiler)
  library(stringr)
  library(ggplot2)
  
  gene.entrez <- bitr(markers$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = orgdb)
  markers.merge <- merge(markers, gene.entrez, by.x = 'gene', by.y = 'SYMBOL')
  sample <- split(markers.merge$ENTREZID, markers.merge$cluster)
  cc <- compareCluster(sample, fun = 'enrichGO', ont = ont, pvalueCutoff = pvaluecutoff, OrgDb = orgdb)
  
  dotplot(cc, showCategory = show.cat) + scale_y_discrete(label = function(x) str_wrap(x, width = width.y)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(y = labs.y, x = labs.x)
}

### reading data --------------------------------------
mat.ctrl <- Read10X('~/filtered_feature_bc_matrix_Ctrl/')
mat.r <- Read10X('~/filtered_feature_bc_matrix_R/')
mat.cb <- Read10X('~/filtered_feature_bc_matrix_CB/')

bm.ctrl <- CreateSeuratObject(mat.ctrl, min.cells = 3, min.features = 200)
bm.r <- CreateSeuratObject(mat.r, min.cells = 3, min.features = 200)
bm.cb <- CreateSeuratObject(mat.cb, min.cells = 3, min.features = 200)
bm.ctrl$group <- 'Ctrl'
bm.r$group <- 'R'
bm.cb$group <- 'CBD'
dim(bm.ctrl); dim(bm.r); dim(bm.cb)

### doublet removing and data integration ------------------------------------
doublets.ctrl <- DoubletEst(bm.ctrl)
doublets.r <- DoubletEst(bm.r)
doublets.cb <- DoubletEst(bm.cb)
save(doublets.ctrl, doublets.r, doublets.cb, file = '~/Doublets.RData')

genes.inter <- intersect(intersect(rownames(bm.ctrl), rownames(bm.r)), rownames(bm.cb))
object.list <- c(bm.ctrl[genes.inter, !colnames(bm.ctrl) %in% doublets.ctrl],
                 bm.r[genes.inter, !colnames(bm.r) %in% doublets.r], 
                 bm.cb[genes.inter, !colnames(bm.cb) %in% doublets.cb])
for (i in 1:3) {
  object.list[[i]] <- NormalizeData(object.list[[i]])
  object.list[[i]] <- FindVariableFeatures(object.list[[i]])
}
anchors <- FindIntegrationAnchors(object.list, dims = 1:20)
bm.int <- IntegrateData(anchors, dims = 1:20)

DefaultAssay(bm.int) <- 'RNA'
bm.int[['percent_mito']] <- PercentageFeatureSet(bm.int, pattern = '^mt-')
VlnPlot(bm.int, features = 'percent_mito') + geom_hline(yintercept = 25)
bm.int <- bm.int[, bm.int@meta.data$percent_mito <= 5]

### quality control------------------------------------------------
qc.mat <- bm.int@meta.data[, c('Group', 'nFeature_RNA', 'nCount_RNA')]
ggplot(qc.mat, aes(x = Group, y = nFeature_RNA)) + geom_violin(aes(fill = Group)) + geom_boxplot(aes(color = Group), fill = 'white', width = 0.2) +
  scale_fill_manual(values = c('Ctrl' = 'forestgreen', 'R' = 'blue', 'CBD' = 'firebrick2')) +
  scale_color_manual(values = c('Ctrl' = 'forestgreen', 'R' = 'blue', 'CBD' = 'firebrick2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA)) +
ggplot(qc.mat, aes(x = Group, y = nCount_RNA)) + geom_violin(aes(fill = Group)) + geom_boxplot(aes(color = Group), fill = 'white', width = 0.2) +
  scale_fill_manual(values = c('Ctrl' = 'forestgreen', 'R' = 'blue', 'CBD' = 'firebrick2')) +
  scale_color_manual(values = c('Ctrl' = 'forestgreen', 'R' = 'blue', 'CBD' = 'firebrick2')) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA))

### clustering and markers -------------------------------------------------------
DefaultAssay(bm.int) <- 'integrated'
noise.genes <- NoiseGeneDetec(bm.int)
save(noise.genes, file = '~/noiesGenes.RData')

bm.int <- ScaleData(bm.int, features = setdiff(rownames(bm.int), noise.genes), vars.to.regress = 'percent_mito')
bm.int <- RunPCA(bm.int, features = setdiff(VariableFeatures(bm.int), noise.genes))
ElbowPlot(bm.int, ndims = 50)
bm.int <- FindNeighbors(bm.int, dims = 1:30)
bm.int <- FindClusters(bm.int, resolution = 0.35)
DimPlot(bm.int, label = TRUE)

bm.int <- RunUMAP(bm.int, dims = 1:20, n.neighbors = 20, seed.use = 1000, n.epochs = 800, min.dist = 0.5, local.connectivity = 2)
DimPlot(bm.int, label = TRUE, cols = color.int) + DimPlot(bm.int, group.by = 'Phase')

DefaultAssay(bm.int) <- 'RNA'
bm.int <- ScaleData(bm.int)
bm.int <- CellCycleScoring(bm.int, s.features = cc.genes.updated.2019$s.genes %>% human2mouse() %>% pull(mouseGene),
                           g2m.features = cc.genes.updated.2019$g2m.genes %>% human2mouse() %>% pull(mouseGene))


new.ident <- c('preB', 'CMP1', 'HSC', 'MPP', 'CMP2', 'Ery1', 'proB', 'GraP', 'immatureB', 'Gra', 'CLP', 'MEP', 'MC', 'pDC', 'Mono/Mac', 'NK', 'Ery2')
names(new.ident) <- levels(bm.int)
bm.int <- RenameIdents(bm.int, new.ident)
Idents(bm.int) <- factor(Idents(bm.int), levels = c('HSC', 'MPP', 'CMP1', 'CMP2', 'Mono/Mac', 'GraP', 'Gra', 'MC', 'MEP', 'Ery1', 'Ery2', 'CLP', 'proB', 'preB', 'immatureB',
                                                    'pDC', 'NK'))
bm.int$cluster <- Idents(bm.int)
bm.int$group <- factor(bm.int$group, levels = c('Ctrl', 'R', 'CBD'))

### Figure 4B Dimplot for all clusters ------------------------------------------------------------------------------
color.int <- c("maroon4","#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#990000","#9900cc","#66FF66","slateblue1","mediumblue","#CC0033","#FF0000")
names(color.int) <- levels(bm.int)
CellDimPlot(bm.int, group.by = 'cluster', palcolor = color.int, theme_use = 'theme_blank', show_stat = F)

### Figure S4B Splited dimplot for groups  ---------------------------------------------------------------------------
DimPlot(bm.int, split.by = 'group', cols = color.int)

### Figure 4C Heatmaps for expression of lineage markers -------------------------------------------------------------
Idents(bm.int) <- 'cluster'
feature.genes <- c('Hlf', 'Procr', 'Flt3', 'Mpo', 'Ms4a3', 'Csf1r', 'Cd68', 'S100a8', 'S100a9', 'Cpa3', 'Fcer1a', 'Vamp5', 'Pf4', 'Car1', 'Gypa', 'Hbb-a2', 'Lef1',
                   'H2ax', 'Cecr2', 'Ms4a1', 'Cd40', 'Siglech', 'Irf8', 'Gzma', 'Ms4a4b')
mat.expr <- AverageExpression(bm.int, features = feature.genes, assays = 'RNA')
mat.expr <- t(scale(t(mat.expr$RNA)))
meta.expr <- data.frame(row.names = levels(bm.int), Cluster = levels(bm.int))
pheatmap(mat.expr, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = meta.expr, annotation_colors = list(Cluster = color.int.new), show_colnames = FALSE, 
         color = colorRampPalette(c('slateblue4', 'grey', 'white', 'orange', 'red', 'firebrick4'))(40))

### Figure S4C Dotplot for top5 DEGs ------------------------------------------------------------------------------
markers.int <- FindAllMarkers(bm.int, only.pos = TRUE, min.pct = 0.1)
DotPlot(bm.int, features = unique(markers.int %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC) %>% pull(gene))) +
  scale_color_gradient2(low = 'darkgrey', high = 'firebrick4', mid = 'grey') + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Figure 4D Barplot showing cluster proportions among groups -----------------------------------------------------
mat.ratio <- table(bm.int@meta.data[, c('Group', 'cluster')])
for (i in 1:3) {
  sum <- sum(mat.ratio[i, ])
  for (q in 1:17) {
    mat.ratio[i, q] <- mat.ratio[i, q]/sum
  }
}
mat.ratio <- round(mat.ratio, digits = 3)
ggtexttable(t(mat.ratio))
mat.ratio <- as.data.frame(mat.ratio)
ggplot(mat.ratio, aes(x = group, y = Freq)) + geom_bar(position = 'stack', stat = 'identity', aes(fill = cluster)) + scale_fill_manual(values = color.int.new) +
  theme(panel.background = element_blank(), axis.line = element_line(), axis.title.x = element_blank(), axis.text = element_text(size = 15), axis.title.y = element_text(size = 20))

### Figure S4D Expression of HSC-specific markers ---------------------------------------------------------
DefaultAssay(hsc) <- 'RNA'
VlnPlot(hsc, features = c('Hlf', 'Procr'), ncol = 2)

### Figure S5A Expressiong of homing capacity of HSCs -----------------------------------
genes.homing <- c('Tfpi', 'Cxcr4', 'Itga4', 'Itga5')
hsc.5 <- AddModuleScore(hsc.5, features = list(genes.homing), name = 'Homing_score', assay = 'RNA')
mat.homing <- hsc.5@meta.data[, 'Homing_score1', drop = FALSE]
mat.homing$Group <- hsc.5@meta.data$Group
ggplot(mat.homing, aes(x = Group, y = Homing_score1)) + geom_boxplot(aes(fill = Group), color = 'black', outlier.shape = NA) +
  geom_signif(comparisons = list(c('Ctrl', 'R'), c('R', 'CBD'), c('Ctrl', 'CBD')), step_increase = 0.12, y_position = c(0.7), map_signif_level = TRUE, textsize = 3) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA)) +
  scale_fill_manual(values = c('Ctrl' = 'forestgreen', 'R' = 'blue', 'CBD' = 'firebrick2'))

### Figure 4E Point plot for DEGs in individual clusters among groups ------------------------------------
DefaultAssay(bm.int) <- 'RNA'
Idents(bm.int) <- 'group'
degs <- lapply(unique(bm.int@meta.data$cluster), function(x) {
  FindMarkers(bm.int[, bm.int$cluster == x], ident.1 = 'Ctrl', ident.2 = 'R')
})
names(degs) <- unique(bm.int$cluster)
for (i in 1:17) {
  degs[[i]]$gene <- rownames(degs[[i]])
  degs[[i]]$cluster[degs[[i]]$avg_log2FC > 0] <- 'Ctrl'
  degs[[i]]$cluster[degs[[i]]$avg_log2FC < 0] <- 'R'
}

df <- data.frame()
for (i in 1:length(degs)) {
  dat <- degs[[i]]
  dat$celltype <- names(degs)[[i]]
  dat$gene <- rownames(dat)
  dat$change <- ifelse(dat$avg_log2FC > 0 & dat $ p_val < 0.05, 'Ctrl',
                       ifelse(dat$avg_log2FC < 0 & dat$p_val < 0.05, 'R', 'nochange'))
  df <- rbind(df, dat)
}

rownames(df) <- 1:nrow(df)
df <- df[df$change != 'nochange', ]

for (i in 1:17) {
  print(range(df[df$celltype %in% levels(bm.int$cluster)[i], ]$avg_log2FC))
}

df$celltype <- factor(df$celltype, levels = levels(bm.int$cluster))
dfbar <- data.frame(x = c(seq(1, 17, 1)),
                    y = c(1.6, 1.65, 0.9, 1.1, 1.55, 0.56, 0.94, 1.17, 1.9, 1, 2.18, 1.42, 0.99, 1.29, 1.17, 1.3, 1.53),
                    celltype = levels(bm.int$Cluster.new))
dfbar1 <- data.frame(x = c(seq(1, 17, 1)),
                     y = c(-3.15, -2.03, -1.05, -1.25, -4.25, -1.97, -3, -1.86, -1.49, -1.3, -1.97, -2.07, -3.26, -2.85, -2.5, -2.84, -4.89),
                     celltype = levels(df$celltype))
dfcol <- data.frame(celltype = levels(df$celltype),
                    y = 0, lable = levels(df$celltype))
table(df[, c('change', 'celltype')])
bar.cr <- as.data.frame(t(table(df[, c('change', 'celltype')])))
bar.cr$celltype <- factor(bar.cr$celltype, levels = bar.cr[bar.cr$change %in% 'Ctrl', ]$celltype[order(bar.cr[bar.cr$change %in% 'Ctrl', ]$Freq, decreasing = F)])
bar.cr$change <- factor(bar.cr$change, levels = c('R', 'Ctrl'))
ggplot(bar.cr, aes(x = celltype, y = Freq)) + geom_bar(aes(fill = change), width = 0.75, stat = 'identity', position = position_dodge(0.75)) +
  scale_fill_manual(values = c('blue', 'forestgreen')) + coord_flip() +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5))+
  labs(x = 'Frequency', y = '', title = 'Number of DEGs')


df$celltype <- factor(df$celltype, levels = levels(bm.int$cluster))
dfnum <- data.frame(celltype = rep(levels(df$celltype), 2),
                    change = c(rep('Ctrl', 17), rep('R', 17)),
                    y = c(rep(2.4, 17), rep(-5.2, 17)),
                    #y = c(c(1.1, 1.4, 1.1, 0.5, 0.6, 0.5, 1.6, 0.7) + 0.8, c(-1.63, -1.25, -1.2, -0.8, -0.91, -0.85, -1.1, -1.2) - 0.8),
                    num = c(382, 431, 34, 33, 268, 22, 93, 223, 789, 31, 52, 305, 61, 77, 287, 623, 318,
                            601, 419, 40, 51, 181, 44, 101, 154, 526, 62, 292, 299, 103, 114, 220, 420, 396 ))

dfbar$celltype <- factor(dfbar$celltype, levels = levels(bm.int$cluster))
dfbar1$celltype <- factor(dfbar1$celltype, levels = levels(bm.int$cluster))

ggplot() + geom_col(data = dfbar, aes(x = celltype, y = y, fill = celltype), alpha = 0.2) + scale_fill_manual(values = color.int.new) +
  geom_col(data = dfbar1, aes(x = celltype, y = y, fill = celltype), alpha = 0.2) + 
  geom_tile(data = dfcol, aes(x = celltype, y = y, fill = celltype), height = 0.3, color = 'white') +
  geom_jitter(data = df, aes(x = celltype, y = avg_log2FC, color = change), width = 0.2, size = 0.3) + 
  scale_color_manual(values = c('forestgreen', 'blue')) +
  #geom_text_repel(data = topn_df, aes(x = celltype, y = avg_log2FC, label = gene), size = 3) +
  theme(panel.background = element_blank(), axis.line.y = element_line(), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x = '', y = 'avg_log2FC') +
  geom_text(data = dfnum, aes(x = celltype, y = y, label = num, color = change), size = 4)

Idents(bm.int) <- 'group'
degs2 <- lapply(unique(bm.int@meta.data$cluster), function(x) {
  FindMarkers(bm.int[, bm.int$cluster == x], ident.1 = 'CBD', ident.2 = 'R')
})
names(degs2) <- unique(bm.int$cluster)
for (i in 1:17) {
  degs2[[i]]$gene <- rownames(degs2[[i]])
  degs2[[i]]$cluster[degs2[[i]]$avg_log2FC > 0] <- 'CBD'
  degs2[[i]]$cluster[degs2[[i]]$avg_log2FC < 0] <- 'R'
}

df2 <- data.frame()
for (i in 1:length(degs2)) {
  dat <- degs2[[i]]
  dat$celltype <- names(degs2)[[i]]
  dat$gene <- rownames(dat)
  dat$change <- ifelse(dat$avg_log2FC > 0 & dat $ p_val < 0.05, 'CBD',
                       ifelse(dat$avg_log2FC < 0 & dat$p_val < 0.05, 'R', 'nochange'))
  df2 <- rbind(df2, dat)
}

rownames(df2) <- 1:nrow(df2)
df2 <- df2[df2$change != 'nochange', ]

for (i in 1:17) {
  print(range(df2[df2$celltype %in% levels(bm.int$cluster)[i], ]$avg_log2FC))
}

df2$celltype <- factor(df2$celltype, levels = levels(bm.int$cluster))
dfbar2 <- data.frame(x = c(seq(1, 17, 1)),
                     y = c(1.07, 1.26, 0.8, 1.32, 1.96, 0.56, 0.89, 0.76, 1.45, 0.72, 1.5, 0.67, 0.76, 1.05, 0.72, 0.85, 0.95),
                     celltype = levels(bm.int$Cluster.new))
dfbar1 <- data.frame(x = c(seq(1, 17, 1)),
                     y = c(-1.13, -1, -1.35, -1.1, -1.98, -2.02, -3, -1.99, -0.94, -0.89, -1.18, -2.45, -2.75, -2.4, -1.84, -2.23, -2.2),
                     celltype = levels(df2$celltype))
dfcol2 <- data.frame(celltype = levels(df2$celltype),
                     y = 0, lable = levels(df2$celltype))
table(df2[, c('change', 'celltype')])
bar.cr2 <- as.data.frame(t(table(df2[, c('change', 'celltype')])))
bar.cr2$celltype <- factor(bar.cr2$celltype, levels = bar.cr2[bar.cr2$change %in% 'CBD', ]$celltype[order(bar.cr2[bar.cr2$change %in% 'CBD', ]$Freq, decreasing = F)])
bar.cr2$change <- factor(bar.cr2$change, levels = c('R', 'CBD'))
ggplot(bar.cr2, aes(x = celltype, y = Freq)) + geom_bar(aes(fill = change), width = 0.75, stat = 'identity', position = position_dodge(0.75)) +
  scale_fill_manual(values = c('blue', 'firebrick3')) + coord_flip() +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5))+
  labs(x = 'Frequency', y = '', title = 'Number of DEGs')


df$celltype <- factor(df$celltype, levels = levels(bm.int$Cluster.new))
dfnum <- data.frame(celltype = rep(levels(df$celltype), 2),
                    change = c(rep('CBD', 17), rep('R', 17)),
                    y = c(rep(2.2, 17), rep(-3.2, 17)),
                    #y = c(c(1.1, 1.4, 1.1, 0.5, 0.6, 0.5, 1.6, 0.7) + 0.8, c(-1.63, -1.25, -1.2, -0.8, -0.91, -0.85, -1.1, -1.2) - 0.8),
                    num = c(169, 14, 8, 12, 153, 7, 42, 83, 41, 9, 23, 133, 36, 43, 95, 168, 111,
                            177, 24, 26, 19, 277, 35, 74, 70, 35, 22, 398, 168, 48, 65, 81, 143, 182))

dfbar$celltype <- factor(dfbar$celltype, levels = levels(bm.int$Cluster.new))
dfbar1$celltype <- factor(dfbar1$celltype, levels = levels(bm.int$Cluster.new))

ggplot() + geom_col(data = dfbar, aes(x = celltype, y = y, fill = celltype), alpha = 0.2) + scale_fill_manual(values = color.int.new) +
  geom_col(data = dfbar1, aes(x = celltype, y = y, fill = celltype), alpha = 0.2) + 
  geom_tile(data = dfcol, aes(x = celltype, y = y, fill = celltype), height = 0.3, color = 'white') +
  geom_jitter(data = df, aes(x = celltype, y = avg_log2FC, color = change), width = 0.2, size = 0.3) + 
  scale_color_manual(values = c('firebrick3', 'blue')) +
  #geom_text_repel(data = topn_df, aes(x = celltype, y = avg_log2FC, label = gene), size = 3) +
  theme(panel.background = element_blank(), axis.line.y = element_line(), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x = '', y = 'avg_log2FC') +
  geom_text(data = dfnum, aes(x = celltype, y = y, label = num, color = change), size = 4)

### Figure 4F Barplot showing the number of DEGs -------------------------------------------
bar.cr <- as.data.frame(t(table(df[, c('change', 'celltype')])))
bar.cr$celltype <- factor(bar.cr$celltype, levels = bar.cr[bar.cr$change %in% 'Ctrl', ]$celltype[order(bar.cr[bar.cr$change %in% 'Ctrl', ]$Freq, decreasing = F)])
bar.cr$change <- factor(bar.cr$change, levels = c('R', 'Ctrl'))
ggplot(bar.cr, aes(x = celltype, y = Freq)) + geom_bar(aes(fill = change), width = 0.75, stat = 'identity', position = position_dodge(0.75)) +
  scale_fill_manual(values = c('blue', 'forestgreen')) + coord_flip() +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5))+
  labs(x = 'Frequency', y = '', title = 'Number of DEGs')

### Figure S4E Barplot showing the number of DEGs -------------------------------------------
bar.cr2 <- as.data.frame(t(table(df2[, c('change', 'celltype')])))
bar.cr2$celltype <- factor(bar.cr2$celltype, levels = bar.cr2[bar.cr2$change %in% 'CBD', ]$celltype[order(bar.cr2[bar.cr2$change %in% 'CBD', ]$Freq, decreasing = F)])
bar.cr2$change <- factor(bar.cr2$change, levels = c('R', 'CBD'))
ggplot(bar.cr2, aes(x = celltype, y = Freq)) + geom_bar(aes(fill = change), width = 0.75, stat = 'identity', position = position_dodge(0.75)) +
  scale_fill_manual(values = c('blue', 'firebrick3')) + coord_flip() +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5))+
  labs(x = 'Frequency', y = '', title = 'Number of DEGs')

### Figure 4G-H Presentation of GO terms enriched in all clusters ------------------------------------------
Idents(bm.int) <- 'group'
markers.group <- FindMarkers(bm.int, ident.1 = 'Ctrl', ident.2 = 'R', logfc.threshold = log2(1.1), return.thresh = 0.05)

go.r <- enrichGO(bitr(markers.group$gene[markers.group$cluster %in% 'R'], OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                    OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
kegg.r <- enrichKEGG(bitr(markers.group$gene[markers.group$cluster %in% 'R'], OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                        organism = 'mmu')
kegg.r <- setReadable(kegg.ctrl, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID')

term.r <- c('response to oxidative stress', 'apoptotic mitochondrial changes', 'intrinsic apoptotic signaling pathway', 'intrinsic apoptotic signaling pathway in response to DNA damage',
            'regulation of mitotic cell cycle')
cnetplot(go.r, showCategory = term.r, colorEdge = T, color_gene = 'deeppink4', color_category = 'darkslateblue',
         color_edge = c('brown4', 'darkgoldenrod3', 'darkolivegreen', 'darkslateblue'))

markers.group.2 <- FindMarkers(bm.int, ident.1 = 'CBD', ident.2 = 'R', logfc.threshold = log2(1.1), return.thresh = 0.05)
go.cbd <- enrichGO(bitr(markers.group.2$gene[markers.group.2$cluster %in% 'CBD'], OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                 OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
term.cbd <- c('leukocyte migration', 'myeloid leukocyte migration', 'neutrophil migration', 'granulocyte migration', 'cell chemotaxis',
              'regulation of leukocyte differentiation', 'regulation of hemopoiesis', 'regulation of lymphocyte differentiation', 'T cell differentiation',
              'mononuclear cell differentiation', 'regulation of inflammatory response', 'humoral immune response', 'negative regulation of immune system process',
              'neutrophil mediated immunity', 'negative regulation of cytokine production')
results.cbd <- go.cbd@result[go.cbd@result$Description %in% term.cbd, ]
results.cbd$Description <- factor(results.cbd$Description, levels = results.cbd$Description[order(results.cbd$pvalue, decreasing = T)])
results.cbd$logp <- -1*log10(results.cbd$pvalue)
ggplot(results.cbd, aes(x = logp, y = Description)) + geom_bar(stat = 'identity', fill = 'firebrick3') +
  theme(axis.line = element_line(), panel.background = element_blank())

### Figure 5A HSC reUMAP-----------------------------------------------------------------------
hsc <- bm.int[, bm.int$cluster %in% 'HSC']
DefaultAssay(hsc) <- 'integrated'
hsc <- RunPCA(hsc, features = rownames(hsc))
ElbowPlot(hsc, ndims = 50)

hsc <- FindNeighbors(hsc, dims = 1:10)
hsc <- RunUMAP(hsc, dims = 1:10)
hsc <- FindClusters(hsc, resolution = 2)
CellDimPlot(hsc.small, group.by = 'group', palcolor =  c('forestgreen', 'blue', 'firebrick3'), pt.size = 1)

### Figure 5B Dotplot for expression of DEGs in HSCs among groups--------------------------------------------------------------
DefaultAssay(hsc) <- 'RNA'
hsc <- ScaleData(hsc, features = rownames(hsc))
Idents(hsc) <- 'group'
table(hsc@meta.data$group)
markers.hsc <- FindAllMarkers(hsc, only.pos = TRUE, return.thresh = 0.05)
table(markers.hsc$cluster)
DotPlot(hsc, features = unique(markers.hsc %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC) %>% pull(gene)), cols = c('darkgrey', 'darkred')) + coord_flip()

### Figure 5C and S5D Enrichplot for GO terms enrichde in HSCs among groups-------------------------------------------------------------
hsc.cbd.r <- hsc[, hsc$group %in% c('R', 'CBD')]
hsc.cbd.r <- ScaleData(hsc.cbd.r, features = rownames(hsc.cbd.r))
hsc.cbd.r <- RunDEtest(hsc.cbd.r, group_by = 'group', fc.threshold = 2^0.25, only.pos = T)
VolcanoPlot(hsc.cbd.r, group_by = 'group')

degs.cbd.r <- hsc.cbd.r@tools$DEtest_Group$AllMarkers_wilcox
degs.cbd.r <- degs.cbd.r[with(degs.cbd.r, avg_log2FC > 1 & p_val_adj < 0.05), ]
hsc.cbd.r <- RunEnrichment(srt = hsc.cbd.r, group_by = 'group', db = 'GO_BP', species = 'Mus_musculus', DE_threshold = 'avg_log2FC > 0.25 & p_val < 0.05')
EnrichmentPlot(hsc.cbd.r, group_by = 'group', plot_type = 'enrichmap', palcolor = c('darkorchid4', 'deeppink3', 'darkgreen'))

hsc.ctrl.r <- hsc[, hsc$group %in% c('R', 'Ctrl')]
hsc.ctrl.r <- ScaleData(hsc.ctrl.r, features = rownames(hsc.ctrl.r))
hsc.ctrl.r <- RunDEtest(hsc.ctrl.r, group_by = 'group', fc.threshold = 2^0.25, only.pos = T)
VolcanoPlot(hsc.ctrl.r, group_by = 'group')

degs.ctrl.r <- hsc.ctrl.r@tools$DEtest_Group$AllMarkers_wilcox
degs.ctrl.r <- degs.ctrl.r[with(degs.ctrl.r, avg_log2FC > 1 & p_val_adj < 0.05), ]
hsc.ctrl.r <- RunEnrichment(srt = hsc.ctrl.r, group_by = 'group', db = 'GO_BP', species = 'Mus_musculus', DE_threshold = 'avg_log2FC > 0.25 & p_val < 0.05')
EnrichmentPlot(hsc.ctrl.r, group_by = 'group', plot_type = 'enrichmap', palcolor = c('darkorchid4', 'deeppink3', 'darkgreen'))

### Figure 5D-E Barplot showing representative GO terms enriched in HSCs among groups -----------------------------------------------------------
hsc.ctrl.r <- hsc[, hsc$group %in% c('Ctrl', 'R')]
markers.ctrl.r <- FindAllMarkers(hsc.ctrl.r, only.pos = T)
genes.ctrl <- markers.ctrl.r$gene[markers.ctrl.r$cluster %in% 'Ctrl']
go.hsc.ctrl <- enrichGO(bitr(genes.ctrl, OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                        OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
go.display.hsc.ctrl <- c('stem cell population maintenance', 'stem cell differentiation', 'lymphocyte differentiation', 'myeloid cell differentiation', 
                         'BMP signaling pathway', 'canonical Wnt signaling pathway', 'cell-cell signaling by wnt', 'negative regulation of cell cycle', 
                         'homeostasis of number of cells', 'regulation of hemopoiesis')
go.hsc.ctrl.result <- go.hsc.ctrl@result[go.hsc.ctrl@result$Description %in% go.display.hsc.ctrl, ]
go.hsc.ctrl.result$logp <- -1*log10(go.hsc.ctrl.result$p.adjust)
go.hsc.ctrl.result$cluster <- 'Ctrl'

genes.r <- markers.ctrl.r$gene[markers.ctrl.r$cluster %in% 'R']
go.hsc.r <- enrichGO(bitr(genes.R, OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                     OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
View(go.hsc.r@result)
go.display.hsc.r <- c('intrinsic apoptotic signaling pathway', 'intrinsic apoptotic signaling pathway in response to DNA damage', 'intrinsic apoptotic signaling pathway by p53 class mediator',
                      'mitotic cell cycle phase transition', 'response to oxidative stress', 'activation of immune response',
                      'regulation of apoptotic signaling pathway', 'negative regulation of cell migration', 'regulation of inflammatory response',
                      'cellular respiration')
go.hsc.r.result <- go.hsc.r@result[go.hsc.r@result$Description %in% go.display.hsc.r, ]
go.hsc.r.result$logp <- log10(go.hsc.r.result$p.adjust)
go.hsc.r.result$cluster <- 'R'

go.hsc.comb <- rbind(go.hsc.ctrl.result, go.hsc.r.result)

go.hsc.comb$Description <- factor(go.hsc.comb$Description,
                                  levels = go.hsc.comb$Description[order(go.hsc.comb$logp, decreasing = F)])
ggplot(go.hsc.comb, aes(x = logp, y = Description)) + geom_bar(stat = 'identity', aes(fill = cluster)) + scale_fill_manual(values = c('forestgreen', 'blue')) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5)) +
  labs(x = '-log10(p.adj)', y = '', title = 'HSC')


hsc.cbd.r <- hsc[, hsc$group %in% c('CBD', 'R')]
markers.cbd.r <- FindAllMarkers(hsc.cbd.r, only.pos = T)

genes.cbd <- markers.cbd.r$gene[markers.cbd.r$cluster %in% 'CBD']
go.hsc.cbd <- enrichGO(bitr(genes.cbd, OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                        OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
go.display.hsc.cbd <- c('stem cell population maintenance', 'stem cell proliferation', 'stem cell differentiation', 
                        'mitotic cell cycle phase transition', 'positive regulation of cell cycle', 'cell cycle G1/S phase transition',
                        'myeloid cell differentiation', 'lymphocyte differentiation', 'stem cell differentiation', 'regulation of leukocyte differentiation')
go.hsc.cbd.result <- go.hsc.cbd@result[go.hsc.cbd@result$Description %in% go.display.hsc.cbd, ]
go.hsc.cbd.result$logp <- -1*log10(go.hsc.cbd.result$p.adjust)
go.hsc.cbd.result$cluster <- 'CBD'

genes.r.2 <- markers.cbd.r$gene[markers.cbd.r$cluster %in% 'R']
go.hsc.r.2 <- enrichGO(bitr(genes.R, OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID') %>% pull(ENTREZID),
                     OrgDb = org.Mm.eg.db, ont = 'BP', readable = T)
go.display.hsc.r.2 <- c('intrinsic apoptotic signaling pathway', 'intrinsic apoptotic signaling pathway in response to DNA damage', 'intrinsic apoptotic signaling pathway by p53 class mediator',
                      'mitotic cell cycle phase transition', 'response to oxidative stress', 'activation of immune response',
                      'regulation of apoptotic signaling pathway', 'negative regulation of cell migration', 'regulation of inflammatory response',
                      'cellular respiration')
go.hsc.r.result.2 <- go.hsc.r.2@result[go.hsc.r.2@result$Description %in% go.display.hsc.r.2, ]
go.hsc.r.result.2$logp <- log10(go.hsc.r.result.2$p.adjust)
go.hsc.r.result.2$cluster <- 'R'

go.hsc.comb <- rbind(go.hsc.cbd.result, go.hsc.r.result.2)

go.hsc.comb$Description <- factor(go.hsc.comb$Description,
                                  levels = go.hsc.comb$Description[order(go.hsc.comb$logp, decreasing = F)])
ggplot(go.hsc.comb, aes(x = logp, y = Description)) + geom_bar(stat = 'identity', aes(fill = cluster)) + scale_fill_manual(values = c('blue', 'firebrick3')) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5)) +
  labs(x = '-log10(p.adj)', y = '', title = 'HSC')

### Figure 5F Heatmap showing expression of genes involved in indicated GO BP in HSC-----------------------------------------------------------
genes.wnt <- strsplit(go.hsc.cbd@result$geneID[go.hsc.cbd@result$Description %in% 'Wnt signaling pathway'], split = '/')[[1]]
genes.bmp <- strsplit(go.hsc.cbd@result$geneID[go.hsc.cbd@result$Description %in% 'BMP signaling pathway'], split = '/')[[1]]
genes.negapop <- strsplit(go.hsc.cbd@result$geneID[go.hsc.cbd@result$Description %in% 'negative regulation of apoptotic signaling pathway'], split = '/')[[1]]
genes.intrinsicapop <- strsplit(go.hsc.r@result$geneID[go.hsc.r@result$Description %in% 'intrinsic apoptotic signaling pathway'], split = '/')[[1]]
genes.responseoxidative <- strsplit(go.hsc.r@result$geneID[go.hsc.r@result$Description %in% 'response to oxidative stress'], split = '/')[[1]]
genes.oxidativephos <- strsplit(go.hsc.r@result$geneID[go.hsc.r@result$Description %in% 'oxidative phosphorylation'], split = '/')[[1]]

data.hsc.feature <- as.data.frame(hsc@assays$RNA@scale.data[c(genes.wnt, genes.bmp, genes.negapop, genes.intrinsicapop, 
                                                              genes.responseoxidative, genes.oxidativephos), ])
data.hsc.feature$Ctrl <- apply(data.hsc.feature[, rownames(hsc@meta.data)[hsc$group %in% 'Ctrl']], 1, mean)
data.hsc.feature$R <- apply(data.hsc.feature[, rownames(hsc@meta.data)[hsc$group %in% 'R']], 1, mean)
data.hsc.feature$CBD <- apply(data.hsc.feature[, rownames(hsc@meta.data)[hsc$group %in% 'CBD']], 1, mean)
dim(data.hsc.feature)
data.hsc.feature <- data.hsc.feature[, 2463:2465]
data.hsc.feature <- scale(t(data.hsc.feature))

annot.features <- data.frame(row.names = colnames(data.hsc.feature),
                             annot = c(rep('Wnt', 22), rep('BMP', 11), rep('negApop', 14), rep('Apoptotic', 16), rep('responseOxi', 21), rep('OxPhos', 12)))
color.annot <- RColorBrewer::brewer.pal('Paired', n = 6)
names(color.annot) <- c('Apoptotic', 'BMP', 'negApop', 'OxPhos', 'responseOxi', 'Wnt')
pheatmap(data.hsc.feature, cluster_cols = T, cluster_rows = F, annotation_col = annot.features, cutree_cols = 5, annotation_colors = list(annot = color.annot),
         color = colorRampPalette(c('slateblue4', 'grey', 'white', 'red', 'firebrick4'))(40))

### Figure S5C-D Volcano plot for DEGs and top10 GO terms enriched by HSCs from Control and IR+CBD -------------------------------------------
markers.hsc.ctrl.r$logp <- -log10(markers.hsc.ctrl.r$p_val)
markers.hsc.ctrl.r$avg_log2FC[markers.hsc.ctrl.r$cluster %in% 'R'] <- -1*markers.hsc.ctrl.r$avg_log2FC[markers.hsc.ctrl.r$cluster %in% 'R']
ggplot(markers.hsc.ctrl.r, aes(x = avg_log2FC, y = logp)) + geom_point(aes(color = cluster)) +
  theme(axis.line = element_line(), panel.background = element_blank()) +
  scale_color_manual(values = c('blue', 'forestgreen')) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(x = 'Average log2(fold change)', y = '-log10(P value)')

compareGO(markers.hsc.ctrl.r, show.cat = 10, width.y = 90, ont = 'CC')

### Figure S5E Expression of Atf2 in clusters -----------------------------------------------------------------
FeatureStatPlot(bm.int, stat.by = 'Atf2', group.by = 'cluster', palcolor = color.int)

### Fgiure 6A Venn plot showing the number of pairwise DEGs enriched by HSC from Control and IR+CBD versus IR+DMSO--------------------------------------------------
library(ggvenn)
deg.venn <- list(Ctrl = c(markers.hsc.ctrl.r$gene[markers.hsc.ctrl.r$cluster %in% 'Ctrl'], markers.hsc.cbd.ctrl$gene[markers.hsc.cbd.ctrl$cluster %in% 'Ctrl']),
                 CBD = c(markers.hsc.cbd.r$gene[markers.hsc.cbd.r$cluster %in% 'CBD'], markers.hsc.cbd.ctrl$gene[markers.hsc.cbd.ctrl$cluster %in% 'CBD']),
                 R = c(markers.hsc.ctrl.r$gene[markers.hsc.ctrl.r$cluster %in% 'R'], markers.hsc.cbd.r$gene[markers.hsc.cbd.r$cluster %in% 'R']))
ggvenn(deg.venn, stroke_color = 'white', show_percentage = FALSE, fill_color = c('forestgreen', 'firebrick2', 'blue'), set_name_color = c('forestgreen', 'firebrick2', 'blue'),
       set_name_size = 4)

### Figure 6B Heatmap showing the expression of 73 genes indicated in Figure 5A -----------------------------------------------------------------
gene.key.1 <- intersect(markers.hsc.ctrl.r$gene[markers.hsc.ctrl.r$cluster %in% 'Ctrl'],
                        markers.hsc.cbd.r$gene[markers.hsc.cbd.r$cluster %in% 'CBD'])
DotPlot(hsc, features = gene.key.1) + coord_flip()
mat.inter.up <- AverageExpression(hsc, features = gene.key.1, assays = 'RNA')
mat.inter.up <- t(scale(t(mat.inter.up$RNA)))
annot.inter.up <- data.frame(row.names = c('Ctrl', 'R', 'CBD'), Group = c('Ctrl', 'R', 'CBD'))
annot.inter.up$Group <- factor(annot.inter.up, levels = c('Ctrl', 'R', 'CBD'))
pheatmap(mat.inter.up, cluster_rows = TRUE, angle_col = 90, cluster_cols = FALSE, annotation_col = annot.inter.up, annotation_colors = list(Group = c('forestgreen', 'blue', 'firebrick2')),
         color = colorRampPalette(c('slateblue4', 'grey', 'white', 'red', 'firebrick4'))(40))

### Figure 5L-M Barplot showing Wnt-related pathways enriched and expression of Wnt-related genes in HSCs ------------------------------------------------
mat.cbd <- go.cbd@result[go.cbd@result$Description %in% c('Wnt signaling pathway', 'canonical Wnt signaling pathway', 'cell-cell signaling by wnt',
                                                          'regulation of Wnt signaling pathway', 'regulation of canonical Wnt signaling pathway',
                                                          'negative regulation of apoptotic signaling pathway', 'regulation of apoptotic signaling pathway'
), ]
mat.cbd$logp <- -1*log10(mat.cbd$pvalue)
mat.cbd$Description <- factor(mat.cbd$Description, levels = mat.cbd$Description[order(mat.cbd$pvalue, decreasing = T)])
ggplot(mat.cbd, aes(x = logp, y = Description)) + geom_bar(fill = 'firebrick3', stat = 'identity') +
  theme(panel.background = element_blank(), axis.line = element_line())

FeatureStatPlot(hsc.5, stat.by = c('Lrp6', 'Emd', 'Pten', 'Atp6ap2', 'Kpna1', 'Ppm1b'), group.by = 'Group', palcolor = c('forestgreen', 'blue', 'firebrick2'),
                comparisons = list(c('Ctrl', 'R'), c('R', 'CBD'), c('Ctrl', 'CBD')), plot_type = 'box', stack = T, flip = T)











