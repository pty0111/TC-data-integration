setwd("~/CCR7_DC/MLN_RORgt_MHCII_multiome/")

library(ArchR) # v1.0.1
library(Seurat) # v4.0.4
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

plot.clusters <- function(sro = NULL, idx = NULL, nudge_y = 0.5, cl.name = "Clusters",
                          vis = sro@reductions$umap@cell.embeddings,
                          clusters = sro$Clusters.annot, clusters.label = sro$Clusters,
                          col = rna.pal, pref.C = T, labels = T, line = F){
  gp <- ggplot() + labs(color = cl.name) + theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(is.null(idx)) idx <- rep(T, nrow(vis))
  gp <- gp + geom_point(aes(x = vis[idx, 1], y = vis[idx, 2], color = clusters[idx]), size = 2, alpha = 0.7)
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters.label))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[idx & clusters == c, 1])
      vis.cl[c, 2] <- median(vis[idx & clusters == c, 2])
    }
    gp + geom_label(aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
                    size = 5, show.legend = F, nudge_y = nudge_y)
  } else gp
}
plot.continuous.value <- function(vis, val, val.name, scale.color, point.size = 3){
  ggplot() +
    geom_point(aes(x = vis[, 1], y = vis[, 2], color = val), size = point.size, alpha = 0.8) +
    scale.color + labs(color = val.name) + theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
}
plot.violin <- function(df, col, ...){
  ggplot(df, aes(...)) +
    # geom_violin(show.legend = F, alpha = 0.2) +
    geom_boxplot(show.legend = F) +
    theme_classic() + scale_color_manual(values = col) + scale_fill_manual(values = col)
}

# rna.pal <- c("#e051bc", "#D790FF", "#BC23FF", "#72418F", "#3A2465", "#7daedb", "#0976de",
#              "#003a85", "#17e81b", "#77d9ae", "#14a38e", "#05a15c", "#4e7506", "#0d6930",
#              "#405933", "#0d3b0e", "#c2b206", "#d67c22", "#e32712", "#A05837", "#5B4534")
# names(rna.pal) <- names(rna.annot)
# rna.pal2 <- rna.pal; names(rna.pal2) <- 1:21
#
# atac.pal <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
#               "#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6",
#               "#63FFAC","#B79762","#004D43")
# names(atac.pal) <- paste0("C", 1:length(atac.pal))
# save(rna.annot, rna.pal, rna.pal2, atac.pal, file = "palette.RData")
load("palette.RData")

pdf("Seurat/plots/QC/cluster-QC.pdf", width = 20, height = 20)
plot.all.QC(seurat.obj, ident = "Clusters", col = rna.pal2)
dev.off()

cl <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1)
cl$Clusters.annot <- factor(cl$Clusters.annot, levels = names(rna.annot))
cl$Clusters <- factor(cl$Clusters, levels = 1:21)

plot.clusters(clusters = cl$Clusters.annot, clusters.label = cl$Clusters,
              vis = read.csv("Seurat/results/UMAP.csv", header = 1, row.names = 1, col.names = c("cell", "UMAP_1", "UMAP_2"))) %>%
  ggsave(filename = "Seurat/plots/clusters.pdf", width = 15, height = 12)
plot.clusters(clusters = cl$Clusters.annot, clusters.label = cl$Clusters,
              vis = read.csv("Seurat/results/FDL.csv", header = 1, row.names = 1, col.names = c("cell", "FDL_1", "FDL_2"))) %>%
  ggsave(filename = "Seurat/plots/clusters-FDL.pdf", width = 15, height = 15)

idx <- rownames(subset(cl, Clusters %in% 1:16))
plot.clusters(clusters = cl[idx, "Clusters.annot"], clusters.label = cl[idx, "Clusters"], col = rna.pal[1:16],
              vis = read.csv("Seurat/results/C1-16_UMAP.csv", header = 1, row.names = 1,
                             col.names = c("cell", "UMAP_1", "UMAP_2"))[idx, ]) %>%
  ggsave(filename = "Seurat/plots/clusters_C1-16.pdf", width = 15, height = 12)
plot.clusters(clusters = cl[idx, "Clusters.annot"], clusters.label = cl[idx, "Clusters"], col = rna.pal[1:16],
              vis = read.csv("Seurat/results/C1-16_FDL.csv", header = 1, row.names = 1,
                             col.names = c("cell", "FDL_1", "FDL_2"))[idx, ]) %>%
  ggsave(filename = "Seurat/plots/clusters-FDL_C1-16.pdf", width = 15, height = 15)

joint.md <- read.csv("Joint/results/meta-data.csv", header = 1, row.names = 1)
joint.md$wsnn_res.1.2 <- factor(joint.md$wsnn_res.1.2); joint.md$wsnn_res.2 <- factor(joint.md$wsnn_res.2)
plot_grid(
  plot.clusters(clusters = joint.md[rownames(joint.umap), "wsnn_res.1.2"],
                clusters.label = joint.md[rownames(joint.umap), "wsnn_res.1.2"],
                vis = joint.umap, col = NULL, line = T, cl.name = "Joint Clusters\nres=1.2"),
  plot.clusters(clusters = joint.md[rownames(joint.umap), "wsnn_res.2"],
                clusters.label = joint.md[rownames(joint.umap), "wsnn_res.2"],
                vis = joint.umap, col = NULL, line = T, cl.name = "Joint Clusters\nres=2.0"),
  ncol = 1) %>% ggsave(filename = "Joint/plots/Joint-clusters.pdf", width = 15, height = 30)

atac.UMAP <- read.csv("ArchR/ArchROutput_SeuratC1-16/Embeddings/UMAP.csv", header = 1, row.names = 1, col.names = c("cell", "UMAP_1", "UMAP_2"))
atac.md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", header = 1, row.names = 1)
atac.md$Clusters <- factor(atac.md$Clusters, levels = paste0("C", 1:13))
rna.md <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1)[gsub("^.+#", "", rownames(atac.md)), ]
pdf(file = "ArchR/ArchROutput_SeuratC1-16/Plots/UMAP-after_filtering_by_scRNA-seq.pdf", width = 15, height = 12)
plot.clusters(clusters = atac.md$Clusters, clusters.label = atac.md$Clusters,
              col = atac.pal, vis = atac.UMAP, pref.C = F)
plot.clusters(clusters = factor(rna.md$Clusters.annot, levels = names(rna.annot)[1:16]),
              clusters.label = factor(rna.md$Clusters, levels = 1:16),
              col = rna.pal[1:16], vis = atac.UMAP) +
  labs(color = "Clusters.scRNAseq")
dev.off()

pdf("Seurat/plots/CellCycleScores.pdf", width = 20, height = 10)
df <- cbind(
  read.csv("Seurat/results/meta-data.csv", row.names = 1, header = T),
  read.csv("Seurat/results/UMAP.csv", row.names = 1, header = T),
  read.csv("Seurat/results/FDL.csv", row.names = 1, header = T, col.names = c("cell", "FDL_1", "FDL_2"))
); df$Clusters <- factor(df$Clusters)
plot_grid(
  plot.violin(df, x = Clusters, y = S.Score, color = Clusters.annot, fill = Clusters.annot, col = rna.pal) + labs(y = "S Phase Score"),
  plot.violin(df, x = Clusters, y = G2M.Score, color = Clusters.annot, fill = Clusters.annot, col = rna.pal) + labs(y = "G2/M Phase Score"),
  ncol = 1
)
gp <- ggplot(df) + scale_color_gradient(high = "red", low = "grey") +
  theme_classic() + theme(axis.ticks = element_blank(), axis.text = element_blank())
plot_grid(
  gp + geom_point(aes(x = UMAP_1, y = UMAP_2, color = S.Score), alpha = (df$S.Score + 0.2) / 1.5, size = 2) + labs(color = "S Phase Score"),
  gp + geom_point(aes(x = UMAP_1, y = UMAP_2, color = G2M.Score), alpha = (df$G2M.Score + 0.2) / 2, size = 2) + labs(color = "G2/M Phase Score"),
  ncol = 2
)
plot_grid(
  gp + geom_point(aes(x = FDL_1, y = FDL_2, color = S.Score), alpha = (df$S.Score + 0.2) / 1.5, size = 2) + labs(color = "S Phase Score"),
  gp + geom_point(aes(x = FDL_1, y = FDL_2, color = G2M.Score), alpha = (df$G2M.Score + 0.2) / 2, size = 2) + labs(color = "G2/M Phase Score"),
  ncol = 2
)
dev.off()
df <- subset(df, Clusters %in% 1:16)
df$Clusters <- ifelse(df$Clusters %in% 9:16, "9-16", df$Clusters) %>% factor(levels = c(1:8, "9-16"))
plot_grid(
  plot.violin(df, x = Clusters, y = S.Score, color = Clusters, fill = Clusters,
              col = c(rna.pal2[1:8], `9-16` = as.character(rna.pal2[11]))) + labs(y = "S Phase Score"),
  plot.violin(df, x = Clusters, y = G2M.Score, color = Clusters, fill = Clusters,
              col = c(rna.pal2[1:8], `9-16` = as.character(rna.pal2[11]))) + labs(y = "G2/M Phase Score"),
  ncol = 1
) %>% ggsave(filename = "Seurat/plots/CellCycleScores-groupedLTis.pdf", width = 10, height = 10)

df <- read.csv("Seurat/results_before_filtering_by_scATAC-seq/UMAP.csv", header = 1, row.names = 1)
df$Clusters <- factor(read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1)[rownames(df), "Clusters"], levels = 1:21)
gp <- ggplot(mapping = aes(x = UMAP_1, y = UMAP_2, color = Clusters)) + theme_classic()
ggsave(
  plot = plot_grid(
    gp + geom_point(data = subset(df, is.na(Clusters)), color = "grey", size = 1)+
      geom_point(data = subset(df, !is.na(Clusters)), size = 2) + scale_color_manual(values = rna.pal2),
    gp + geom_point(data = subset(df, !is.na(Clusters)), color = "blue", size = 2) + scale_color_manual(values = rna.pal2) +
      geom_point(data = subset(df, is.na(Clusters)), color = "grey", size = 1),
    ncol = 1),
  filename = "Seurat/plots_before_filtering_by_scATAC-seq/clusters_after_filtering_by_scATAC-seq.pdf", width = 15, height = 20
)

atac.sample <- "3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19"
df <- read.csv("Seurat/results_before_filtering_by_scATAC-seq/meta-data.csv", row.names = 1)
df <- cbind(
  rna.before = df, rna.after = read.csv("Seurat/results/meta-data.csv", row.names = 1)[rownames(df), ],
  atac = read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", row.names = 1)[paste(atac.sample, rownames(df), sep = "#"), ],
  joint = joint.md[rownames(df), c("wsnn_res.1.2", "wsnn_res.2")]
)
df$atac.Clusters <- factor(as.integer(gsub("^C", "", df$atac.Clusters)))
table(df[, c("rna.before.Clusters", "rna.after.Clusters")])

colnames(df)[10] <- "rna.Clusters"
table(df[, c("rna.Clusters", "atac.Clusters")])
table(df[, c("rna.Clusters", "joint.wsnn_res.1.2")])
table(df[, c("rna.Clusters", "joint.wsnn_res.2")])
table(df[, c("atac.Clusters", "joint.wsnn_res.1.2")])
table(df[, c("atac.Clusters", "joint.wsnn_res.2")])

genes <- subset(seurat.obj@assays$RNA@meta.features, toupper(Symbol) %in% toupper(c("Rorc", "Ptprc", "RORa", "Rag1")))
genes$ID <- rownames(genes); rownames(genes) <- genes$Symbol
pdf("Seurat/plots/violin-plots.pdf", width = 15, height = 7)
df <- cbind(
  seurat.obj@meta.data,
  Ptprc = seurat.obj@assays$RNA@data[genes["Ptprc", "ID"], ],
  Rag1  = seurat.obj@assays$RNA@data[genes["Rag1",  "ID"], ],
  Rorc  = seurat.obj@assays$RNA@data[genes["Rorc",  "ID"], ],
  Rora  = seurat.obj@assays$RNA@data[genes["Rora",  "ID"], ]
) %>% subset(Clusters %in% 1:16)
df$Clusters <- ifelse(df$Clusters %in% 9:16, "9-16", df$Clusters)
df$Clusters.annot <- gsub("LTi.+", "LTi", df$Clusters.annot)
plot.violin(df, x = Clusters, y = Ptprc, color = Clusters.annot,
            fill = Clusters.annot, col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
plot.violin(df, x = Clusters, y = Rag1, color = Clusters.annot,
            fill = Clusters.annot, col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
plot.violin(df, x = Clusters, y = Rorc, color = Clusters.annot,
            fill = Clusters.annot, col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
plot.violin(df, x = Clusters, y = Rora, color = Clusters.annot,
            fill = Clusters.annot, col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
dev.off()

genes <- genes[c(2, 4), ]
seurat.obj$Clusters2 <- ifelse(seurat.obj$Clusters %in% 9:16, "9-16", seurat.obj$Clusters)
ggsave(filename = "~/Desktop/dotplot.pdf", width = 5, height = 5,
       DotPlot(subset(seurat.obj, Clusters %in% 1:16), group.by = "Clusters2",
               assay = "RNA", features = genes$ID,
               col.min = -1, col.max = 1, scale = 15) +
         labs(x = "", y = "Cluster", size = "Percent\nExpressed") +
         scale_x_discrete(breaks = genes$ID, labels = genes$Symbol) +
         scale_size_continuous(limits = c(0, 100), range = c(0.1, 8)) +
         scale_color_gradient2(low = "#f9eee8", high = "#5a1419", mid = "#e46c50", midpoint = 0) +
         guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3, title = "Average\nExpression")) +
         theme(axis.text.x.bottom = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
               axis.title.x = element_blank(),
               legend.position = "right", legend.justification = "center",
               legend.key.width = unit(0.5, "inch"), legend.box = "vertical",
               legend.title = element_text(vjust = 1, hjust = 1, size = 5),
               legend.text = element_text(size = 3), plot.margin = unit(rep(1, 4), "mm"))
)

df$Clusters.annot <- factor(gsub("LTi.+", "LTi", df$Clusters.annot),
                            levels = c(levels(df$Clusters.annot)[1:8], "LTi"))
ggsave(
  filename = "~/Desktop/Rorc-expression-violin-2.pdf", width = 5, height = 5,
  plot =
plot.violin(na.omit(df),
            x = Clusters.annot, y = Rorc,
            color = Clusters.annot, #fill = Clusters.annot,
            col = c(rna.pal[1:8], LTi = as.character(rna.pal[11]))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
)

sro <- subset(seurat.obj, Clusters %in% setdiff(1:16, c(1, 6, 7)))
sro$Clusters <- as.character(sro$Clusters); sro$Clusters[sro$Clusters %in% 9:16] <- "9-16"
# g.list <- list(
#   `Transcription-Factor` = c("Zbtb16", "Rora", "Tox", "Tox2", "Ets1", "Ikzf2",
#                              "Bach2", "Tbx21", "Aire", "Irf8", "Etv3", "Prdm16"),
#   `Cell-Surface` = c("Cxcr6", "Il7r", "Il2ra", "Il2rb", "Il23r", "Il22", "Cd4", "Tnfsf11",
#                      "Ncr1", "Ccr6", "Zbtb46", "Flt3", "Itgax", "Ccr7", "Ly75", "Marcks",
#                      "Mreg", "Socs2", "Dnase1l3", "Alcam", "Col17a1", "Il15"),
#   `Neuronal-Gene` = c("Gal", "Ncam1", "Nrxn1", "Nrxn2", "Nlgn2", "Nrg1"),
#   `Antigen-Presenting` = c("Cd74", "Ciita", "Nlrc5", "H2-DMa", "H2-DMb1",
#                            "H2-DMb2", "H2-M2", "Wdfy4", "Cd80", "Cd86", "Cd40"),
#   `TGFb-Gene` = c("Itgav", "Itgb8", "Nrros", "Lrrc32")
# )
gene.info <- sro@assays$RNA@meta.features
g.list <- data.frame(readxl::read_xlsx("Dotplot gene lists.xlsx"), check.names = F)
colnames(g.list) <- gsub("\\/", "-", colnames(g.list))
W <- c(7.1, 3, 3.5, 3.1, 3, 3.5) * 1.1; names(W) <- colnames(g.list)
H <- c(3.1, 2.5, 3, 3, 1.9, 2.7) * 0.9; names(H) <- colnames(g.list)
dp.list <- sapply(colnames(g.list), function(list.name){
  gene.list <- na.omit(g.list[, list.name])
  gene.list <- gene.info[match(toupper(gene.list), toupper(gene.info$Symbol)), ]
  gene.list <- gene.list[rownames(gene.list) %in% rownames(sro), ]
  if(!all(toupper(na.omit(g.list[, list.name])) %in% toupper(gene.list$Symbol))) print(list.name)
  ggsave(filename = paste0("Seurat/plots/dotplots/", list.name, ".pdf"),
         width = W[list.name], height = H[list.name],
  DotPlot(sro, assay = "RNA", features = rownames(gene.list), group.by = "Clusters",
               col.min = -1, col.max = 1, scale = 15) +
    labs(x = "", y = "Cluster", size = "Percent\nExpressed") +
    scale_x_discrete(breaks = rownames(gene.list), labels = gene.list$Symbol) +
    scale_size_continuous(limits = c(0, 100), range = c(0.1, 8)) +
    scale_color_gradient2(low = "#f9eee8", high = "#5a1419", mid = "#e46c50", midpoint = 0) +
    guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3, title = "Average\nExpression")) +
    theme(axis.text.x.bottom = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
          axis.title.x = element_blank(),
          legend.position = "right", legend.justification = "center",
          legend.key.width = unit(0.5, "inch"), legend.box = "vertical",
          legend.title = element_text(vjust = 1, hjust = 1, size = 5),
          legend.text = element_text(size = 3), plot.margin = unit(rep(1, 4), "mm"))
  )
})

imp.expr <- readRDS("Seurat/results/imputed-expr.rds")
umap <- read.csv("Seurat/results/UMAP_C1-16.csv", row.names = 1, header = T)
id <- rownames(subset(gene.info, Symbol == "Irf8"))
eb <- element_blank()
df <- cbind(umap, expr = imp.expr[id, rownames(umap)])
ggsave(filename = "~/Desktop/Irf8-expr-overlay.pdf", width = 7, height = 5,
ggplot(df[order(df$expr, decreasing = F), ]) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = expr), size = 1, alpha = 0.8) +
  labs(title = gene.info[id, "Symbol"], color = "imputed\nlog-normalized\nexpression") +
  scale_color_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    limits = c(0, max(df$expr))
  ) + theme_classic() +
  theme(axis.ticks = eb, axis.text = eb, axis.title = eb,
        title = element_text(face = "italic"))
)

mid <- subset(cisbp, TF_Name == "Irf8")$Motif_ID
sample <- "3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19#"
df <- cbind(df, t(dev@assays@data$z[mid, paste0(sample, rownames(df))]))
m.ord <- order(sapply(mid, function(id){ cor(df$expr, df[, id]) }), decreasing = T)
plot_grid(
  plotlist = lapply(mid[m.ord], function(id){
    # ci <- order(df[, id], decreasing = F)
    ci <- sample(1:nrow(df), size = nrow(df), replace = F)
    ggplot(df[ci, ]) +
      geom_point(aes(x = UMAP_1, y = UMAP_2, color = df[ci, id]), size = 1, alpha = 0.8) +
      labs(title = id, color = "chromVAR\ndeviation\nz-score") +
      scale_color_gradientn(
        colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))
      ) + theme_classic() +
      theme(axis.ticks = eb, axis.text = eb, axis.title = eb,
            title = element_text(face = "italic"))
  }), ncol = 4
) %>% ggsave(filename = "~/Desktop/Irf8-acc-overlay.pdf", width = 7 * 4, height = 5 * 4)

rna.annot <- c(names(rna.annot)[1:8], rep("LTi", 8))
cl <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1) %>% subset(Clusters %in% 1:16)
cl$Clusters <- ifelse(cl$Clusters %in% 9:16, "9-16", as.character(cl$Clusters)) %>% factor(levels = c(as.character(1:8), "9-16"))
cl$Clusters.annot <- factor(gsub("LTi Variation \\d+$", "LTi", cl$Clusters.annot), levels = rna.annot[1:9])
cl <- cl[order(cl$Clusters, decreasing = T), ]

atac.UMAP <- read.csv("ArchR/ArchROutput_SeuratC1-16/Embeddings/UMAP.csv", header = 1, row.names = 1, col.names = c("cell", "UMAP_1", "UMAP_2"))
atac.md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", header = 1, row.names = 1)
atac.md$Clusters <- as.integer(gsub("C", "", atac.md$Clusters))
atac.md$Clusters <- paste0("C", atac.md$Clusters) %>%
  ifelse(test = atac.md$Clusters == 7, yes = "C6") %>%
  ifelse(test = atac.md$Clusters %in% c(6, 8:13), yes = "C7-13") %>%
  factor(levels = paste0("C", c(1:6, "7-13")))
rna.md <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1)[gsub("^.+#", "", rownames(atac.md)), ]
rna.md$Clusters <- cl[rownames(rna.md), "Clusters"]
rna.md$Clusters.annot <- cl[rownames(rna.md), "Clusters.annot"]

atac.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(atac.pal) <- paste0("C", c(1:6, "7-13"))

umap <- read.csv("Seurat/results/UMAP_C1-16.csv", header = 1, row.names = 1)
cl <- cl[rownames(umap), ]
ci <- order(cl$Clusters, decreasing = T)
plot.clusters(clusters = cl$Clusters.annot[ci], clusters.label = cl$Clusters[ci],
              col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])), vis = umap[ci, ]) %>%
  ggsave(filename = "Seurat/plots/clusters_C1-16.pdf", width = 15, height = 10)
pdf(file = "ArchR/ArchROutput_SeuratC1-16/Plots/UMAP.pdf", width = 15, height = 10)
plot.clusters(clusters = atac.md$Clusters, clusters.label = atac.md$Clusters,
              col = atac.pal, vis = atac.UMAP, pref.C = F)
plot.clusters(clusters = rna.md$Clusters.annot, clusters.label = rna.md$Clusters,
              col = c(rna.pal[1:8], LTi = as.character(rna.pal[11])), vis = atac.UMAP) +
  labs(color = "Clusters.scRNAseq")
dev.off()

eb <- element_blank()
rna.annot <- c(names(rna.annot)[1:8], rep("LTi", 8))
df <- read.csv("Seurat/results/C1-16_UMAP.csv", header = 1, row.names = 1)
df <- cbind(df, seurat.obj@meta.data[rownames(df), c("TRA.score", "Clusters", "Clusters.annot")])
df$Clusters <- ifelse(df$Clusters %in% 9:16, "9-16", as.character(df$Clusters)) %>%
  factor(levels = c(as.character(1:8), "9-16"))
df$Clusters.annot <- factor(gsub("LTi Variation \\d+$", "LTi", df$Clusters.annot), levels = rna.annot[1:9])
pdf("Seurat/plots/TRA-score.pdf", width = 15, height = 12)
ggplot(df[order(!is.na(df$TRA.score), df$TRA.score, decreasing = F), ]) +
  labs(title = "TRA enrichment (4587 genes)", color = "z-score") +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = TRA.score), size = 2, alpha = 0.8) +
  scale_color_distiller(palette = "Spectral", breaks = seq(-2, 8, 2)) +
  theme_classic() + theme(axis.ticks = eb, axis.text = eb, axis.title = eb)
ggplot(subset(df, Clusters %in% c(as.character(c(2:5, 8)), "9-16")),
       aes(x = Clusters, y = TRA.score, fill = Clusters.annot)) + geom_boxplot() +
  scale_fill_manual(values = c(rna.pal[1:8], LTi = as.character(rna.pal[11]))) +
  theme_bw() + labs(title = "TRA enrichment (4587 genes)", y = "z-score", fill = "")
dev.off()

seurat.obj <- readRDS("Seurat/results/SRO.rds")
load("palette.RData")
rna.annot <- c(names(rna.annot)[1:8], rep("LTi", 8))
imp.expr <- readRDS("Seurat/results/imputed-expr.rds")
ci <- seurat.obj$Clusters %in% 1:16
ca <- columnAnnotation(
  Cluster = ifelse(seurat.obj$Clusters %in% 9:16, "9-16", as.character(seurat.obj$Clusters))[ci] %>%
    factor(levels = c(as.character(1:8), "9-16")),
  col = list(Cluster = c(rna.pal2[1:8], `9-16` = as.character(rna.pal2[11])))
)
expr.col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                            c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                              "#f09173", "#d46052", "#b32339", "#690822", "#690822"))

genes <- c("Tnfrsf11a", "Cd40", "Tnfrsf1b", "Tnfrsf1a", "Traf1", "Traf2", "Traf3", "Traf6",
           "Map3k1", "Map3k14", "Ikbkb", "Ikbkg", "Rel", "Rela", "Relb", "Birc2", "Birc3",
           "Tnfaip3", "Nfkb1", "Nfkb2", "Nfkbia", "Nfkbib", "Nfkbiz", "Bmp2k", "Tnfrsf11b",
           "Ltb", "Lta", "Ltbr")
genes <- subset(seurat.obj@assays$RNA@meta.features, Symbol %in% genes)
pdf("Seurat/plots/NFKB-heatmap.pdf", width = 20, height = 10)
Heatmap(t(scale(t(imp.expr[rownames(genes), ci]))), name = "scaled\nimputed\nexpression", col = expr.col.ramp,
        top_annotation = ca, row_labels = genes$Symbol,
        row_names_side = "right", row_names_gp = gpar(fontface = "italic"),
        column_split = factor(gsub("LTi Variation \\d+$", "LTi", seurat.obj$Clusters.annot[ci]), levels = unique(rna.annot)),
        show_column_names = F, cluster_rows = T, cluster_columns = F, column_title_rot = 90) %>%
  draw(padding = unit(c(0.5, 0.5, 1, 0.5), "inch"))
dev.off()

# genes <- subset(seurat.obj@assays$RNA@meta.features, grepl("^KRT", toupper(Symbol)))
# pdf("Seurat/plots/KRT-genes-heatmap.pdf", width = 20, height = 10)
gene.set <- readxl::read_xlsx("../mTEC-multiome/GO_term_summary_20220615_213137.xlsx", sheet = 1, col_names = T)
genes <- subset(seurat.obj@assays$RNA@meta.features, toupper(Symbol) %in% toupper(gene.set$Symbol))
pdf("Seurat/plots/neurotransmitter-genes-heatmap-2.pdf", width = 20, height = 15)
Heatmap(t(scale(t(imp.expr[rownames(genes), ci]))), name = "scaled\nimputed\nexpression", col = expr.col.ramp,
        top_annotation = ca, row_labels = genes$Symbol,
        row_names_side = "right", row_names_gp = gpar(fontface = "italic"),
        column_split = factor(gsub("LTi Variation \\d+$", "LTi", seurat.obj$Clusters.annot[ci]), levels = unique(rna.annot)),
        show_column_names = F, cluster_rows = T, cluster_columns = F, column_title_rot = 90) %>%
  draw(padding = unit(c(0.5, 0.5, 1, 0.5), "inch"))
dev.off()

pdf("~/Desktop/ATAC-vs-RNA.pdf", width = 7, height = 5)
atac.md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", row.names = 1, header = T)
# atac.md$Clusters <- as.integer(gsub("C", "", atac.md$Clusters))
atac.md$Clusters[atac.md$Clusters == "C6"] <- "X"
atac.md$Clusters[atac.md$Clusters == "C7"] <- "C6"
atac.md$Clusters[atac.md$Clusters == "X"] <- "C7"
atac.md$Clusters <- factor(atac.md$Clusters, levels = paste0("C", c(3, 4, 2, 1, 5:13)))
# atac.md$Clusters.scRNAseq <- ifelse(atac.md$Clusters.scRNAseq %in% paste0("C", 9:16),
#                                     "C9-16", atac.md$Clusters.scRNAseq)
atac.md$Clusters.scRNAseq <- factor(atac.md$Clusters.scRNAseq, levels = paste0("C", 1:16))
freq <- as.matrix(table(atac.md[, c("Clusters", "Clusters.scRNAseq")]))
Heatmap(t(apply(freq, 1, function(x){ x / sum(x) })), name = "number\nof cells",
        col = colorRamp2(breaks = c(0, 0.5, 1), colors = c("#dedede", "#a11515", "#780000")),
        cluster_columns = F, cluster_rows = F,
        row_title = "scATAC-seq Cluster", column_title = "scRNA-seq Cluster")
dev.off()

temp.pal <- rna.pal2[c(5, 4, 2, 3)]
names(temp.pal) <- 1:4
fdl <- read.csv("ArchR/ArchROutput_SeuratC1-16/Embeddings/TC-FDL.csv", header = T, row.names = 1)
pt <- read.csv("palantir/results/pseudotime.csv", header = T, row.names = 1)
md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", row.names = 1, header = T)[rownames(fdl), ]
plot_grid(
  ggplot(cbind(fdl, pseudotime = pt[rownames(fdl), 1])) +
    geom_point(aes(x = x, y = y, color = pseudotime)) +
    scale_color_gradientn(colours = viridis::plasma(20), values = seq(0, 1, length.out = 20)) +
    theme_classic() + labs(x = "FDL_1", y = "FDL_2") +
    theme(axis.ticks = element_blank(), axis.text = element_blank()),
  ggplot(cbind(fdl, Cluster = gsub("C", "", md$Clusters))) +
    geom_point(aes(x = x, y = y, color = Cluster)) +
    theme_classic() + labs(x = "FDL_1", y = "FDL_2") +
    scale_color_manual(values = temp.pal) +
    theme(axis.ticks = element_blank(), axis.text = element_blank()),
  ncol = 1
) %>% ggsave(filename = "palantir/pseudotime-on-FDL.pdf", width = 7, height = 15)

atac.md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", header = 1, row.names = 1)
atac.md$Clusters <- atac.md$Clusters %>%
  ifelse(test = atac.md$Clusters == "C7", yes = "C6") %>%
  ifelse(test = atac.md$Clusters %in% paste0("C", c(6, 8:13)), yes = "C7-13") %>%
  factor(levels = c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13"))
rownames(atac.md) <- gsub("^.+#", "", rownames(atac.md))
group.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(group.pal) <- paste0("C", c(1:6, "7-13"))
sro <- subset(seurat.obj, Clusters %in% 1:16)
Itgb8.ID <- subset(sro@assays$RNA@meta.features, Symbol == "Itgb8") %>% rownames()
eb <- element_blank()
ggsave(filename = "~/Desktop/Itgb8-violin.pdf", width = 7, height = 9,
data.frame(Itgb8 = sro@assays$RNA@data[Itgb8.ID, ], Cluster = atac.md[colnames(sro), "Clusters"]) %>%
  ggplot(mapping = aes(x = Itgb8, y = 0)) +
  geom_violin(aes(color = Cluster, fill = Cluster), alpha = 0.9, show.legend = F) +
  facet_wrap(~Cluster, ncol = 1, strip.position = "left") +
  scale_color_manual(values = group.pal) + scale_fill_manual(values = group.pal) +
  theme_classic() + theme(axis.title = eb, axis.text.y = eb, axis.ticks.y = eb,
                          panel.spacing = unit(0, "mm"),
                          panel.border = element_rect(color = "black", fill = NA, size = 0.1))
)
ggsave(filename = "~/Desktop/names.pdf", width = 3, height = 9,
  ggplot(data.frame(names = c("TC I", "TC II", "TC III", "TC IV", "ILC3p", "NCR+", "LTi")),
         aes(x = 0, y = 1:7)) +
    geom_text(aes(label = names, color = c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13")), show.legend = F) +
    scale_color_manual(values = group.pal) + scale_y_reverse() + theme_classic()
)

load("palette.RData")
sro <- readRDS("Seurat/results/SRO.rds")
imp.expr <- readRDS("Seurat/results/imputed-expr.rds")
gene.info <- sro@assays$RNA@meta.features
gene.list <- openxlsx::read.xlsx(
  xlsxFile = "../mTEC-multiome/GO_term_summary_20220615_133242.xlsx",
  sheet = 1, colNames = T, rowNames = F
)
gene.list <- subset(gene.info, toupper(Symbol) %in% toupper(gene.list$Symbol))
cells <- sro$Clusters %in% 1:16
ca <- columnAnnotation(
  Celltype = sro@meta.data[cells, "Clusters.annot"] %>%
    gsub(pattern = "^LTi.+", replacement = "LTi") %>%
    factor(levels = c(names(rna.pal[1:8]), "LTi")),
  col = list(Celltype = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
)
col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                       c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                         "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
pdf(width = 20, height = 15, file = "Seurat/plots/neurotransmitters-receptors.pdf")
Heatmap(
  matrix = t(scale(t(imp.expr[rownames(gene.list), cells]))),
  name = "scaled\nimputed\nexpression", col = col.ramp,
  top_annotation = ca, show_column_names = F,
  column_split = sro@meta.data[cells, "Clusters"],
  cluster_rows = T, row_names_side = "right", row_labels = gene.list$Symbol,
  cluster_columns = T, cluster_column_slices = F,
  column_title_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontface = "italic", fontsize = 5), use_raster = T
)
dev.off()

md <- read.csv("Seurat/results/meta-data.csv", row.names = 1, header = T)
atac.md <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", row.names = 1, header = T)
atac.md$Clusters <- ifelse(test = atac.md$Clusters == "C7", yes = "6", no = atac.md$Clusters) %>%
  ifelse(test = atac.md$Clusters == "C6", yes = "7") %>%
  gsub(pattern = "C", replacement = "") %>% as.integer()
atac.md$Clusters.annot <- c("TC 4", "TC 3", "TC 1", "TC 2", "ILC3p", "NCR+ ILC3", rep("LTi", 7))[atac.md$Clusters]
pdf(file = "Seurat/plots/cluster-breakdown.pdf", width = 15, height = 10)
ggplot(md, aes(x = factor(Clusters), y = ..count..)) +
  geom_bar(aes(fill = Clusters.annot)) + scale_fill_manual(values = rna.pal) +
  geom_label(aes(label = paste0(round(..count.. / sum(..count..), 3) * 100, "%")), stat = "count") +
  geom_text(aes(label = gsub("LTi.+", "LTi", Clusters.annot)), stat = "count", angle = 90, hjust = 0, nudge_y = 50) +
  scale_y_continuous(name = "number of cells", n.breaks = 10, limits = c(0, 2000)) +
  labs(x = "RNA cluster") + theme_bw() + theme(legend.position = "none")
ggplot(atac.md, aes(x = factor(Clusters), y = ..count..)) +
  geom_bar(aes(fill = Clusters.annot)) + scale_fill_manual(values = c(rna.pal[1:8], LTi = as.character(rna.pal[11]))) +
  geom_label(aes(label = paste0(round(..count.. / sum(..count..), 3) * 100, "%")), stat = "count") +
  geom_text(aes(label = Clusters.annot), stat = "count", angle = 90, hjust = 0, nudge_y = 50) +
  scale_y_continuous(name = "number of cells", n.breaks = 10) +
  labs(x = "ATAC cluster") + theme_bw() + theme(legend.position = "none")
dev.off()

pal.res <- cbind(
  read.csv(
    file = "palantir/results/pseudotime.csv",
    header = T, row.names = 1, col.names = c("id", "pseudotime")
  ),
  read.csv(
    file = "palantir/results/entropy.csv",
    header = T, row.names = 1, col.names = c("id", "entropy")
  ),
  read.csv(
    file = "palantir/results/branch_probs.csv",
    header = T, row.names = 1, check.names = F
  )
)
vis <- read.csv(
  file = "Seurat/results/FDL.csv",
  header = T, row.names = 1, col.names = c("cell", "UMAP_1", "UMAP_2")
)
md <- read.csv("Seurat/results/meta-data.csv", header = T, row.names = 1)
md$Clusters <- factor(md$Clusters)
md$Clusters.annot <- factor(md$Clusters.annot, levels = names(rna.annot))
md <- md[order(as.integer(md$Clusters), decreasing = T), ]
vis <- vis[order(as.integer(md[rownames(vis), "Clusters"]), decreasing = T), ]
pal.res <- pal.res[rownames(vis), ]
start <- "GTATGTTCATCCAGGT-1"
ggsave(
  filename = "palantir/pseudotime-1.pdf", width = 12 * 2, height = 10 * 2,
  plot =
    plot_grid(
      plot.continuous.value(
        vis = vis, val = pal.res$pseudotime, val.name = "Palantir\npseudotime",
        scale.color = scale_color_gradient2(low = "#fcff9c", high = "#6e0000",
                                            mid = "#ff4900", midpoint = 0.5)) +
        geom_point(aes(x = vis[start, 1], y = vis[start, 2]), size = 3, shape = 8) +
        geom_text_repel(aes(x = vis[start, 1], y = vis[start, 2], label = "start"),
                        min.segment.length = 0, segment.linetype = 2, size = 10,
                        nudge_x = 0, nudge_y = 1000),
      plot.continuous.value(
        vis = vis, val = pal.res$entropy, val.name = "Palantir\nentropy",
        scale.color = scale_color_gradient(low = "#9df5ef", high = "#011463")) +
        geom_point(aes(x = vis[start, 1], y = vis[start, 2]), size = 3, shape = 8) +
        geom_text_repel(aes(x = vis[start, 1], y = vis[start, 2], label = "start"),
                        min.segment.length = 0, segment.linetype = 2, size = 10,
                        nudge_x = 0, nudge_y = 1000),
      plot.clusters(
        clusters = gsub("LTi.+", "LTi", md$Clusters.annot) %>% factor(levels = c(names(rna.pal[1:8]), "LTi")),
        vis = vis[rownames(md), ], labels = F, col = c(rna.pal[1:8], LTi = as.character(rna.pal[11]))
      ) + theme(axis.line = element_line())
    )
)
plot_grid(
  plotlist = lapply(3:ncol(pal.res), function(i){
    end.point <- colnames(pal.res)[i]
    nx <- c(1, 0, 3)[i - 2] * 500
    ny <- c(2, 8, 4)[i - 2] * 500
    plot.continuous.value(
      vis = vis, val = pal.res[, end.point], val.name = "branch\nprobability",
      scale.color = scale_color_distiller(palette = "Spectral")) +
      geom_point(aes(x = vis[start, 1], y = vis[start, 2]), size = 3, shape = 8) +
      geom_point(aes(x = vis[end.point, 1], y = vis[end.point, 2]), size = 3, shape = 8) +
      geom_text_repel(aes(x = vis[start, 1], y = vis[start, 2], label = "start"),
                      min.segment.length = 0, segment.linetype = 2, size = 5,
                      nudge_x = 0, nudge_y = 1000) +
      geom_text_repel(aes(x = vis[end.point, 1], y = vis[end.point, 2], label = "terminal state"),
                      min.segment.length = 0, segment.linetype = 2, size = 5,
                      nudge_x = nx, nudge_y = ny)
  }), ncol = 3
) %>% ggsave(filename = "palantir/pseudotime-2.pdf", width = 12 * 3, height = 10)


