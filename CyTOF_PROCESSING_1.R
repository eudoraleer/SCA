######### Individual CSV ############
library("dplyr")
library("limma")
library("corrplot")
library("gplots")
library("plot3D")
library("ggplot2")
library("edgeR")
library("org.Hs.eg.db")
library("ggthemes")
library("patchwork")
library("Seurat")
library("Matrix")
library("ggridges")
library("ggrepel")

source('~/Dropbox/KI/Studies/LXXLZ/Website/Web_Scripts/SCA_Analysis_Functions_V1.0.0.R')
# source('/proj/store_eudoraleer/PROJECTS/scRNASEQ/SCRIPTS/SCA_Analysis_Functions_V1.0.0.R')
color_conditions <- color_ini()

# output_dir <- paste(system("pwd", intern = T),"/", sep = "")
output_dir <- "~/Dropbox/KI/Studies/LXXLZ/Website/CyTOF/CyTOF_OUTPUT/P131_Mendeley_Data_PBMC_CyTOF/"
input_dir <- "~/Desktop/Desktop/CyTOF_ALL/P131_Mendeley_Data_PBMC_CyTOF/"
log_dir <- "~/Dropbox/KI/Studies/LXXLZ/Website/CyTOF/CyTOF_OUTPUT/P131_Mendeley_Data_PBMC_CyTOF/"

output_dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/CyTOF/OUTPUT/"
log_dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/CyTOF/LOG/"
setwd(output_dir)
args <- commandArgs(TRUE)
file <- as.character(args[1])

cname <- gsub("\\.csv","",file, ignore.case = T)
project_name <- gsub(".*\\/(.*)","\\1",system("pwd", intern = T))

x <- read.csv(file)
colnames(x)

# Input to seurat
meta <- data.frame(CELL_ID = paste("CELL_",1:nrow(x),sep = ""),
                   SAMPLE_ID = x$SAMPLE_ID,
                   TISSUE = x$TISSUE,
                   SUBTISSUE = x$SUBTISSUE)

x <- x[,grep("cell_id|sample_id|TISSUE|SUBTISSUE|gate_source", colnames(x), ignore.case = T, invert = T)]
x <- t(x)
dim(x)
colnames(x) <- meta$CELL_ID
# markers <- toupper(row.names(x))
# markers <- gsub("\\.|\\-|\\s+","_",markers)

x <- CreateSeuratObject(x)
head(x@meta.data)

x@meta.data <- cbind(x@meta.data, meta)
x$PROJECT <- project_name
x$orig.ident <- x$SAMPLE_ID
x@assays$RNA@data <- asinh(x@assays$RNA@counts/5)
x <- ScaleData(x)

somePNGPath = paste(output_dir,"0CyTOF_RNA_INFO_RAW_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=5000, height=3000, units = "px", res = 300)
print(VlnPlot(x,features="nCount_RNA", group.by = "SAMPLE_ID", pt.size = 0,
              cols = gen_colors(color_conditions$colorful, length(unique(x$SAMPLE_ID)))))
dev.off()

current <- data.frame(CELL_ID = x$CELL_ID, SAMPLE_ID = x$SAMPLE_ID, nCount_RNA = x$nCount_RNA)
# current$isOUTLIER <- "NO"
plotx <- split(current,current$SAMPLE_ID)
plotx <- lapply(plotx, function(x){
  outliers <- NULL
  outliers <- boxplot(x$nCount_RNA, plot=FALSE)$out
  x <- x[which(!x$nCount_RNA %in% outliers),]
  return(x)
})

plotx <- do.call(rbind.data.frame, plotx)

head(row.names(x@meta.data))
# x$isOUTLIER <- "YES"
# x@meta.data[which(x$CELL_ID %in% plotx$CELL_ID),"isOUTLIER"] <- "NO"
# x <- subset(x, subset = isOUTLIER == "NO")
x <- subset(x, cells = plotx$CELL_ID)

somePNGPath = paste(output_dir,"0CyTOF_RNA_INFO_FILTERED_OUTLIERS_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=5000, height=3000, units = "px", res = 300)
print(VlnPlot(x,features="nCount_RNA", group.by = "SAMPLE_ID", pt.size = 0,
              cols = gen_colors(color_conditions$colorful, length(unique(x$SAMPLE_ID)))))
dev.off()

x <- subset(x, subset = nCount_RNA < 3000)
m <- data.frame(FILE = file,NUM_CELLS_POST_FILTER = ncol(x))
write.table(m,paste(log_dir,"/POST_FILTER_CELL_NUM_3K",cname,".txt", sep = ""), quote=F, row.names=F)

somePNGPath = paste(output_dir,"0CyTOF_FILTERED_RNA_3K_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=5000, height=3000, units = "px", res = 300)
print(VlnPlot(x,features="nCount_RNA", group.by = "SAMPLE_ID", pt.size = 0,
              cols = gen_colors(color_conditions$colorful, length(unique(x$SAMPLE_ID)))))
dev.off()

# saveRDS(x, gsub("\\.csv",".RDS",file, ignore.case = T))

plotx <- data.frame(t(x@assays$RNA@data))
plotx$SAMPLE_ID <- paste("S",x$SAMPLE_ID, sep = "")
plotx$SUBTISSUE <- x$SUBTISSUE
ggdf <- reshape2::melt(plotx)
colnames(ggdf) <- c("SAMPLE_ID","SUBTISSUE","MARKER","ASINH")
head(ggdf)

somePNGPath = paste(output_dir,"0CyTOF_RIDGE_MARKERS_FILTERED_3K_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=8000, height=6000, units = "px", res = 300)
p <- ggplot(ggdf,
            aes(x = ASINH, y = SAMPLE_ID, fill = SAMPLE_ID)) +
  ggtitle(paste(cname, ": NO FILTERING"))+
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  facet_wrap(~MARKER, nrow = 6, scales = "free") +
  theme_hc() + ylab("DENSITY") + xlab("ASINH EXPRESSION")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        axis.text.x = element_text(angle = 0, size = 25),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 35, face=2, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 35, face=2, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        # legend.key.size = unit(1, "cm"),
        legend.position = "none",
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(size = 40, face=2, hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_fill_manual(values = gen_colors(color_conditions$tableau20, length(unique(ggdf$MARKER))))
print(p)
dev.off()

# Check if there is batch effect
x <- FindVariableFeatures(x)
x <- RunPCA(x, features = VariableFeatures(x))
x <- RunUMAP(x, reduction = "pca", dims = 1:30)
x <- RunTSNE(x, reduction = "pca", dims = 1:30, check_duplicates = FALSE)

saveRDS(x,paste(cname,"_POST_tSNE.RDS", sep = ""))

plotx <- gen10x_plotx(x)
plotx$SAMPLE_ID <- paste("S",x$SAMPLE_ID, sep = "")

p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", title = paste(cname,"\n\nUMAP - CyTOF SAMPLES"),isfacet = T,
                        col=color_conditions$tenx, color_by = "SAMPLE_ID", group_by = "SAMPLE_ID",
                        xlabel = "UMAP1", ylabel = "UMAP2")

somePNGPath <- paste(output_dir,cname,"_UMAP_LABELED_BEFORE_INTEGRATION.png", sep = "")
png(somePNGPath, width = 4000, height =3200, units = "px", res = 200)
print(p1)
dev.off()

p2 <- own_facet_scatter(plotx, "tSNE_1", "tSNE_2", title = paste(cname,"\n\ntSNE - CyTOF SAMPLES"),isfacet = T,
                        col=color_conditions$tenx, color_by = "SAMPLE_ID", group_by = "SAMPLE_ID",
                        xlabel = "tSNE1", ylabel = "tSNE2")

somePNGPath <- paste(output_dir,cname,"_TSNE_LABELED_BEFORE_INTEGRATION.png", sep = "")
png(somePNGPath, width = 4000, height =3200, units = "px", res = 200)
print(p2)
dev.off()

x$SAMPLE_ID <- paste("S",x$SAMPLE_ID, sep = "")
data <- SplitObject(x, split.by = "SAMPLE_ID")
data_features <- SelectIntegrationFeatures(object.list = data)
data <- lapply(data, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- FindVariableFeatures(x)
  x <- ScaleData(x, verbose=F, features = data_features)
  x <- RunPCA(x, npcs = 30, verbose = FALSE, features = data_features)
})
data_anchors <- FindIntegrationAnchors(object.list=data, reduction = "rpca", anchor.features = data_features)
data_integrated <- IntegrateData(data_anchors, dims = 1:30)
DefaultAssay(data_integrated) <- "integrated"
data_integrated <- ScaleData(data_integrated, verbose = FALSE)
data_integrated <- RunPCA(data_integrated, verbose = FALSE)
data_integrated <- RunTSNE(data_integrated, reduction = "pca", dims = 1:30)
data_integrated <- RunUMAP(data_integrated, reduction = "pca", dims = 1:30)
data_integrated <- FindNeighbors(data_integrated, reduction = "pca", dims = 1:30)
data_integrated <- FindClusters(data_integrated, resolution = 1)
saveRDS(data_integrated,paste(cname,"_POST_INTEGRATION_POST_FINDCLUSTERS.RDS", sep = ""))

DefaultAssay(x) <- "RNA"
somePNGPath <- paste(output_dir,"Preintegration_PC12_",cname, "_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
DimPlot(object = x, reduction = "pca", pt.size = .1, group.by = "SAMPLE_ID")
dev.off()

# somePNGPath <- paste(output_dir,"Preintegration_PC12_BY_BATCH_",project_name,".png", sep = "")
# png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
# DimPlot(object = x, reduction = "pca", pt.size = .1, group.by = "BATCH")
# dev.off()

somePNGPath <- paste(output_dir,"Preintegration_PC1_BY_SAMPLES_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
VlnPlot(object = x, features = "PC_1", group.by = "SAMPLE_ID", pt.size = .1)
dev.off()

# somePNGPath <- paste(output_dir,"Preintegration_PC1_BY_BATCH_",project_name,".png", sep = "")
# png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
# VlnPlot(object = x, features = "PC_1", group.by = "BATCH", pt.size = .1)
# dev.off()

DefaultAssay(data_integrated) <- "RNA"
data_integrated <- SCTransform(data_integrated)
DefaultAssay(data_integrated) <- "SCT"
data_integrated <- RunHarmony(data_integrated, "SAMPLE_ID", plot_convergence = TRUE, assay.use = "SCT")

somePNGPath <- paste(output_dir,"Postintegration_HARMONY12_BY_SAMPLE_ID_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
DimPlot(object = data_integrated, reduction = "harmony", pt.size = .1, group.by = "SAMPLE_ID")
dev.off()

DefaultAssay(data_integrated) <- "integrated"
somePNGPath <- paste(output_dir,"Postintegration_PC12_BY_SAMPLE_ID_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
DimPlot(object = data_integrated, reduction = "pca", pt.size = .1, group.by = "SAMPLE_ID")
dev.off()

# somePNGPath <- paste(output_dir,"Postintegration_HARMONY1_BY_BATCH_",project_name,".png", sep = "")
# png(somePNGPath, width = 1000, height = 800, units = "px", res = 100)
# VlnPlot(object = x, features = "harmony_1", group.by = "BATCH", pt.size = .1)
# dev.off()

DefaultAssay(data_integrated) <- "RNA"
somePNGPath <- paste(output_dir,"Postintegration_HARMONY_1_TO_15_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 2000, units = "px", res = 100)
DimHeatmap(data_integrated,reduction = "harmony", dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
dev.off()

DefaultAssay(data_integrated) <- "integrated"
somePNGPath <- paste(output_dir,"Postintegration_PCA_1_TO_15_",project_name,".png", sep = "")
png(somePNGPath, width = 1000, height = 2000, units = "px", res = 100)
DimHeatmap(data_integrated,reduction = "pca", dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE, assays = "integrated")
dev.off()

DefaultAssay(data_integrated) <- "RNA"
plotx <- data.frame(t(data_integrated@assays$RNA@data))
plotx$SAMPLE_ID <- data_integrated$SAMPLE_ID
plotx$SUBTISSUE <- data_integrated$SUBTISSUE
ggdf <- reshape2::melt(plotx)
colnames(ggdf) <- c("SAMPLE_ID","SUBTISSUE","MARKER","ASINH")
head(ggdf)

somePNGPath = paste(output_dir,"0CyTOF_RIDGE_MARKERS_RAW_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=8000, height=6000, units = "px", res = 300)
p <- ggplot(ggdf,
            aes(x = ASINH, y = SAMPLE_ID, fill = SAMPLE_ID)) +
  ggtitle(paste(cname, ": NO FILTERING"))+
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  facet_wrap(~MARKER, nrow = 6, scales = "free") +
  theme_hc() + ylab("DENSITY") + xlab("ASINH EXPRESSION")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        axis.text.x = element_text(angle = 0, size = 25),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 35, face=2, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 35, face=2, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        # legend.key.size = unit(1, "cm"),
        legend.position = "none",
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(size = 40, face=2, hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_fill_manual(values = gen_colors(color_conditions$tableau20, length(unique(ggdf$MARKER))))
print(p)
dev.off()

# remove outliers
markers <- unique(ggdf$MARKER)
samples <- unique(ggdf$SAMPLE_ID)
filtered_ggdf <- NULL
for(i in 1:length(markers)){
  for(j in 1:length(samples)){
    print(paste("Current sample: ", samples[j], " & ", markers[i], sep = ""))
    current <- ggdf[which(ggdf$SAMPLE_ID == samples[j] & ggdf$MARKER == markers[i]),]
    outliers <- NULL
    outliers <- boxplot(current$ASINH, plot=FALSE)$out
    if(length(outliers) > 0){
      current <- current[which(!current$ASINH %in% outliers),]
    }
    filtered_ggdf <- rbind(filtered_ggdf, current)
    # current <- subset(current, current$ASINH > (Q[1] - 1.5*iqr) & warpbreaks$breaks < (Q[2]+1.5*iqr))
  }
}

saveRDS(filtered_ggdf,paste("FILTERED_OUTLIERS_ASINH_EXPR_",cname,".RDS",sep = ""))
# ggdf <- filtered_ggdf
# max(filtered_ggdf[filtered_ggdf$MARKER == "CD116","ASINH"])

somePNGPath = paste(output_dir,"0CyTOF_RIDGE_MARKERS_FILTERED_OUTLIERS_",project_name,"_",cname,".png", sep = "")
png(somePNGPath, width=8000, height=6000, units = "px", res = 300)
p <- ggplot(filtered_ggdf, #ggdf[sample(1:nrow(ggdf), size = 10000, replace = F),],  #ggdf[ggdf$MARKER == "CD4",],
            aes(x = ASINH, y = SAMPLE_ID, fill = SAMPLE_ID)) +
  ggtitle(cname)+
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  facet_wrap(~MARKER, nrow = 6, scales = "free") +
  theme_hc() + ylab("DENSITY") + xlab("ASINH EXPRESSION")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 35, face=2, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 35, face=2, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        # legend.key.size = unit(1, "cm"),
        legend.position = "none",
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(size = 40, face=2, hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(ggdf$MARKER))))
print(p)
dev.off()

saveRDS(x, gsub("\\.csv",".RDS",file, ignore.case = T))

samples <- unique(plotx$SAMPLE_ID)
markers <- colnames(plotx)[grep("file|peak|sne|population|SAMPLE_ID|tissue",colnames(plotx), invert = T, ignore.case = T)]

plot_median <- NULL
expr <- NULL

for(i in 1:length(samples)){
  current <- plotx[which(plotx$SAMPLE_ID == samples[i]),]
  current <- current[,grep("file|peak|sne|population|SAMPLE_ID|tissue",colnames(current), invert = T, ignore.case = T)]
  colnames(current)
  for (j in 1:ncol(current)){
    current[,j] <- as.numeric(as.character(current[,j]))
    expr <- c(expr,median(current[current[,j] > 0,j]))
  }
  plot_median <- rbind(plot_median, expr)
  expr <- NULL
}

plot_median <- data.frame(t(plot_median))
row.names(plot_median) <- samples
colnames(plot_median) <- markers
saveRDS(plot_median,paste(cname,"_MEDIAN_ASINH_EXPRESSION.RDS",sep = ""))

plot_median <- scale(plot_median)
plot_median <- t(scale(t(plot_median)))
somePNGPath = paste(output_dir,cname,"_MEDIAN_ASINH_EXPR.png", sep = "")
png(somePNGPath, width=4000, height=3100, units = "px", res = 300)
heatmap.2(as.matrix(t(plot_median)),margin=c(15, 15), trace="none",key=T, keysize=1,col=jet2.col (n = 100, alpha = 1),srtCol=45, scale="none", Colv = T,
          Rowv = T,main = paste("CyTOF - ", cname, " MEDIAN ASINH EXPR"),
          density.info="none", cexCol=0.8,cexRow=0.8)
dev.off()

save.image(paste(cname,"RData", sep = ""))

mds <- plotMDS(plot_median, plot = FALSE)
pca_out <- prcomp(t(plot_median), center = TRUE, scale. = FALSE)
ggdf <- data.frame(SAMPLE_ID = colnames(plot_median), MDS1 = mds$x, MDS2 = mds$y, PCA1 = pca_out$x[,1], PCA2 = pca_out$x[,2])
ggdf <- plyr::join(ggdf, pheno_data, by = "File")
save.image(paste(cname,"RData", sep = ""))












