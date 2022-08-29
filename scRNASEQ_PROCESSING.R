library("SingleR")
hpca.se <- HumanPrimaryCellAtlasData()
library("Seurat", lib="/proj/eudoraleer/")
Packages <- c("monocle3","dplyr","patchwork", #"scCATCH",
              "reshape2","SeuratData","Signac","reticulate","plotly",
              "viridis","cowplot","plot3D","DDRTree","pheatmap",
              "scater","ggbeeswarm","ggthemes","M3Drop","gplots",
              "gridExtra","ggpubr","plotly","lattice","ReactomePA",
              "HGNChelper","DOSE","enrichplot","harmony","bbknn",
              "ggrepel","akmedoids","grid","multtest","metap",
              "akmedoids","htmlwidgets","monocle","celldex",
              "ggnewscale","clusterProfiler","ggupset","europepmc",
              "dendextend","cicero","org.Hs.eg.db","scRNAseq")

suppressMessages(lapply(Packages, library, character.only = TRUE))
source('/proj/store_eudoraleer/PROJECTS/scRNASEQ/SCRIPTS/SCA_Analysis_Functions_V1.0.0.R')
color_conditions <- color_ini()

mito_threshold <- 30
hs <- org.Hs.eg.db
hgnc.table <- getCurrentHumanMap()
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

files <- read.table("/proj/lxx_storage/TISSUE_OMICS/GEO/LOG/ALL_scRNASEQ_RDS.txt", header = F)
colnames(files) <- "FILE"
files$NORMALIZED <- ""
output_dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/"
sout <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/SAMPLEWISE/"
allout <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/INTEGRATED/"
rdsout <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/RDS/"
project_name <- "HUMAN TISSUE scRNASEQ"
data_current <- NULL

for(j in 1:nrow(files)){
  data_current[[j]] <- readRDS(files[j,"FILE"])
  if(length(grep("\\.",(data_current[[j]]@meta.data$nCount_RNA)))>0){
    files[j,"NORMALIZED"] <- "YES"
  }else{
    files[j,"NORMALIZED"] <- "NO"
  }
  
  sample_name <- gsub(".*\\/(.*)\\.RDS","\\1",files[j,"FILE"])

  data_current[[j]][["Percent_Mito"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
  if(max(data_current[[j]][["Percent_Mito"]]) == 0){
    plot_mito <- FALSE
  }else{
    plot_mito <- TRUE
  }
  
  data_current[[j]] <- add_names(data_current[[j]], data_current[[j]]$SAMPLE_ID, current_ident = "orig.ident")
  
  plotx <- data_current[[j]]@meta.data
  
  somePDFPath = paste(sout,"01.1SCA_PLOT_PREFILTERING_RNA_INFO_VIOLIN_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=8,pointsize=10)
  
  p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 18)
  p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 18)
  if(plot_mito == TRUE){
    p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 18)
    p <- p1+p2+p3
  }else{p <- p1+p2}
  print(p+plot_annotation(title = sample_name, theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
  data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > 200 & nFeature_RNA <= 20000 & Percent_Mito < mito_threshold)
  
  DefaultAssay(data_current[[j]]) <- 'RNA'
  
  if(files[j,"NORMALIZED"] == "NO"){
    data_current[[j]] <- NormalizeData(data_current[[j]], verbose = TRUE)
  }
  data_current[[j]] <- FindVariableFeatures(object = data_current[[j]])
  data_current[[j]] <- suppressWarnings(ScaleData(data_current[[j]], verbose = FALSE))
  data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
  data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
  data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
  

  if(plot_mito == TRUE){
    metrics <-  c("nCount_RNA", "nFeature_RNA", "Percent_Mito")
  }else{
    metrics <-  c("nCount_RNA", "nFeature_RNA")
  }  
  
  somePDFPath = paste(sout,"01.2SCA_PLOT_POST_FILTERING_UMAP_METRICS_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=ifelse(length(metrics) <=2, 6, ceiling(length(metrics))/3*8),pointsize=10) 
  print(FeaturePlot(data_current[[j]], 
                    reduction = "umap", 
                    features = metrics,
                    pt.size = 0.4, 
                    order = TRUE,
                    min.cutoff = 'q10',
                    label = FALSE, cols = c("green","blue"))+
          plot_annotation(title = paste("UMAP Features - ",sample_name, sep = ""),
                          theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))))
  dev.off()
  
  # head(data_current[[j]]@meta.data)
  data_current[[j]] <- suppressWarnings(SCTransform(data_current[[j]], verbose = FALSE))
  plotx <- data_current[[j]]@meta.data
  somePDFPath = paste(sout,"02.0SCA_PLOT_POSTFILTERING_RNA_INFO_VIOLIN_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=8,pointsize=10)
  p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 15)
  p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 15)
  if(plot_mito == TRUE){
    p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 15)
    p <- p1+p2+p3
  }else{
    p <- p1+p2
  }
  print(p+plot_annotation(title = sample_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  
  dev.off()
  
  n <- 10
  var <- c('gmean', 'variance', 'residual_variance')
  data_current[[j]][["SCT"]]@meta.features <- SCTResults(data_current[[j]][["SCT"]], slot = "feature.attributes")[, var]
  data_current[[j]][["SCT"]]@meta.features$variable <- FALSE
  data_current[[j]][["SCT"]]@meta.features[VariableFeatures(data_current[[j]][["SCT"]] ), "variable"] <- TRUE
  colnames(data_current[[j]][["SCT"]]@meta.features) <- paste0("sct.", colnames(data_current[[j]][["SCT"]]@meta.features))
  
  DefaultAssay(data_current[[j]]) <- "RNA"
  
  data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
  
  somePDFPath = paste(sout,"02.3SCA_PLOT_TOP_FEATURES_ASSOC_PC12_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=7,pointsize=10)
  print(VizDimLoadings(data_current[[j]], dims = 1:2, reduction = "pca",
                       col = color_conditions$mark[sample(1:length(color_conditions$mark), size = 1)])+
          plot_annotation(title = sample_name, theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
  
  somePDFPath = paste(sout,"02.6SCA_PLOT_HEATMAP_TOP_PCS_TOP_FEATURES_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=18,pointsize=10)
  print(DimHeatmap(data_current[[j]], dims = 1:15, cells = 500, 
                   # assays = "RNA", 
                   balanced = TRUE, fast = FALSE)+plot_annotation(title = paste(sample_name, sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  dev.off()
  
  plotx <- data.frame(PC = 1:30, SD = data_current[[j]]@reductions$pca@stdev[1:30])
  plotx$Label <- ""
  num_pcs <- ifelse(ceiling(elbow_point(plotx$PC,plotx$SD)$x)<10, 10,
                    ceiling(elbow_point(1:30,data_current[[j]]@reductions$pca@stdev[1:30])$x))

  data_current[[j]] <- FindNeighbors(data_current[[j]], dims = 1:num_pcs)
  data_current[[j]] <- FindClusters(data_current[[j]], resolution = 0.8)
  data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30), check_duplicates = FALSE)
  data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
  
  plotx <- gen10x_plotx(data_current[[j]])
  plotx$CLUSTER <- Idents(data_current[[j]])
  plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
  
  p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER", plot_title = "UMAP",
                     col = NULL, annot = TRUE, legend_position = "right", numeric = TRUE)
  p2 <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CLUSTER", plot_title = "tSNE",
                     col = NULL, annot = TRUE, legend_position = "right", numeric = TRUE)
  
  somePDFPath = paste(sout,"02.8SCA_PLOT_UMAP_AUTOCLUSTER_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=9,pointsize=10)
  print(p1+plot_annotation(title = sample_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  dev.off()
  
  somePDFPath = paste(sout,"02.9SCA_PLOT_TSNE_AUTOCLUSTER_",sample_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=9,pointsize=10)
  print(p2+plot_annotation(title = sample_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  dev.off()
  
  saveRDS(data_current[[j]], paste(rdsout,"/POST_FILTER_",sample_name,".RDS",sep = ""))
}

data_current <- lapply(data_current, function(x){
  # x <- ScaleData(x, verbose=F, features = data_features, vars.to.regress = c("nCount_RNA", "Percent_Mito"))
  # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = data_features)
  # x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
})

integrative_features <- SelectIntegrationFeatures(object.list = data_current)
data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                       reduction = "rpca", anchor.features = integrative_features)
data <- IntegrateData(anchorset = data_anchors)
DefaultAssay(data) <- "integrated"
data <- ScaleData(data, verbose = FALSE)

saveRDS(data, paste(allout,"/INTEGRATED_ALL_TISSUES_scRNASEQ.RDS",sep = ""))

DefaultAssay(data) <- "RNA"
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
data <- RunTSNE(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
data_dim <- data.frame(gen10x_plotx(data), DATA_TYPE = "BEFORE_INTEGRATION", SAMPLE_ID = data$orig.ident)
write.csv(data_dim, paste(allout, "05.7SCA_TABLE_DIM_PARAMETERS_BEFORE_INTEGRATION.csv", sep = ""), quote = F, row.names = T)

integration_method <- "SEURAT"
integration_name <- "SEURAT_INTEGRATED"
DefaultAssay(data) <- "integrated"
reduction_method <- "pca"

data <- RunPCA(data, verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
data <- RunTSNE(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
current <- cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
write.csv(current, paste(allout, "05.8SCA_TABLE_DIM_PARAMETERS_POST_",integration_method, "_INTEGRATION.csv", sep = ""), quote = F, row.names = T)

data_dim <- rbind(data_dim, cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident)))

p <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                         data_dim[data_dim$DATA_TYPE == integration_name,],
                         dim1= "UMAP_1", dim2 = "UMAP_2", group = "SAMPLE_ID",
                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                         maintitle = "HUMAN TISSUES", titlesize = 35, col = color_conditions$tenx)

somePDFPath = paste(allout,"05.9SCA_PLOT_UMAP_BEFORE_POST_",integration_method, "_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=30, height=15,pointsize=10)
print(p)
dev.off()

p <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                         data_dim[data_dim$DATA_TYPE == integration_name,],
                         dim1= "tSNE_1", dim2 = "tSNE_2", group = "SAMPLE_ID",
                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                         maintitle = "HUMAN TISSUES", titlesize = 35, col = color_conditions$tenx)

somePDFPath = paste(allout,"06.0SCA_PLOT_TSNE_BEFORE_POST_",integration_method, "_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=30, height=15,pointsize=10)
print(p)
dev.off()

p <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                         data_dim[data_dim$DATA_TYPE == integration_name,],
                         dim1= "PC_1", dim2 = "PC_2", group = "SAMPLE_ID",
                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                         maintitle = "HUMAN TISSUES", titlesize = 35, col = color_conditions$tenx)

somePDFPath = paste(allout,"06.1SCA_PLOT_PCA_BEFORE_POST_",integration_method, "_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=30, height=15,pointsize=10)
print(p)
dev.off()

data_dim$TISSUE <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"TISSUE"]
data_dim$SUBTISSUE <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"GROUP"]
data_dim$BATCH <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"PROJECT_NAME"]

p <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                         data_dim[data_dim$DATA_TYPE == integration_name,],
                         dim1= "UMAP_1", dim2 = "UMAP_2",group = "TISSUE",
                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                         maintitle = "ALL TISSUES", titlesize = 35, col = color_conditions$warm)

somePDFPath = paste(allout,"06.2SCA_PLOT_BY_GROUP_UMAP_BEFORE_AFTER_",integration_method, "_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=30, height=15,pointsize=10)
print(p)
dev.off()


p <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                         data_dim[data_dim$DATA_TYPE == integration_name,],
                         dim1= "tSNE_1", dim2 = "tSNE_2",group = "TISSUE",
                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                         maintitle = "ALL TISSUES", titlesize = 35, col = color_conditions$warm)

somePDFPath = paste(allout,"06.2SCA_PLOT_BY_GROUP_TSNE_BEFORE_AFTER_",integration_method, "_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=30, height=15,pointsize=10)
print(p)
dev.off()

p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,],
                   x = "UMAP_1", y = "UMAP_2", group = "TISSUE", plot_title = "UMAP",
                   col = NULL, annot = TRUE, legend_position = "right", numeric = TRUE)
p2 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,],
                   x = "tSNE_1", y = "tSNE_2", group = "TISSUE", plot_title = "tSNE",
                   col = NULL, annot = TRUE, legend_position = "right", numeric = TRUE)

somePDFPath = paste(allout,"06.3SCA_PLOT_UMAP_AUTOCLUSTER_POST_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=12, height=9,pointsize=10)
print(p1+plot_annotation(title = "ALL TISSUES - POST INTEGRATION UMAP", theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
dev.off()

somePDFPath = paste(allout,"06.4SCA_PLOT_TSNE_AUTOCLUSTER_POST_INTEGRATION.pdf", sep = "")
pdf(file=somePDFPath, width=12, height=9,pointsize=10)
print(p2+plot_annotation(title = "ALL TISSUES - POST INTEGRATION tSNE", theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
dev.off()

save.image(paste(allout,"/INTEGRATED_ALL_TISSUES_scRNASEQ_POST_06.4.RData",sep = ""))

DefaultAssay(data) <- "integrated"

somePDFPath = paste(allout,"07.1SCA_PLOT_HEATMAP_TOP_PCS_TOP_FEATURES_INTEGRATED.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=18,pointsize=10)
print(DimHeatmap(data, assays = 'integrated', reduction = "pca", dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)+
        plot_annotation(title = paste("All TIssues - TOP 15 PCs", sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
dev.off()

plotx <- data.frame(PC = 1:30, SD = data@reductions$pca@stdev[1:30])
elbow_y <- "PC"
plotx$Label <- ""
num_pcs <- ifelse(ceiling(elbow_point(plotx[,elbow_y],plotx$SD)$x)<10, 10,
                  ceiling(elbow_point(1:30,data@reductions[[reduction_method]]@stdev[1:30])$x))

plotx[plotx[,elbow_y] == num_pcs,"Label"] <- paste("Selected No. of ",elbow_y,"s = ", num_pcs, sep = "")
p <- ggplot(plotx, aes(plotx[,elbow_y], SD, label = Label, color = Label))+
  geom_point(size = 5)+
  theme_classic() +
  geom_point(data=plotx[plotx[,elbow_y] == num_pcs,], aes(x=plotx[plotx[,elbow_y] == num_pcs,elbow_y], y=SD))+
  scale_color_manual(values = color_conditions$tenx)+
  geom_text_repel(box.padding = 10, max.overlaps = Inf, color = color_conditions$tenx[2], size = 10)+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
  xlab(elbow_y)+
  ggtitle("All Tissues")

somePDFPath = paste(zip_out,"07.2SCA_PLOT_PC_VS_SD_",integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=12, height=8,pointsize=10)
print(p)
dev.off()

data <- FindNeighbors(data, reduction = reduction_method, dims = 1:num_pcs)
data <- FindClusters(data, resolution = 0.8)

integration_cluster <- "integrated_snn_res.0.8"

plotx <- gen10x_plotx(data)
plotx$CLUSTER <- Idents(data)
plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
plotx$TISSUE <- data$TISSUE
plotx$SUBTISSUE <- data$GROUP
current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))),decreasing = F)

p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", title = project_name,isfacet = F,
                        col=color_conditions$tenx, color_by = "CLUSTER", group_by = NULL,
                        xlabel = "UMAP_1", ylabel = "UMAP_2")
p2 <- own_facet_scatter(plotx, "tSNE_1", "tSNE_2", title = project_name,isfacet = F,
                        col=color_conditions$tenx, color_by = "CLUSTER", group_by = NULL,
                        xlabel = "UMAP_1", ylabel = "UMAP_2")

somePDFPath = paste(allout,"07.3SCA_PLOT_UMAP_AUTOCLUSTER_",integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=14, height=7,pointsize=10)
print(p1)
dev.off()

somePDFPath = paste(allout,"07.4SCA_PLOT_TSNE_AUTOCLUSTER_",integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=14, height=7,pointsize=10)
print(p2)
dev.off()

clu_ann <- SingleR(test = as.SingleCellExperiment(data), ref = hpca.se, assay.type.test=1,labels = hpca.se$label.main)
write.csv(data.frame(clu_ann), paste(allout, "08.0SCA_TABLE_AUTO_ANNOTATION_CELL_TYPE_CLUSTERS_",integration_name,".csv", sep = ""), quote = F, row.names = F)

data$Cell_Type <- clu_ann$pruned.labels
data@meta.data[which(is.na(data$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
Idents(data) <- data$Cell_Type
plotx <- gen10x_plotx(data)
plotx$CELL_TYPE <- data$Cell_Type

somePDFPath = paste(allout,"08.1SCA_PLOT_UMAP_AUTO_CELL_TYPE_IDENTIFICATION_",integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=24, height=13,pointsize=10)
print(plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = project_name,
col = color_conditions$tenx, annot = TRUE, legend_position = "right"))
dev.off()

somePDFPath = paste(allout,"08.2SCA_PLOT_TSNE_AUTO_CELL_TYPE_IDENTIFICATION_",integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=24, height=13,pointsize=10)
print(plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CELL_TYPE", plot_title = project_name,
                   col = color_conditions$tenx, annot = TRUE, legend_position = "right"))
dev.off()

current_out <- deanalysis(data, current_clusters, plot_title = project_name,group=NULL,de_analysis = "findallmarkers")
current_data_markers <- current_out$current_data_markers
de_type <- current_out$de_type
de_name <- current_out$de_name
top1 <- current_out$top1
topn <- current_out$topn
current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
write.csv(current_data_markers, paste(allout, "07.5SCA_TABLE_AUTO_CLUSTER_TOP_MARKERS_",integration_name,".csv", sep = ""), quote = F, row.names = F)

current <- group_medianexpr(current_data_markers, data, ref_group = integration_cluster, group = "seurat_clusters", cell_type = F)
plot_median <- current$plot_median
top_markers <- current$top_markers

somePDFPath = paste(allout,"08.3SCA_PLOT_HEATMAP_MEDIAN_EXPR_CLUSTER_TOP_FOLDCHANGE_SIG_GENES_", integration_name, ".pdf", sep = "")
pdf(file=somePDFPath, width=10, height=n*1.2,pointsize=10)
heatmap.2(plot_median,margin=c(15, 15), trace="none",key=T, keysize=1,
          col=jet2.col (n = 100, alpha = 1),main = project_name,
          srtCol=0,
          scale="none", Colv = T, Rowv = T,
          density.info="none", cexCol=2,cexRow=2)
dev.off()

current <- data@meta.data
total_counts <- data.frame(table(current[,integration_cluster]))

current <- data.frame(table(current[,c("SAMPLE_ID",integration_cluster)]))
current <- current[current$Freq != 0,]
colnames(current) <- c("SAMPLE_ID","CLUSTER","COUNT")
current$CLUSTER_COUNT <- total_counts[match(current$CLUSTER, total_counts$Var1),"Freq"]
current$PROPORTION <- current$COUNT/current$CLUSTER_COUNT

node_proportion <- current

p <- ggplot(node_proportion, aes(CLUSTER, PROPORTION, fill = SAMPLE_ID))+
  geom_bar(stat="identity", alpha=0.8)+
  coord_polar()+
  scale_fill_viridis(discrete = T)+
  ggtitle(paste("Frequency of Samples in Each Node\n", project_name, sep = ""))+
  theme_bw(base_size = 28)+
  theme(plot.margin = unit(c(3,3,3,3), "cm"),
        plot.title = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title=element_text(size=30, face = "bold"), 
        legend.text=element_text(size=20),
        legend.key.size = unit(2, 'lines'),
        axis.title.x = element_text(colour="black", size = 30, face = "bold", vjust = -10),
        axis.title.y = element_text(colour="black", size = 30, face = "bold", vjust = 10),
        strip.text = element_text(size = 30, face = "bold"),
        axis.text.x=element_text(colour="black", size = 30),
        axis.text.y=element_text(colour="black", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

somePDFPath = paste(allout,"08.7SCA_PLOT_NODE_SUMMARY_SAMPLE_PROPORTION_",integration_name, ".pdf", sep = "")
pdf(somePDFPath, width = 25, height = 20, pointsize = 10)
print(p)
dev.off()

current <- group_medianexpr(current_data_markers, data, group = "Cell_Type", cell_type = T)
plot_median_cell_type <- current$plot_median
top_markers_cell_type <- current$top_markers

somePDFPath = paste(allout,"08.8SCA_PLOT_HEATMAP_MEDIAN_EXPR_CELL_TYPE_TOP_FOLDCHANGE_SIG_GENES_", integration_name,".pdf", sep = "")
pdf(file=somePDFPath, width=12, height=n*1.4,pointsize=10)
heatmap.2(plot_median_cell_type,margin=c(30, 15), trace="none",key=T, keysize=1,
          col=jet2.col (n = 100, alpha = 1),main = project_name,
          srtCol=45,
          scale="none", Colv = T, Rowv = T,
          density.info="none", cexCol=2,cexRow=2)
dev.off()

current <- data@meta.data
total_counts <- data.frame(table(current$Cell_Type))

current <- data.frame(table(current[,c("SAMPLE_ID","Cell_Type")]))
current <- current[current$Freq != 0,]
colnames(current) <- c("SAMPLE_ID","CELL_TYPE","COUNT")
current$CELL_TYPE_COUNT <- total_counts[match(current$CELL_TYPE, total_counts$Var1),"Freq"]
current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

node_proportion <- current

p <- ggplot(node_proportion, aes(CELL_TYPE, PROPORTION, fill = SAMPLE_ID))+
  geom_bar(stat="identity", alpha=0.8)+
  coord_polar()+
  scale_fill_viridis(discrete = T)+
  ggtitle(paste("Frequency of Samples in Each Node\n", project_name, sep = ""))+
  theme_bw(base_size = 28)+
  theme(plot.margin = unit(c(3,3,3,3), "cm"),
        plot.title = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title=element_text(size=30, face = "bold"), 
        legend.text=element_text(size=20),
        legend.key.size = unit(2, 'lines'),
        axis.title.x = element_text(colour="black", size = 30, face = "bold", vjust = -10),
        axis.title.y = element_text(colour="black", size = 30, face = "bold", vjust = 10),
        strip.text = element_text(size = 30, face = "bold"),
        axis.text.x=element_text(colour="black", size = 30),
        axis.text.y=element_text(colour="black", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

somePDFPath = paste(allout,"09.2SCA_PLOT_CELL_TYPE_SAMPLE_PROPORTION_",integration_name, ".pdf", sep = "")
pdf(somePDFPath, width = 25, height = 20, pointsize = 10)
print(p)
dev.off()

Idents(data) <- integration_cluster

save.image(paste(allout,"/INTEGRATED_POST_09.2.RData",sep = ""))


























