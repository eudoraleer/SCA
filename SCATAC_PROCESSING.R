# module load miniconda3/4.5.4

library("Seurat", lib="/proj/eudoraleer/")
Packages <- c("Signac","patchwork","ggpubr","dplyr","liftOver",
              "ggplot2","GenomeInfoDb", "gridExtra","motifmatchr",
              "JASPAR2020","TFBSTools","cowplot","SeuratWrappers",
              "EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38",
              "EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19",
              "monocle3","cicero")

suppressMessages(lapply(Packages, library, character.only = TRUE))
source('/proj/store_eudoraleer/PROJECTS/scRNASEQ/SCRIPTS/SCA_Analysis_Functions_V1.0.0.R')
# source('~/Dropbox/KI/Studies/LXXLZ/Website/Web_Scripts/SCA_Analysis_Functions_V1.0.0.R')
# macs_path <- "/Users/luadmpan/MyPythonEnv/bin/macs3"

color_conditions <- color_ini()

library(Seurat, lib="/proj/eudoraleer/")
library(Signac, lib="/proj/eudoraleer/")
library("liftOver")

# Summary + Filtering
input_dir <- "/proj/lxxlxx_storage/PROJECTS/TISSUEOMICS/scATAC/RDS/"
setwd(input_dir)

files <- list.files(input_dir, pattern = "\\.RDS")
summary <- NULL

for(i in 1:length(files)){
  x <- readRDS(files[i])
  DefaultAssay(x) <- "peaks"
  print(paste(i,":",nrow(x@meta.data)))
  
  summary <- rbind(summary, data.frame(FILE = files[i],
                                       ASSEMBLY = unique(x$REF_GENOME),
                                       NUM_CELLS_RAW = nrow(x@meta.data)))
}

write.table(summary,"/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/scATAC/DATA/SUMMARY_RAW_CELL_COUNTS.txt", quote=F, row.names=F)








project_name <- "scATAC_ALL_TISSUES"
input_dir <- "/proj/lxxlxx_storage/PROJECTS/TISSUEOMICS/scATAC/RDS/"
email <- "lu.pan@ki.se"
output_dir = "~/Desktop/atatc_pbmcs/Output/"
processed <- "YES"
phenodata_dir <- "~/Desktop/atatc_pbmcs/SCA_METADATA_atac_mix.csv"
batch_correction_method <- "SEURAT_V3"
time_series <- "YES"
pseudo_time_para <- "cell_type2"
integrate_rnaseq <- "YES"
rna_seq_dir <- "~/Desktop/atatc_pbmcs/"

dir.create(output_dir, recursive = T)
prefix_name <- gsub("(.*)_20[0-9]+","\\1",project_name)

# sink(paste(output_dir,"/",project_name,".log",sep = ""), append=TRUE)
Sys.setenv(TZ="Europe/Stockholm")

# check_zip(input_dir, type = "ATAC")

pheno_data <- pheno_ini(phenodata_dir, pipeline = "ATAC")
#####################################################################################################

pheno_singlecell <- list.files(input_dir, pattern = ".*singlecell\\.csv$|.*per_barcode_metrics\\.csv$", full.names = T, ignore.case = T, recursive = T)
pheno_fragments <- list.files(input_dir, pattern = ".*fragment.*tsv.*", full.names = T, ignore.case = T, recursive = T)
pheno_fragments <- pheno_fragments[grep("\\.tbi", pheno_fragments, ignore.case = T, invert = T)]

error_files <- list.files(input_dir, pattern = "gz\\.tsv$", full.names = T, ignore.case = T)

if(length(error_files) > 0){
  for(i in 1:length(error_files)){
    current_corrected_file <- gsub("gz\\.tsv","gz",error_files[i])
    system(paste("mv ", error_files[i], " ",current_corrected_file, sep = ""))
  }
}

data <- NULL
sample_names <- NULL

for(i in 1:nrow(pheno_data)){
  
  current <- Read10X_h5(filename = paste(input_dir, "/", pheno_data[i,"FILE"], sep = ""))
  current_singlecell <- NULL
  current_singlecell <- pheno_singlecell[grep(gsub("(.*)_filtered_.*|(.*)_raw.*","\\1", pheno_data[i,"FILE"],ignore.case = T), pheno_singlecell, ignore.case = T)]
  if(length(current_singlecell) > 0){
    current_meta <- "YES"
    current_singlecell <- read.csv(file = current_singlecell, header = TRUE, row.names = 1)
  }else{
    current_meta <- "NO"
  }
  
  if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    ref_genome <- "hg19"
    load("/Users/luadmpan/Dropbox/KI/Studies/LXXLZ/Website/10X/ATAC/Data/hg19_EnsDb.Hsapiens.v75.RData")
    seqlevelsStyle(ref_annot) <- "UCSC"
    genome(ref_annot) <- "hg19"
  }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    ref_genome <- "GRCh38.p12"
    load("/Users/luadmpan/Dropbox/KI/Studies/LXXLZ/Website/10X/ATAC/Data/hg38_EnsDb.Hsapiens.v86.RData")
    seqlevelsStyle(ref_annot) <- "UCSC"
    genome(ref_annot) <- "GRCh38.p12"
  }else{
    stop("You are supplying non-human samples. Please try again.")
  }
  
  
  ################################################################################################################
  
  DefaultAssay(current_seurat) <- 'peaks'
  Annotation(current_seurat) <- ref_annot
  current_seurat <- NucleosomeSignal(object = current_seurat)
  current_seurat <- TSSEnrichment(object = current_seurat, fast = FALSE)
  current_seurat$high.tss <- ifelse(current_seurat@meta.data[,grep("TSS.enrichment",colnames(current_seurat@meta.data), ignore.case = T)] > 2, 'High', 'Low')
  
  somePDFPath = paste(output_dir,"01.3SCA_PLOT_NORM_TSS_ENRICHMENT_SCORES_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=5,pointsize=10) 
  print(TSSPlot(current_seurat, group.by = 'high.tss') + NoLegend()+ scale_color_manual(values = color_conditions$alternate))
  dev.off()
  
  somePDFPath = paste(output_dir,"01.4SCA_PLOT_FRAGMENT_LENGTH_PERIODICITY_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=5,pointsize=10)
  current_col <- grep("nucleosome_signal",colnames(current_seurat@meta.data), ignore.case = T)
  current_seurat$nucleosome_group <- ifelse(current_seurat@meta.data[,current_col] > 4, 'NS > 4', 'NS < 4')
  current_groups <- unique(current_seurat$nucleosome_group)
  if(length(current_groups) == 1){
    print(FragmentHistogram(object = current_seurat) + 
            scale_fill_manual(values = color_conditions$tenx)+ggtitle(paste("FRAGMENT LENGTH PERIODICITY: ", current_groups, sep = "")))
  }else{
    print(FragmentHistogram(object = current_seurat, group.by = 'nucleosome_group') + 
            scale_fill_manual(values = color_conditions$tenx)+ggtitle("FRAGMENT LENGTH PERIODICITY"))
  }
  dev.off()
  
  if(current_meta == "YES" & toupper(data_type) != toupper("multiome")){
    colnames(current_seurat@meta.data)[grep("pct_reads_in_peaks", colnames(current_seurat@meta.data), ignore.case = T)] <- "PCT_Reads_in_Peaks"
    colnames(current_seurat@meta.data)[grep("blacklist_ratio", colnames(current_seurat@meta.data), ignore.case = T)] <- "Blacklist_Ratio"
  }
  
  colnames(current_seurat@meta.data)[grep("peak_region_fragments", colnames(current_seurat@meta.data), ignore.case = T)] <- "Peak_Region_Fragments"
  colnames(current_seurat@meta.data)[grep("TSS.enrichment", colnames(current_seurat@meta.data), ignore.case = T)] <- "TSS_Enrichment"
  colnames(current_seurat@meta.data)[grep("nucleosome_signal", colnames(current_seurat@meta.data), ignore.case = T)] <- "Nucleosome_Signal"
  
  if(current_meta == "YES" & toupper(data_type) != toupper("multiome")){
    current_features = c('PCT_Reads_in_Peaks', 'Peak_Region_Fragments','TSS_Enrichment', 'Blacklist_Ratio', 'Nucleosome_Signal')
  }else{
    current_features = c('Peak_Region_Fragments','TSS_Enrichment', 'Nucleosome_Signal')
  }
  
  somePDFPath = paste(output_dir,"01.5SCA_PLOT_QUALITY_CONTROL_PEAKS_LOGSCALE_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=16, height=8,pointsize=10)
  print(violin_plot(current_seurat, features = current_features, ncol = length(current_features),
                    col = color_conditions$tenx, x_lab = "SAMPLE_ID", log_status = TRUE))
  dev.off()
  
  if(toupper(data_type) == toupper("multiome")){
    current_seurat <- subset(
      x = current_seurat,
      subset = nCount_peaks < 100000 & nCount_peaks > 1000 &
        nCount_RNA < 25000 & nCount_RNA > 1000 &
        nFeature_RNA > 200 & nFeature_RNA < 2500 & 
        Nucleosome_Signal < 4 &
        TSS_Enrichment > 2 &
        Peak_Region_Fragments > 3000 &
        Peak_Region_Fragments < 20000 &
        percent.mt < 5
    )
    
    DefaultAssay(current_seurat) <- 'RNA'
    somePDFPath = paste(output_dir,"01.6SCA_PLOT_POSTQC_RNASEQ_LOGSCALE_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width = 16, height = 8,pointsize=10) 
    print(violin_plot(current_seurat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), col = color_conditions$tenx, x_lab = "SAMPLE_ID", log_status = TRUE))
    dev.off()
    
    somePDFPath = paste(output_dir,"01.7SCA_PLOT_POSTQC_RNASEQ_SCATTERED_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=10,pointsize=10) 
    plot1 <- FeatureScatter(current_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = color_conditions$warm)
    plot2 <- FeatureScatter(current_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = color_conditions$cold)
    print(plot1 + plot2)
    dev.off()
    
  }else{
    current_seurat <- subset(
      x = current_seurat,
      subset = Peak_Region_Fragments > 3000 &
        Peak_Region_Fragments < 20000 &
        nCount_peaks < 100000 &
        nCount_peaks > 1000 &
        PCT_Reads_in_Peaks > 15 &
        Blacklist_Ratio < 0.05 &
        Nucleosome_Signal < 4 &
        TSS_Enrichment > 2)
  }
  
  DefaultAssay(current_seurat) <- "peaks"
  peaks <- CallPeaks(current_seurat, macs2.path = macs_path)
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg19, invert = TRUE)
  }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  }
  
  macs_counts <- FeatureMatrix(fragments = Fragments(current_seurat), features = peaks, cells = colnames(current_seurat))
  rm(peaks)
  
  current_seurat[["peaks"]] <- CreateChromatinAssay(
    counts = macs_counts,
    fragments = pheno_fragments[grep(gsub("(.*)_filtered_.*|(.*)_raw.*","\\1",
                                          pheno_data[i,"FILE"],ignore.case = T), pheno_fragments, ignore.case = T)],
    annotation = ref_annot)
  
  rm(macs_counts)
  rm(ref_annot)
  
  DefaultAssay(current_seurat) <- "peaks"
  current_seurat <- RunTFIDF(current_seurat)
  current_seurat <- FindTopFeatures(current_seurat, min.cutoff = 'q0')
  current_seurat <- RunSVD(current_seurat)
  
  somePDFPath = paste(output_dir,"01.8SCA_PLOT_CORRELATION_SEQUENCING_DEPTH_LSI_COMPONENTS_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=5,pointsize=10) 
  print(DepthCor(current_seurat))
  dev.off()
  
  dim_i <- 2 # Skip first component with high correlation
  current_seurat <- RunUMAP(object = current_seurat, reduction = 'lsi', dims = dim_i:30)
  current_seurat <- RunTSNE(object = current_seurat, reduction = 'lsi', dims = dim_i:30, check_duplicates = FALSE)
  current_seurat <- FindNeighbors(object = current_seurat, reduction = 'lsi', dims = dim_i:30)
  current_seurat <- FindClusters(object = current_seurat, verbose = FALSE, algorithm = 3)
  
  plotx <- data.frame(UMAP_1 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_1"],
                      UMAP_2 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_2"],
                      tSNE_1 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_1"],
                      tSNE_2 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_2"],
                      Cluster = current_seurat$seurat_clusters)
  
  plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.numeric(as.character(unlist(plotx$Cluster))))))
  p1_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "scATAC UMAP", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
  p1_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "scATAC tSNE", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
  
  p_thresh <- 0.05
  DefaultAssay(current_seurat) <- "peaks"
  da_peaks <- FindAllMarkers(
    object = current_seurat,
    min.pct = 0.2,
    test.use = 'LR',
    latent.vars = 'Peak_Region_Fragments', return.thresh = p_thresh)
  da_peaks <- da_peaks[da_peaks$p_val_adj < p_thresh,]
  da_peaks <- da_peaks[order(da_peaks$p_val_adj, decreasing = FALSE),]
  write.csv(da_peaks, paste(output_dir, "01.9SCA_TABLE_TOP_PEAKS_IN_EACH_CLUSTERS_PADJUST_",p_thresh, "_", pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
  
  n <- 6
  wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]
  topn <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:n) # %>% filter(avg_log2FC > 1)
  
  clusters <- sort(as.numeric(as.character(unique(topn$cluster))),decreasing = F)
  somePDFPath = paste(output_dir,"02.0SCA_PLOT_ATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=16, height=10,pointsize=10)
  
  for(k in 1:length(clusters)){
    
    print(VlnPlot(
      object = current_seurat,
      features = unlist(topn[topn$cluster == clusters[k],"gene"]),
      pt.size = 0.1,
      idents = sort(as.numeric(as.character(current_seurat$seurat_clusters))),
      cols = gen_colors(color_conditions$tableau20, length(unique(current_seurat$seurat_clusters)))) & xlab("Clusters") &
        plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))
  }
  
  dev.off()
  
  somePDFPath = paste(output_dir,"02.1SCA_PLOT_ATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_COLOR_PEAKS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=16, height=10,pointsize=10)
  
  for(k in 1:length(clusters)){
    print(VlnPlot(
      object = current_seurat,
      features = unlist(topn[topn$cluster == clusters[k],"gene"]),
      pt.size = 0.1, stack = T,
      idents = sort(as.numeric(as.character(current_seurat$seurat_clusters))),
      flip = T, cols = gen_colors(color_conditions$tableau20, n)) +
        theme(legend.position = "none") +
        ggtitle(paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k],sep = "")) & xlab("Clusters"))
  }
  
  dev.off()
  
  somePDFPath = paste(output_dir,"02.2SCA_PLOT_ATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_FEATURE_PLOT_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=10,pointsize=10)
  
  for(k in 1:length(clusters)){
    print(FeaturePlot(
      object = current_seurat, label = T,
      features = unlist(topn[topn$cluster == clusters[k],"gene"]),
      cols = c("green", "blue"),
      pt.size = 0.1)+
        plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))
  }
  
  dev.off()
  
  top1 <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
  top1 <- top1[order(top1$cluster, decreasing = F),]
  
  current_plots <- NULL
  
  for(k in 1:length(top1$gene)){
    
    current_plots[[k]] <- FeaturePlot(current_seurat, features = top1$gene[k], cols = c("green", "blue"),
                                      label = T, pt.size = 0.1, max.cutoff = 'q95')+
      ggtitle(paste("PEAK: ",top1$gene[k],"\nCLUSTER: ", top1$cluster[k], sep = ""))
  }
  
  somePDFPath = paste(output_dir,"02.3SCA_PLOT_ATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_FEATURE_PLOT_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=ifelse(length(current_plots) < 10,8,14),pointsize=10)
  print(do.call("grid.arrange", c(current_plots, ncol=4)))
  dev.off()
  
  top_FC <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5)))
  
  for(k in 1:length(clusters)){
    current <- unlist(top_FC[top_FC$cluster == clusters[k],"gene"])
    closest_genes <- ClosestFeature(current_seurat, regions = current)
    write.csv(closest_genes, paste(output_dir, "02.4SCA_TABLE_TOP_PEAKS_CLOSEST_GENOMIC_REGIONS_CLUSTER_",clusters[k],"_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
  }
  
  ###########################################################################################################
  
  if((length(grep("DIFFERENT.*CELL", pheno_data[i,"DATA_TYPE"], ignore.case = T)) > 0) |
     pheno_data[i,"DATA_TYPE"] %in% toupper(c("","-","NA","NULL")) |
     is.null(pheno_data[i,"DATA_TYPE"]) | is.na(pheno_data[i,"DATA_TYPE"])){
    
    if(pheno_data[i,"DATA_TYPE"] %in% toupper(c("","-","NA","NULL")) |
       is.null(pheno_data[i,"DATA_TYPE"]) | is.na(pheno_data[i,"DATA_TYPE"])){
      somePDFPath = paste(output_dir,"02.5SCA_PLOT_UMAP_tSNE_scATAC-SEQ_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
      print(p1_umap+p1_tsne)
      dev.off()
    }
    
    gene_activities <- GeneActivity(current_seurat)
    current_seurat[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
    DefaultAssay(current_seurat) <- "ACTIVITY"
    current_seurat <- NormalizeData(current_seurat)
    current_seurat <- ScaleData(current_seurat, features = rownames(current_seurat))
    current_data_markers <- FindAllMarkers(current_seurat, min.pct = 0.25, logfc.threshold = 0.25)
    p_thresh <- 0.05
    current_data_markers <- current_data_markers[current_data_markers$p_val_adj < p_thresh,]
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
    write.csv(current_data_markers, paste(output_dir, "02.6SCA_TABLE_GENE_ACTIVITY_MATRIX_CLUSTERS_COMPARISON_SIG_PADJ",p_thresh,"_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
    
    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    top1 <- current_data_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    current_plots <- NULL
    
    for(k in 1:length(top1$gene)){
      
      current_plots[[k]] <- FeaturePlot(current_seurat, features = top1$gene[k], cols = c("green", "blue"),
                                        label = T, pt.size = 0.1, max.cutoff = 'q95')+
        ggtitle(paste("CLUSTER: ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""))
    }
    
    somePDFPath = paste(output_dir,"02.7SCA_PLOT_GENE_ACTIVITY_MATRIX_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=ifelse(length(top1$gene) < 10, 16, 20), height=ifelse(length(top1$gene) < 10, 8, 16),pointsize=10)
    print(do.call("grid.arrange", c(current_plots, ncol=4)))
    dev.off()
    
    n <- 10
    top_n <- current_data_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    
    current_clusters <- as.numeric(as.character(sort(unique(current_data_markers$cluster))))
    somePDFPath = paste(output_dir,"02.8SCA_PLOT_GENE_ACTIVITY_MATRIX_TOP_",n,"_MARKERS_IN_EACH_CLUSTER_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=14,pointsize=10)
    
    for(j in 1:length(current_clusters)){
      current <- top_n[which(top_n$cluster == current_clusters[j]),]
      current <- current[order(current$p_val_adj, decreasing = F),]
      print(VlnPlot(current_seurat, features = current$gene, pt.size = 0,
                    # log = TRUE,
                    ncol = 4, cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) & xlab("Clusters") &
              plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", current_clusters[j], sep = ""),
                              theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    
    DefaultAssay(current_seurat) <- "ACTIVITY"
    n <- 10
    top1 <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    top1_derivedrna <- current_data_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1_derivedrna <- top1_derivedrna[order(top1_derivedrna$cluster, decreasing = F),]
    
    somePDFPath = paste(output_dir,"02.9SCA_PLOT_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_ACTIVITY_GENE_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=7,pointsize=10)
    
    DefaultAssay(current_seurat) <- "peaks"
    
    for(k in 1:length(clusters)){
      print(CoveragePlot(
        object = current_seurat,
        region = unlist(top1$gene[k]), features = unlist(top1_derivedrna[which(top1_derivedrna$cluster == clusters[k]),"gene"]),
        extend.upstream = 10000,
        extend.downstream = 10000,
        expression.assay = "ACTIVITY")+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    
    DefaultAssay(current_seurat) <- "peaks"
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19)
    }else if(length(grep("hg38|grch38|38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
    }
    
    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    top1 <- current_data_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    current_seurat <- LinkPeaks(
      object = current_seurat,
      peak.assay = "peaks",
      expression.assay = "ACTIVITY",
      genes.use = unlist(top1$gene))
    
    clusters <- sort(as.numeric(as.character(unique(current_seurat$seurat_clusters))))
    somePDFPath = paste(output_dir,"03.0SCA_PLOT_ATAC_Tn5_INSERTION_COVERAGE_LINK_PEAK_TOP1_ACTIVITY_GENE_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=9,pointsize=10)
    
    DefaultAssay(current_seurat) <- "peaks"
    
    for(k in 1:length(clusters)){
      print(CoveragePlot(
        object = current_seurat,
        region = unlist(top1$gene[k]), features = unlist(top1$gene[k]),
        expression.assay = "ACTIVITY",
        extend.upstream = 10000,
        extend.downstream = 10000)+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    
    if(!((pheno_data[i,"RNASEQ_FILE"] %in% toupper(c("","-","NA","NULL"))) |
         is.null(pheno_data[i,"RNASEQ_FILE"]) | is.na(pheno_data[i,"RNASEQ_FILE"]))){
      
      current_rnaseq <- Read10X_h5(filename = paste(input_dir, "/", pheno_data[i,"RNASEQ_FILE"], sep = ""))
      
      if((length(current_rnaseq) == 1) | (length(current_rnaseq) > 10)){
        current_rnaseq <- CreateSeuratObject(counts = current_rnaseq, project = prefix_name, min.cells = 3, min.features = 200)
      }else{
        current_rnaseq <- CreateSeuratObject(counts = current_rnaseq$`Gene Expression`, project = prefix_name, min.cells = 3, min.features = 200)
      }
      
      current_rnaseq[["percent.mt"]] <- PercentageFeatureSet(current_rnaseq, pattern = "^MT-")
      
      somePDFPath = paste(output_dir,"03.1SCA_PLOT_QC_RNASEQ_LOGSCALE_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=10,pointsize=10)
      print(violin_plot(current_rnaseq, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        col = color_conditions$tenx, x_lab = "SAMPLE_ID", log_status = TRUE))
      dev.off()
      
      somePDFPath = paste(output_dir,"03.2SCA_PLOT_QC_RNASEQ_SCATTERED_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=10,pointsize=10) 
      plot1 <- FeatureScatter(current_rnaseq, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = color_conditions$warm)
      plot2 <- FeatureScatter(current_rnaseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = color_conditions$cold)
      print(plot1 + plot2)
      dev.off()
      
      current_rnaseq <- subset(
        x = current_rnaseq,
        nCount_RNA < 25000 & nCount_RNA > 1000 &
          nFeature_RNA > 200 & nFeature_RNA < 2500 &
          percent.mt < 5)
      
      somePDFPath = paste(output_dir,"03.3SCA_PLOT_POSTQC_RNASEQ_LOGSCALE_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=6,pointsize=10) 
      print(violin_plot(current_rnaseq, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        col = color_conditions$tenx, x_lab = "SAMPLE_ID", log_status = TRUE))
      dev.off()
      
      somePDFPath = paste(output_dir,"03.4SCA_PLOT_POSTQC_RNASEQ_SCATTERED_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=10,pointsize=10) 
      plot1 <- FeatureScatter(current_rnaseq, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = color_conditions$warm)
      plot2 <- FeatureScatter(current_rnaseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = color_conditions$cold)
      print(plot1 + plot2)
      dev.off()
      
      current_rnaseq <- suppressWarnings(NormalizeData(current_rnaseq, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE))
      current_rnaseq <- FindVariableFeatures(current_rnaseq)
      current_rnaseq <- ScaleData(current_rnaseq, verbose = FALSE, features = VariableFeatures(current_rnaseq), vars.to.regress = c("nCount_RNA", "percent.mt"))
      current_rnaseq <- RunPCA(current_rnaseq, verbose = FALSE)
      current_rnaseq <- RunTSNE(current_rnaseq, reduction = "pca", dims = 1:30, check_duplicates = FALSE)
      current_rnaseq <- RunUMAP(current_rnaseq, reduction = "pca", dims = 1:30)
      current_rnaseq <- FindNeighbors(current_rnaseq, dims = 1:10)
      current_rnaseq <- FindClusters(current_rnaseq, resolution = 0.5, n.iter = 1000)
      current_rnaseq_markers <- FindAllMarkers(current_rnaseq, min.pct = 0.25, logfc.threshold = 0.25)
      
      p_thresh <- 0.05
      current_rnaseq_markers <- current_rnaseq_markers[current_rnaseq_markers$p_val_adj < p_thresh,]
      current_rnaseq_markers <- current_rnaseq_markers[order(current_rnaseq_markers$p_val_adj, decreasing = F),]
      write.csv(current_rnaseq_markers, paste(output_dir, "03.5SCA_TABLE_GIVEN_RNA_SEQ_AUTO_CLUSTER_TOP_MARKERS_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
      
      wt <- colnames(current_rnaseq_markers)[grep("log.*FC", colnames(current_rnaseq_markers), ignore.case = T)]
      top1 <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
      top1 <- top1[order(top1$cluster, decreasing = F),]
      
      current_plots <- NULL
      
      for(k in 1:length(top1$gene)){
        
        current_plots[[k]] <- FeaturePlot(current_rnaseq, features = top1$gene[k], cols = c("green", "blue"),
                                          label = T, pt.size = 0.1, max.cutoff = 'q95')+
          ggtitle(paste("CLUSTER: ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""))
      }
      
      somePDFPath = paste(output_dir,"03.6SCA_PLOT_GIVEN_RNA_SEQ_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=length(unique(current_rnaseq$seurat_clusters))/4*5,pointsize=10)
      print(do.call("grid.arrange", c(current_plots, ncol=4)))
      dev.off()
      
      n <- 10
      top_n <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
      
      current_clusters <- as.numeric(as.character(sort(unique(current_rnaseq_markers$cluster))))
      somePDFPath = paste(output_dir,"03.7SCA_PLOT_GIVEN_RNA_SEQ_TOP_",n,"_MARKERS_IN_EACH_CLUSTER_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=14,pointsize=10)
      
      for(j in 1:length(current_clusters)){
        current <- top_n[which(top_n$cluster == current_clusters[j]),]
        current <- current[order(current$p_val_adj, decreasing = F),]
        print(VlnPlot(current_rnaseq, features = current$gene, pt.size = 0,
                      # log = TRUE,
                      ncol = 4, cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) & xlab("Clusters") &
                plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", current_clusters[j], sep = ""),
                                theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"))))      }
      
      dev.off()
      
      plotx <- data.frame(UMAP_1 = current_rnaseq@reductions$umap@cell.embeddings[,"UMAP_1"],
                          UMAP_2 = current_rnaseq@reductions$umap@cell.embeddings[,"UMAP_2"],
                          tSNE_1 = current_rnaseq@reductions$tsne@cell.embeddings[,"tSNE_1"],
                          tSNE_2 = current_rnaseq@reductions$tsne@cell.embeddings[,"tSNE_2"],
                          Cluster = current_rnaseq$seurat_clusters,
                          Type = "scRNA-Seq")
      
      plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.numeric(as.character(unlist(plotx$Cluster))))))
      p2_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "scRNA-Seq UMAP", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
      p2_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "scRNA-Seq tSNE", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
      
      somePDFPath = paste(output_dir,"03.8SCA_PLOT_UMAP_GIVEN_scRNASEQ_VS_scATAC-SEQ_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
      print(p1_umap+p2_umap)
      dev.off()
      
      somePDFPath = paste(output_dir,"03.9SCA_PLOT_tSNE_GIVEN_scRNASEQ_VS_scATAC-SEQ_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
      print(p1_tsne+p2_tsne)
      dev.off()
      
    }
  }
  
  if(toupper(data_type) == toupper("multiome")){
    
    DefaultAssay(current_seurat) <- "RNA"
    current_seurat <- suppressWarnings(SCTransform(current_seurat, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE))
    current_seurat <- FindVariableFeatures(current_seurat)
    current_seurat <- ScaleData(current_seurat, verbose = FALSE, features = VariableFeatures(current_seurat), vars.to.regress = c("nCount_RNA", "percent.mt"))
    current_seurat <- RunPCA(current_seurat, verbose = FALSE)
    current_rnaseq_markers <- FindAllMarkers(current_seurat, min.pct = 0.25, logfc.threshold = 0.25)
    
    p_thresh <- 0.05
    current_rnaseq_markers <- current_rnaseq_markers[current_rnaseq_markers$p_val_adj < p_thresh,]
    current_rnaseq_markers <- current_rnaseq_markers[order(current_rnaseq_markers$p_val_adj, decreasing = F),]
    write.csv(current_rnaseq_markers, paste(output_dir, "04.0SCA_TABLE_GIVEN_RNA_SEQ_AUTO_CLUSTER_TOP_MARKERS_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
    
    wt <- colnames(current_rnaseq_markers)[grep("log.*FC", colnames(current_rnaseq_markers), ignore.case = T)]
    top1 <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    current_plots <- NULL
    
    for(k in 1:length(top1$gene)){
      
      current_plots[[k]] <- FeaturePlot(current_seurat, features = top1$gene[k], cols = c("green", "blue"),
                                        label = T, pt.size = 0.1, max.cutoff = 'q95')+
        ggtitle(paste("CLUSTER: ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""))
    }
    
    somePDFPath = paste(output_dir,"04.1SCA_PLOT_GIVEN_RNA_SEQ_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=length(unique(current_seurat$seurat_clusters))/4*5,pointsize=10)
    print(do.call("grid.arrange", c(current_plots, ncol=4)))
    dev.off()
    
    n <- 10
    top_n <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    
    current_clusters <- as.numeric(as.character(sort(unique(current_rnaseq_markers$cluster))))
    somePDFPath = paste(output_dir,"04.2SCA_PLOT_GIVEN_RNA_SEQ_TOP_",n,"_MARKERS_IN_EACH_CLUSTER_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=14,pointsize=10)
    
    for(j in 1:length(current_clusters)){
      current <- top_n[which(top_n$cluster == current_clusters[j]),]
      current <- current[order(current$p_val_adj, decreasing = F),]
      print(VlnPlot(current_seurat, features = current$gene, pt.size = 0,
                    # log = TRUE,
                    ncol = 4, cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) +
              plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", current_clusters[j], sep = ""),
                              theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"))))      }
    
    dev.off()
    
    somePDFPath = paste(output_dir,"04.3SCA_PLOT_UMAP_tSNE_scATAC-SEQ_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p1_umap+p1_tsne)
    dev.off()
    
    dim_i <- 2
    current_seurat <- FindMultiModalNeighbors(
      object = current_seurat,
      reduction.list = list("pca", "lsi"), 
      dims.list = list(1:50, dim_i:40),
      modality.weight.name = "RNA.weight",
      verbose = TRUE
    )
    
    current_seurat <- RunUMAP(
      object = current_seurat,
      nn.name = "weighted.nn",
      assay = "RNA",
      verbose = TRUE)
    
    current_seurat <- RunTSNE(
      object = current_seurat,
      nn.name = "weighted.nn",
      assay = "RNA",
      check_duplicates = FALSE,
      verbose = TRUE)
    
    plotx <- data.frame(UMAP_1 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_2"],
                        tSNE_1 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_1"],
                        tSNE_2 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_2"],
                        Cluster = current_seurat$seurat_clusters,
                        Type = "Joint")
    
    p2_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "Joint scATAC-RNA UMAP", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    p2_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "Joint scATAC-RNA tSNE", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    
    somePDFPath = paste(output_dir,"04.4SCA_PLOT_UMAP_tSNE_Joint_scATAC_RNA_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p2_umap+p2_tsne)
    dev.off()
    
    DefaultAssay(current_seurat) <- "peaks"
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19)
    }else if(length(grep("hg38|grch38|38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
    }
    
    wt <- colnames(current_rnaseq_markers)[grep("log.*FC", colnames(current_rnaseq_markers), ignore.case = T)]
    top1 <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    current_seurat <- LinkPeaks(
      object = current_seurat,
      peak.assay = "peaks",
      expression.assay = "SCT",
      genes.use = unlist(top1$gene))
    
    clusters <- sort(as.numeric(as.character(unique(current_seurat$seurat_clusters))))
    somePDFPath = paste(output_dir,"04.5SCA_PLOT_ATAC_Tn5_INSERTION_COVERAGE_LINK_PEAK_TOP1_RNASEQ_GENE_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=7,pointsize=10)
    
    DefaultAssay(current_seurat) <- "peaks"
    
    for(k in 1:length(clusters)){
      print(CoveragePlot(
        object = current_seurat,
        region = unlist(top1$gene[k]), features = unlist(top1$gene[k]),
        expression.assay = "SCT",
        extend.upstream = 10000,
        extend.downstream = 10000)+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    
    DefaultAssay(current_seurat) <- "RNA"
    n <- 10
    top1 <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    top1_rna <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1_rna <- top1_rna[order(top1_rna$cluster, decreasing = F),]
    
    somePDFPath = paste(output_dir,"04.6SCA_PLOT_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_RNASEQ_GENE_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=7,pointsize=10)
    
    DefaultAssay(current_seurat) <- "peaks"
    
    for(k in 1:length(clusters)){
      print(CoveragePlot(
        object = current_seurat,
        region = unlist(top1$gene[k]), features = unlist(top1_rna[which(top1_rna$cluster == clusters[k]),"gene"]),
        extend.upstream = 10000,
        extend.downstream = 10000)+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    
  }
  
  ######################################### Motifs #############################################
  ######################################## Analysis ############################################
  n <- 6
  current_assay <- "peaks"
  
  if(toupper(data_type) == toupper("multiome")){
    topn_genes <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    used_type <- "RNASEQ"
  }else if((length(grep("DIFFERENT.*CELL", pheno_data[i,"DATA_TYPE"], ignore.case = T)) > 0) |
           pheno_data[i,"DATA_TYPE"] %in% toupper(c("","-","NA","NULL")) |
           is.null(pheno_data[i,"DATA_TYPE"]) | is.na(pheno_data[i,"DATA_TYPE"])){
    topn_genes <- current_rnaseq_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    used_type <- "RNASEQ"
    if(pheno_data[i,"DATA_TYPE"] %in% toupper(c("","-","NA","NULL")) |
       is.null(pheno_data[i,"DATA_TYPE"]) | is.na(pheno_data[i,"DATA_TYPE"])){
      topn_genes <- current_data_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
      used_type <- "ACTIVITY"
    }
    
  }
  DefaultAssay(current_seurat) <- "peaks"
  topn_genes <- topn_genes[order(topn_genes$cluster, decreasing = F),]
  position_freq <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
  
  if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = position_freq)
    motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
    if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
      chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
      chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
      current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                  genome = BSgenome.Hsapiens.UCSC.hg19)
    }
    
  }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = position_freq)
    motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
    if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
      chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
      chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
      current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                  genome = BSgenome.Hsapiens.UCSC.hg38)
    }
  }
  
  somePDFPath = paste(output_dir,"04.7SCA_PLOT_",toupper(current_assay),"_TOP_",used_type,"_GENES_FOOTPRINTING_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=18/1.25,pointsize=10)
  p1 <- PlotFootprint(current_seurat, features = chosen_genes, group.by	= "seurat_clusters")
  print(p1 + patchwork::plot_layout(ncol = 2)+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": TOP ",used_type," GENES \n(Top group with highest accessibility in motif flanking region are labeled)", sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
  dev.off()
  
  p_thresh <- 0.05
  DefaultAssay(current_seurat) <- "peaks"
  
  da_peaks <- FindAllMarkers(
    object = current_seurat,
    only.pos = TRUE,
    test.use = 'LR',
    latent.vars = 'nCount_peaks')
  
  top.da.peak <- unique(da_peaks[da_peaks$p_val_adj < p_thresh,"gene"])
  enriched.motifs <- FindMotifs(object = current_seurat,features = top.da.peak)
  write.csv(enriched.motifs, paste(output_dir, "02.7SCA_TABLE_ENRICHED_MOTIFS_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
  
  n <- 10
  somePDFPath = paste(output_dir,"04.8SCA_PLOT_TOP_GENES_ENRICHED_MOTIFS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8,pointsize=10)
  print(MotifPlot(object = current_seurat, motifs = rownames(enriched.motifs)[1:n])+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], "\nPOSITION WEIGHT MATRICES OF TOP ENRICHED MOTIFS", sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
  dev.off()
  
  # Per-cell motif activity score
  if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg19)
  }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
    current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg38)
  }
  
  n <- 6
  DefaultAssay(current_seurat) <- 'chromvar'
  p2_chromvar <- FeaturePlot(
    object = current_seurat,
    features = enriched.motifs$motif[1:n],
    min.cutoff = 'q10',
    max.cutoff = 'q90',
    pt.size = 0.1, label = T,
    cols = c("green", "blue")) +
    plot_annotation(title = "Motif Activities",
                    theme = theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold")))
  
  somePDFPath = paste(output_dir,"04.9SCA_PLOT_UMAP_VS_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=8,pointsize=10)
  print(p1_umap + p2_chromvar +
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], "\nPeaks VS Top Enriched Motifs", sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
  dev.off()
  
  write.csv(current_seurat@meta.data, paste(output_dir, "05.0SCA_TABLE_METADATA_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
  
  # DefaultAssay(current_seurat) <- "peaks"
  # current_cicero <- as.cell_data_set(current_seurat)
  # current_cicero <- cluster_cells(cds = current_cicero, reduction_method = "UMAP")
  # current_cicero <- learn_graph(current_cicero, use_partition = TRUE)
  
  ######################################### Integrate #############################################
  ######################################## Peaks + RNA-Seq ########################################  
  
  
  if((toupper(data_type) == toupper("individual")) &
     !((is.null(pheno_data[i,"DATA_TYPE"])) | (is.na(pheno_data[i,"DATA_TYPE"])) |
       pheno_data[i,"DATA_TYPE"] %in% toupper(c("","-","NA","NULL")))){
    
    DefaultAssay(current_seurat) <- "ACTIVITY"
    DefaultAssay(current_rnaseq) <- "RNA"
    transfer.anchors <- FindTransferAnchors(reference = current_rnaseq, query = current_seurat,
                                            features = VariableFeatures(object = current_rnaseq), 
                                            reference.assay = "RNA", query.assay = "ACTIVITY",
                                            reduction = "cca")
    predicted_labels <- TransferData(
      anchorset = transfer.anchors,
      refdata = current_rnaseq$seurat_clusters,
      weight.reduction = current_seurat[['lsi']],
      dims = 2:30)
    
    current_seurat <- AddMetaData(object = current_seurat, metadata = predicted_labels)
    current_seurat$annotation_correct <- current_seurat$predicted.id == current_seurat$seurat_clusters
    
    plotx <- data.frame(UMAP_1 = current_rnaseq@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = current_rnaseq@reductions$umap@cell.embeddings[,"UMAP_2"],
                        tSNE_1 = current_rnaseq@reductions$tsne@cell.embeddings[,"tSNE_1"],
                        tSNE_2 = current_rnaseq@reductions$tsne@cell.embeddings[,"tSNE_2"],
                        Cluster = current_rnaseq$seurat_clusters,
                        Type = "scRNA-seq")
    
    plotx <- rbind(plotx, data.frame(UMAP_1 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_1"],
                                     UMAP_2 = current_seurat@reductions$umap@cell.embeddings[,"UMAP_2"],
                                     tSNE_1 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                     tSNE_2 = current_seurat@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                     Cluster = current_seurat$predicted.id,
                                     Type = "scATAC-seq"))
    
    plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.numeric(as.character(unlist(plotx$Cluster))))))
    
    p3_umap <- plot_bygroup(plotx[grep("scRNA-seq", plotx$Type, ignore.case = T),], "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "scRNA-seq UMAP", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    p3_tsne <- plot_bygroup(plotx[grep("scRNA-seq", plotx$Type, ignore.case = T),], "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "scRNA-seq tSNE", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    p4_umap <- plot_bygroup(plotx[grep("scATAC-seq", plotx$Type, ignore.case = T),], "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "scATAC-seq UMAP\n(Predicted Clusters from scRNA-Seq)", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    p4_tsne <- plot_bygroup(plotx[grep("scATAC-seq", plotx$Type, ignore.case = T),], "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "scATAC-seq tSNE\n(Predicted Clusters from scRNA-Seq)", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    
    somePDFPath = paste(output_dir,"05.1SCA_PLOT_UMAP_scRNASEQ_VS_scATAC-SEQ_INTEGRATE_CLUST_PREDICTION_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p3_umap | p4_umap)
    dev.off()
    
    somePDFPath = paste(output_dir,"05.2SCA_PLOT_tSNE_scRNASEQ_VS_scATAC-SEQ_INTEGRATE_CLUST_PREDICTION_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p3_tsne | p4_tsne)
    dev.off()
    
    predictions <- table(current_seurat$seurat_clusters, current_seurat$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)
    p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", 
                                                                                                low = "green", high = "blue") + xlab("True Clusters") + ylab("Predicted Clusters") + 
      theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    correct <- length(which(current_seurat$seurat_clusters == current_seurat$predicted.id))
    incorrect <- length(which(current_seurat$seurat_clusters != current_seurat$predicted.id))
    fetchdata <- FetchData(current_seurat, vars = c("prediction.score.max", "annotation_correct"))
    p2 <- ggplot(fetchdata, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) + 
      geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct", 
                                                                        labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + 
      scale_color_discrete(name = "Annotation Correct", 
                           labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + 
      xlab("Prediction Score") +
      scale_color_manual(values = color_conditions$tenx) +
      scale_fill_manual(values = color_conditions$tenx)
    
    somePDFPath = paste(output_dir,"05.3SCA_PLOT_scRNASEQ_VS_scATAC-SEQ_INTEGRATE_PREDICTION_SCORES_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p1 + p2)
    dev.off()
    
    genes_used <- VariableFeatures(current_rnaseq)
    refdata <- GetAssayData(current_rnaseq, assay = "RNA", slot = "data")[genes_used, ]
    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.
    # imputation (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = current_seurat[["lsi"]], 
                               dims = 2:30)
    rm(refdata)
    current_seurat[["RNA"]] <- imputation
    rm(imputation)
    coembed <- merge(x = current_rnaseq, y = current_seurat)
    coembed <- ScaleData(coembed, features = genes_used, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes_used, verbose = FALSE)
    coembed <- RunUMAP(coembed, dims = 1:30)
    coembed <- RunTSNE(coembed, dims = 1:30, check_duplicates = FALSE)
    
    plotx <- data.frame(UMAP_1 = coembed@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = coembed@reductions$umap@cell.embeddings[,"UMAP_2"],
                        tSNE_1 = coembed@reductions$tsne@cell.embeddings[,"tSNE_1"],
                        tSNE_2 = coembed@reductions$tsne@cell.embeddings[,"tSNE_2"],
                        Data_Type = coembed$orig.ident,
                        Clusters = coembed$seurat_clusters)
    plotx$Clusters <- as.numeric(as.character(plotx$Clusters))
    plotx$Clusters <- factor(plotx$Clusters, levels = sort(unique(as.numeric(as.character(unlist(plotx$Clusters))))))
    
    p5_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Data_Type", plot_title = "Data Type - UMAP", col = color_conditions$tableau20[grep("black",color_conditions$tableau20, invert = TRUE)], annot = TRUE, legend_position = "right")
    p5_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Data_Type", plot_title = "Data Type - tSNE", col = color_conditions$tableau20[grep("black",color_conditions$tableau20, invert = TRUE)], annot = TRUE, legend_position = "right")
    p6_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Clusters", plot_title = "Coembed scATAC + scRNA-SEQ: UMAP", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    p6_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Clusters", plot_title = "Coembed scATAC + scRNA-SEQ: tSNE", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
    
    somePDFPath = paste(output_dir,"05.4SCA_PLOT_UMAP_GIVEN_scRNASEQ_VS_scATAC-SEQ_INTEGRATIVE_COEMBED_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p5_umap | p6_umap)
    dev.off()
    
    somePDFPath = paste(output_dir,"05.5SCA_PLOT_tSNE_GIVEN_scRNASEQ_VS_scATAC-SEQ_INTEGRATIVE_COEMBED_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
    print(p5_tsne | p6_tsne)
    dev.off()
    
    rm(coembed)
    saveRDS(current_rnaseq, paste(output_dir,"RNASEQ_SEURAT_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".RDS", sep = ""))
    
    n <- 10
    top1 <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
    top1 <- top1[order(top1$cluster, decreasing = F),]
    
    DefaultAssay(current_seurat) <- "RNA"
    current_imput_rnaseq_markers <- FindAllMarkers(current_seurat, min.pct = 0.25, logfc.threshold = 0.25)
    top1_rna <- current_imput_rnaseq_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1_rna <- top1_rna[order(top1_rna$cluster, decreasing = F),]
    
    clusters <- sort(as.numeric(as.character(unique(current_seurat$seurat_clusters))))
    somePDFPath = paste(output_dir,"05.6SCA_PLOT_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_RNASEQ_GENE_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=7,pointsize=10)
    
    DefaultAssay(current_seurat) <- "peaks"
    
    for(k in 1:length(clusters)){
      print(CoveragePlot(
        object = current_seurat,
        region = unlist(top1$gene[k]), features = unlist(top1_rna[which(top1_rna$cluster == clusters[k]),"gene"]),
        expression.assay = "RNA",
        extend.upstream = 10000,
        extend.downstream = 10000)+
          plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
    }
    dev.off()
    rm(transfer.anchors)
    rm(genes_used)
    rm(current_rnaseq)
  }
  
  # saveRDS(current_seurat, paste(output_dir,"SEURAT_OBJECT_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".RDS", sep = ""))
  
  DefaultAssay(current_seurat) <- "peaks"
  n <- 1000
  topn <- da_peaks %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:n) # %>% filter(avg_log2FC > 1)
  n <- 5000
  topn <- unique(unlist(topn$gene))[1:ifelse(nrow(topn) < n, nrow(topn), n)]
  current_seurat <- subset(current_seurat, features = topn)
  data[[i]] <- current_seurat
  saveRDS(current_seurat, paste(output_dir,"ATAC_SEURAT_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".RDS", sep = ""))
  sample_names <- c(sample_names, pheno_data[i,"SAMPLE_ID"])
  rm(current_seurat)
  rm(closest_genes)
  rm(position_freq)
  rm(top_FC)
  rm(top.da.peak)
}

data_ref <- unique(toupper(pheno_data[,"REF_GENOME"]))

if((length(data) > 1) & (length(data_ref) == 1)){
  
  combined_peaks <- NULL
  for(i in 1:(length(data) -1)){
    combined_peaks <- reduce(x = c(data[[i]]@assays$peaks@ranges, data[[i+1]]@assays$peaks@ranges))
  }
  
  peakwidths <- width(combined_peaks)
  combined_peaks <- combined_peaks[peakwidths  < 10000 & peakwidths > 20]
  
  for(i in 1:length(data)){
    
    current <- FeatureMatrix(
      fragments = data[[i]]@assays$peaks@fragments,
      features = combined_peaks,
      cells = colnames(data[[i]]))
    
    current <- CreateChromatinAssay(current, fragments = data[[i]]@assays$peaks@fragments)
    data[[i]] <- CreateSeuratObject(current, assay = "peaks")
    data[[i]]$Sample <- sample_names[i]
    
  }
  
  rm(current)
  
  data <- merge(x = data[[1]], y = data[2:length(data)],
                add.cell.ids = sample_names)
  
  n <- ifelse((length(data@assays$peaks@var.features)-1) < 50, (length(data@assays$peaks@var.features) - 1), 50)
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 20)
  data <- RunSVD(data, n = n)
  data <- RunUMAP(data, dims = 2:n, reduction = 'lsi')
  data <- RunTSNE(data, reduction = 'lsi', dims = 2:n, check_duplicates = FALSE)
  data <- FindNeighbors(data, dims = 2:n, reduction = 'lsi')
  data <- FindClusters(data, resolution = 0.5, n.iter = 1000)
  data_markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)
  
  p_thresh <- 0.05
  data_markers <- data_markers[data_markers$p_val_adj < p_thresh,]
  data_markers <- data_markers[order(data_markers$p_val_adj, decreasing = F),]
  write.csv(data_markers, paste(output_dir, "05.7SCA_TABLE_INTEGRATED_ATAC_SIG_PEAKS_PADJUST",p_thresh,"_",pheno_data[i,"SAMPLE_ID"],project_name,".csv", sep = ""), quote = F, row.names = T)
  
  plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                      UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                      tSNE_1 = data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                      tSNE_2 = data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                      Cluster = data$Sample)
  
  plotx$Cluster <- factor(plotx$Cluster)
  
  p1_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "All Samples Integrated scATAC UMAP\nBy Samples", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
  p1_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "All Samples Integrated scATAC tSNE\nBy Samples", col = color_conditions$tenx[grep("black",color_conditions$tenx, invert = TRUE)], annot = TRUE, legend_position = "right")
  
  plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                      UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                      tSNE_1 = data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                      tSNE_2 = data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                      Cluster = data$seurat_clusters)
  
  plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.numeric(as.character(plotx$Cluster)))))
  
  p2_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "Cluster", plot_title = "All Samples Integrated scATAC UMAP\nBy Clusters", col = color_conditions$tableau20[grep("black",color_conditions$tableau20, invert = TRUE)], annot = TRUE, legend_position = "right")
  p2_tsne <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", group = "Cluster", plot_title = "All Samples Integrated scATAC tSNE\nBy Clusters", col = color_conditions$tableau20[grep("black",color_conditions$tableau20, invert = TRUE)], annot = TRUE, legend_position = "right")
  
  somePDFPath = paste(output_dir,"05.8SCA_PLOT_UMAP_INTEGRATED_ATAC_SAMPLES_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
  print(p1_umap+p2_umap)
  dev.off()
  
  somePDFPath = paste(output_dir,"05.9SCA_PLOT_tSNE_INTEGRATED_ATAC_SAMPLES_AUTOCLUST_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=6,pointsize=10) 
  print(p1_tsne+p2_tsne)
  dev.off()
  
  wt <- colnames(data_markers)[grep("log.*FC", colnames(data_markers), ignore.case = T)]
  top1 <- data_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
  top1 <- top1[order(top1$cluster, decreasing = F),]
  
  current_plots <- NULL
  
  for(k in 1:length(top1$gene)){
    
    current_plots[[k]] <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),
                                      label = T, pt.size = 0.1, max.cutoff = 'q95')+
      ggtitle(paste("CLUSTER: ", top1$cluster[k], "\n",top1$gene[k],sep = ""))
  }
  
  somePDFPath = paste(output_dir,"06.0SCA_PLOT_INTEGRATED_ATAC_SAMPLES_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=20, height=length(unique(data$seurat_clusters))/4*5,pointsize=10)
  print(do.call("grid.arrange", c(current_plots, ncol=4)))
  dev.off()
  
  n <- 10
  top_n <- data_markers %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
  
  current_clusters <- as.numeric(as.character(sort(unique(data_markers$cluster))))
  somePDFPath = paste(output_dir,"06.1SCA_PLOT_INTEGRATED_ATAC_SAMPLES_TOP_",n,"_MARKERS_IN_EACH_CLUSTER_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=length(unique(data$seurat_clusters))/8*5, height=8,pointsize=10)
  
  for(j in 1:length(current_clusters)){
    current <- top_n[which(top_n$cluster == current_clusters[j]),]
    current <- current[order(current$p_val_adj, decreasing = F),]
    print(VlnPlot(data, features = current$gene, pt.size = 0,
                  # log = TRUE,
                  cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) +
            plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", current_clusters[j], sep = ""),
                            theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"))))      }
  
  dev.off()
  
  if(length(grep("hg19",data_ref, ignore.case = T)) > 0){
    ref_genome <- "hg19"
    load("/Users/luadmpan/Dropbox/KI/Studies/LXXLZ/Website/10X/ATAC/Data/hg19_EnsDb.Hsapiens.v75.RData")
    seqlevelsStyle(ref_annot) <- "UCSC"
    genome(ref_annot) <- "hg19"
  }else if(length(grep("hg38|grch38",data_ref, ignore.case = T)) > 0){
    ref_genome <- "GRCh38.p12"
    load("/Users/luadmpan/Dropbox/KI/Studies/LXXLZ/Website/10X/ATAC/Data/hg38_EnsDb.Hsapiens.v86.RData")
    seqlevelsStyle(ref_annot) <- "UCSC"
    genome(ref_annot) <- "GRCh38.p12"
  }
  
  Annotation(data) <- ref_annot
  rm(ref_annot)
  gene_activities <- GeneActivity(data)
  data[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
  rm(gene_activities)
  
  DefaultAssay(data) <- "ACTIVITY"
  
  data_activity_markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)
  p_thresh <- 0.05
  data_activity_markers <- data_activity_markers[data_activity_markers$p_val_adj < p_thresh,]
  data_activity_markers <- data_activity_markers[order(data_activity_markers$p_val_adj, decreasing = F),]
  top1_activities <- data_activity_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
  top1_activities <- top1_activities[order(top1_activities$cluster, decreasing = F),]
  
  top1 <- data_markers %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
  top1 <- top1[order(top1$cluster, decreasing = F),]
  
  clusters <- unique(sort(as.numeric(as.character(data$seurat_clusters))))
  somePDFPath = paste(output_dir,"06.2SCA_PLOT_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_GENE_ACTIVITY_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=16,pointsize=10)
  
  DefaultAssay(data) <- "peaks"
  
  for(k in 1:length(clusters)){
    print(CoveragePlot(
      object = data,
      region = unlist(top1$gene[k]),
      # features = unlist(top1_activities[which(top1_activities$cluster == clusters[k]),"gene"]),
      extend.upstream = 10000,
      extend.downstream = 10000)+
        plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", clusters[k], sep = ""),
                        theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
  }
  
  dev.off()
  
  somePDFPath = paste(output_dir,"06.3SCA_PLOT_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP1_PEAK_BY_SAMPLES_IN_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=14,pointsize=10)
  
  DefaultAssay(data) <- "peaks"
  Idents(data) <- data$Sample
  
  for(k in 1:length(clusters)){
    print(CoveragePlot(
      object = data,
      region = unlist(top1$gene[k]),
      # features = unlist(top1_activities[which(top1_activities$cluster == clusters[k]),"gene"]),
      extend.upstream = 10000,
      extend.downstream = 10000)+
        plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", clusters[k], sep = ""),
                        theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))
  }
  
  dev.off()
  
}



















