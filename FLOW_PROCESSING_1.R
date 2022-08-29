Packages <- c("flowCore","igraph","robustbase","SDMTools","vite","rhandsontable",
              "ggcyto","flowStats","gridExtra","reshape2","ggridges","openCyto",
              "flowViz","patchwork","limma","gplots","plot3D","Biobase",
              "FlowSOM","ConsensusClusterPlus","Rtsne","umap","viridis",
              "flowWorkspace", "openCyto","FLOWMAPR") #CytoExploreR

lapply(Packages, library, character.only = TRUE)

Packages <- c("shiny","shinythemes","shinyFiles",
              "shinydashboard","shinyalert")
lapply(Packages, library, character.only = TRUE)

library("CATALYST")

dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/Bone_Marrow/P146_SDY1389_Thymus_Flow/"
# dir <- "~/Desktop/Desktop/CyTOF_ALL/P146_SDY1389_Thymus_Flow/"
# output_dir <- "/Users/luadmpan/Dropbox/KI/Studies/LXXLZ/Website/CyTOF/CyTOF_OUTPUT/P146_SDY1389_Thymus_Flow/"
output_dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/FLOW/OUTPUT/P146_SDY1389_Thymus_Flow/"
dir.create(output_dir, recursive = T)

# setwd("~/Desktop/Desktop/CyTOF_ALL/P146_SDY1389_Thymus_Flow/")
project_name <- gsub(".*\\/(.*)\\/","\\1",output_dir)
system(paste("mkdir ",output_dir, sep = ""))
x <- read.FCS()
x@parameters@data
x@parameters@data$desc[grep("Beads",x@parameters@data$desc, ignore.case = T)]

############## Convert FCS to csv back to FCS again
x <- read.FCS(files[1], transformation=FALSE)
current_markers <- x@parameters@data$desc
current_markers[which(is.na(current_markers))] <- x@parameters@data$name[which(is.na(current_markers))]
colnames(x@exprs) <- current_markers
x <- x@exprs

# files <- list.files(dir, pattern = "fcs", ignore.case = T, full.names = T)
# data <- lapply(files, function(x){
#   x <- read.FCS(x, transformation=FALSE)
# })
# names(data) <- gsub(".*\\/(.*)","\\1",files)
# data <- as(data,"flowSet")
# sce <- prepData(data)
# table(sce$sample_id)
# names(int_colData(sce))

# # construct SCE
# sce <- prepData(data) # raw_sce <- prepData(raw_data)
# # apply normalization; keep raw data
# res <- normCytof(sce, beads = "dvs", k = 50, 
#                  assays = c("counts", "exprs"), overwrite = FALSE)
# # check number & percentage of bead / removed events
# n <- ncol(sce); ns <- c(ncol(res$beads), ncol(res$removed))
# data.frame(
#   check.names = FALSE, 
#   "#" = c(ns[1], ns[2]), 
#   "%" = 100*c(ns[1]/n, ns[2]/n),
#   row.names = c("beads", "removed"))
# 
# save.image("sce_norm_20210621.RData")

#### CSV to FCS ###########
library("flowCore")
library("Biobase")
library("flowWorkspace")
library("ggcyto")
library("openCyto")
library("flowStats")
files <- list.files(pattern=".csv")
tissue <- NULL
subtissue <- NULL
data <- NULL

for(i in 1:length(files)){
x <- read.csv(files[i], header = T)
tissue <- c(tissue,unique(x$TISSUE))
subtissue <- c(subtissue,unique(x$SUBTISSUE))
colnames(x) <- gsub("\\.|-","_",colnames(x))
x <- x[,grep("TISSUE", colnames(x), ignore.case = T, invert = T)] # bc_separation_dist|mahalanobis_dist|
meta <- data.frame(name=colnames(x),desc=colnames(x))
meta$range <- apply(apply(x,2,range),2,diff)
meta$minRange <- apply(x,2,min)
meta$maxRange <- apply(x,2,max)
x <- new("flowFrame",exprs=as.matrix(x),parameters=AnnotatedDataFrame(meta))
x@parameters@data
data[[i]] <- x
}
fs_data = as(data,"flowSet")
sample_names <- gsub("(.*)\\.(csv|fcs)","\\1",files, ignore.case = T)
fs_data@phenoData@data$name <- cols <- NULL
data_summary <- NULL

for(i in 1:length(files)){
  x <- read.FCS(files[i])
  x@parameters@data
  cols[[i]] <- x@parameters@data
  
  # fs_data <- read.flowSet(path = dir, pattern = ".fcs", alter.names = TRUE, transformation = FALSE)
  fs_data <- as(x, "flowSet")
  fs_data[[1]]@parameters@data
  fs_data@phenoData@data$name
  sample_names <- gsub("^(.*)\\.fcs","\\1",files[i], ignore.case = T)
  sample_names <- gsub("\\s+","_",sample_names)
  fs_data@phenoData@data$name <- sample_names
  fs_data@phenoData@data$name
  
  gs_data <- NULL
  for(i in 1:length(fs_data)){
    gs_data[[i]] <- GatingSet(fs_data[i])
  }
  
  somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_FSC-H_GATING_",sample_names[i],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=8, height=6,pointsize=10)
  for(i in 1:length(gs_data)){
    current_ff <- gh_pop_get_data(gs_data[[i]])
    chnl <- colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
    mylimits <- ggcyto_par_set(limits = "instrument")
    if(is.null(range_fsc_ssc)){
      g <- openCyto:::.flowClust.2d(current_ff, channels = chnl, K=3, target=c(1e5,1e4), quantile=0.95, filterId = "SelectedCells")
    }else{
      g <- openCyto:::.flowClust.2d(current_ff, channels = chnl, K=3, target=c(as.numeric(as.character(range_fsc_ssc))), quantile=0.95, filterId = "SelectedCells")
    }
    # print(autoplot(current_ff, x = chnl[1], y = chnl[2], bin = 300) + mylimits + geom_gate(g) + geom_stats(adjust = 1))
    gs_pop_add(gs_data[[i]], g)
    recompute(gs_data[[i]])
    print(ggcyto(gs_data[[1]],aes(x=chnl[1],y=chnl[2]),subset="SelectedCells")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
    
  }
  dev.off()
  
  somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_SSC-A_GATING_",sample_names[i],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=8, height=6,pointsize=10)
  for(i in 1:length(gs_data)){
    local_min <- density(fs_data[[i]]@exprs[,"FSC-A"])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,"FSC-A"])$y)))==+2)+1]
    current_ff <- gh_pop_get_data(gs_data[[i]], "defaultEllipsoidGate")
    chnl <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
    # g <- polygonGate(filterId = "Singlets","FSC-A"=c(5e4,24e4,24e4,8e4),"SSC-A"=c(0,2e4,15e4,5e4)) # define gate
    # g <- gate_singlet(current_ff, area = chnl[1], height = chnl[2], filterId = "Singlets")
    g <- openCyto:::.mindensity(current_ff, channels = chnl, filterId = "Singlets", gate_range=c(local_min[1]-500,local_min[2]-1000))
    mylimits <- ggcyto_par_set(limits = "instrument")
    gs_pop_add(gs_data[[i]], g, parent = "defaultEllipsoidGate")
    recompute(gs_data[[i]])
    print(ggcyto(gs_data[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
    
    # print(autoplot(current_ff, chnl, strip.text = "gate") + mylimits + geom_density(fill="forestgreen") + geom_gate(g))
    # print(autoplot(fs_data[[i]], chnl, colnames(fs_data[[i]])[grep("SSC.*A", colnames(fs_data[[i]]), ignore.case = T)], bin = 200, strip.text = "gate") + geom_gate(g) + geom_stats(adjust = 0.5))
  }
  dev.off()
  
  dir.create(paste(output_dir, "Filtered_FCS_Input/", sep = ""), recursive = T)
  filtered_data <- NULL
  for(i in 1:length(gs_data)){
    ps <- data.frame(gs_pop_get_count_fast(gs_data[[i]]))
    ps$Percent_of_parent <- ps$Count/ps$ParentCount
    data_summary <- rbind(data_summary, ps)
    filtered_data[[i]] <- gh_pop_get_data(gs_data[[i]], y = "Singlets")
    write.FCS(filtered_data[[i]], paste(output_dir, "Filtered_FCS_Input/Filtered_",sample_names[i],".fcs", sep = ""))
  }
  
  print(paste("Done: ", files[i]))
}
cols <- NULL
data_summary <- NULL

for(i in 1:length(files)){
  x <- read.FCS(files[i])
  x@parameters@data
  cols[[i]] <- x@parameters@data
  
  # fs_data <- read.flowSet(path = dir, pattern = ".fcs", alter.names = TRUE, transformation = FALSE)
  fs_data <- as(x, "flowSet")
  fs_data[[1]]@parameters@data
  fs_data@phenoData@data$name
  sample_names <- gsub("^(.*)\\.fcs","\\1",files[i], ignore.case = T)
  sample_names <- gsub("\\s+","_",sample_names)
  fs_data@phenoData@data$name <- sample_names
  fs_data@phenoData@data$name
  
  gs_data <- NULL
  for(i in 1:length(fs_data)){
    gs_data[[i]] <- GatingSet(fs_data[i])
  }
  
  somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_FSC-H_GATING_",sample_names[i],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=8, height=6,pointsize=10)
  for(i in 1:length(gs_data)){
    current_ff <- gh_pop_get_data(gs_data[[i]])
    chnl <- colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
    mylimits <- ggcyto_par_set(limits = "instrument")
    if(is.null(range_fsc_ssc)){
      g <- openCyto:::.flowClust.2d(current_ff, channels = chnl, K=3, target=c(1e5,1e4), quantile=0.95, filterId = "SelectedCells")
    }else{
      g <- openCyto:::.flowClust.2d(current_ff, channels = chnl, K=3, target=c(as.numeric(as.character(range_fsc_ssc))), quantile=0.95, filterId = "SelectedCells")
    }
    # print(autoplot(current_ff, x = chnl[1], y = chnl[2], bin = 300) + mylimits + geom_gate(g) + geom_stats(adjust = 1))
    gs_pop_add(gs_data[[i]], g)
    recompute(gs_data[[i]])
    print(ggcyto(gs_data[[1]],aes(x=chnl[1],y=chnl[2]),subset="SelectedCells")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
    
  }
  dev.off()
  
  somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_SSC-A_GATING_",sample_names[i],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=8, height=6,pointsize=10)
  for(i in 1:length(gs_data)){
    local_min <- density(fs_data[[i]]@exprs[,"FSC-A"])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,"FSC-A"])$y)))==+2)+1]
    current_ff <- gh_pop_get_data(gs_data[[i]], "defaultEllipsoidGate")
    chnl <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
    # g <- polygonGate(filterId = "Singlets","FSC-A"=c(5e4,24e4,24e4,8e4),"SSC-A"=c(0,2e4,15e4,5e4)) # define gate
    # g <- gate_singlet(current_ff, area = chnl[1], height = chnl[2], filterId = "Singlets")
    g <- openCyto:::.mindensity(current_ff, channels = chnl, filterId = "Singlets", gate_range=c(local_min[1]-500,local_min[2]-1000))
    mylimits <- ggcyto_par_set(limits = "instrument")
    gs_pop_add(gs_data[[i]], g, parent = "defaultEllipsoidGate")
    recompute(gs_data[[i]])
    print(ggcyto(gs_data[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
    
    # print(autoplot(current_ff, chnl, strip.text = "gate") + mylimits + geom_density(fill="forestgreen") + geom_gate(g))
    # print(autoplot(fs_data[[i]], chnl, colnames(fs_data[[i]])[grep("SSC.*A", colnames(fs_data[[i]]), ignore.case = T)], bin = 200, strip.text = "gate") + geom_gate(g) + geom_stats(adjust = 0.5))
  }
  dev.off()
  
  dir.create(paste(output_dir, "Filtered_FCS_Input/", sep = ""), recursive = T)
  filtered_data <- NULL
  for(i in 1:length(gs_data)){
    ps <- data.frame(gs_pop_get_count_fast(gs_data[[i]]))
    ps$Percent_of_parent <- ps$Count/ps$ParentCount
    data_summary <- rbind(data_summary, ps)
    filtered_data[[i]] <- gh_pop_get_data(gs_data[[i]], y = "Singlets")
    write.FCS(filtered_data[[i]], paste(output_dir, "Filtered_FCS_Input/Filtered_",sample_names[i],".fcs", sep = ""))
  }
  
  print(paste("Done: ", files[i]))
}

fs_data@phenoData@data$name
############################
files <- list.files(pattern=".fcs", ignore.case = T)

cols <- NULL
data_summary <- NULL

for(k in 1:length(files)){
  x <- read.FCS(files[k])
  x@parameters@data
  current <- data.frame(x@exprs)
  cols[[i]] <- x@parameters@data

# fs_data <- read.flowSet(path = dir, pattern = ".fcs", alter.names = TRUE, transformation = FALSE)
fs_data <- as(x, "flowSet")
fs_data[[1]]@parameters@data
fs_data@phenoData@data$name
sample_names <- gsub("^(.*)\\.fcs","\\1",files[k], ignore.case = T)
sample_names <- gsub("\\s+","_",sample_names)
fs_data@phenoData@data$name <- sample_names
fs_data@phenoData@data$name

gs_data <- NULL
for(i in 1:length(fs_data)){
  gs_data[[i]] <- GatingSet(fs_data[i])
}

somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_FSC-H_GATING_",sample_names[k],"_",project_name,".pdf", sep = "")
pdf(file=somePDFPath, width=8, height=6,pointsize=10)
for(i in 1:length(gs_data)){
  current_ff <- gh_pop_get_data(gs_data[[i]])
  chnl <- colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
  mylimits <- ggcyto_par_set(limits = "instrument")
    g <- openCyto:::.flowClust.2d(current_ff, channels = chnl, K=3, target=c(1e5,1e4), quantile=0.95, filterId = "SelectedCells")
  # print(autoplot(current_ff, x = chnl[1], y = chnl[2], bin = 300) + mylimits + geom_gate(g) + geom_stats(adjust = 1))
  gs_pop_add(gs_data[[i]], g)
  recompute(gs_data[[i]])
  print(ggcyto(gs_data[[i]],aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
  
}
dev.off()

somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_SSC-A_GATING_",sample_names[k],"_",project_name,".pdf", sep = "")
pdf(file=somePDFPath, width=8, height=6,pointsize=10)
for(i in 1:length(gs_data)){
  local_min <- density(fs_data[[i]]@exprs[,"FSC-A"])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,"FSC-A"])$y)))==+2)+1]
  current_ff <- gh_pop_get_data(gs_data[[i]], "defaultEllipsoidGate")
  chnl <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
  # g <- polygonGate(filterId = "Singlets","FSC-A"=c(5e4,24e4,24e4,8e4),"SSC-A"=c(0,2e4,15e4,5e4)) # define gate
  # g <- gate_singlet(current_ff, area = chnl[1], height = chnl[2], filterId = "Singlets")
  g <- openCyto:::.mindensity(current_ff, channels = chnl, filterId = "Singlets", gate_range=c(local_min[1]-500,local_min[2]-1000))
  mylimits <- ggcyto_par_set(limits = "instrument")
  gs_pop_add(gs_data[[i]], g, parent = "defaultEllipsoidGate")
  recompute(gs_data[[i]])
  print(ggcyto(gs_data[[i]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g)+ggcyto_par_set(limits = "instrument")) # check gate
  
  # print(autoplot(current_ff, chnl, strip.text = "gate") + mylimits + geom_density(fill="forestgreen") + geom_gate(g))
  # print(autoplot(fs_data[[i]], chnl, colnames(fs_data[[i]])[grep("SSC.*A", colnames(fs_data[[i]]), ignore.case = T)], bin = 200, strip.text = "gate") + geom_gate(g) + geom_stats(adjust = 0.5))
}
dev.off()

dir.create(paste(output_dir, "Filtered_FCS_Input/", sep = ""), recursive = T)
filtered_data <- NULL
for(i in 1:length(gs_data)){
  ps <- data.frame(gs_pop_get_count_fast(gs_data[[i]]))
  ps$Percent_of_parent <- ps$Count/ps$ParentCount
  data_summary <- rbind(data_summary, ps)
  filtered_data[[i]] <- gh_pop_get_data(gs_data[[i]], y = "Singlets")
  write.FCS(filtered_data[[i]], paste(output_dir, "Filtered_FCS_Input/Filtered_",sample_names[k],".fcs", sep = ""))
}

print(paste("Done: ", files[i]))
}

write.csv(data_summary, paste(output_dir, "CyTOF_TABLE_GATING_SUMMARY_",project_name,".csv",sep = ""),row.names = F)

cofactor <- 5
filtered_data <- as(filtered_data,"flowSet")
sampleNames(filtered_data) <- sampleNames(fs_data)
pData(filtered_data)$name <- sampleNames(fs_data)
data <- NULL
for(i in 1:length(filtered_data)){
  # data <- rbind(data,data.frame(SAMPLE_ID = pData(filtered_data)$name[i], asinh(exprs(filtered_data[[i]])/cofactor)))
  data <- rbind(data,data.frame(SAMPLE_ID = pData(filtered_data)$name[i], exprs(filtered_data[[i]])))
}

fs_data[[1]]@parameters@data
colnames(data) <- c("SAMPLE_ID",fs_data[[1]]@parameters@data$desc)
data <- data[,grep("P74S|P74B|P76S|P76B|^[0-9]+[A-Z]+$|Beads|DNA|Viability|^NA$|^NA\\.[0-9]+$", colnames(data), ignore.case=T, invert=T)]
data <- data[,which(!is.na(colnames(data)))]
data <- data[,grep("P74S|P74B|P76S|P76B|^[0-9]+[A-Z]+$|Beads|DNA|Viability|^NA$|^NA\\.[0-9]+$", colnames(data), ignore.case=T, invert=T)]
colnames(data)
data$SAMPLE_ID <- gsub("^_(.*)\\.[0-9]+\\.fcs","\\1",data$SAMPLE_ID, ignore.case = T)
data$TISSUE <- tissue
data$SUBTISSUE <- subtissue
summary(data[,-1])
write.csv(data, paste(output_dir, "CyTOF_POST_GATING_",project_name,".csv",sep = ""),row.names = F)





