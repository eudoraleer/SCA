library("flowCore")
library("Biobase")
library("flowWorkspace")
library("ggcyto")
library("openCyto")
library("flowStats")

output_dir <- "/proj/lxx_storage/TISSUE_OMICS/GEO/OUTPUT/FLOW/OUTPUT/P144_SDY1389_Bone_Marrow_Flow/"
dir.create(output_dir, recursive = T)
project_name <- gsub(".*\\/(.*)\\/","\\1",output_dir)
system(paste("mkdir ",output_dir, sep = ""))

files <- list.files(pattern=".fcs", ignore.case = T)

data_summary <- NULL

for(k in 1:length(files)){
  x <- read.FCS(files[k])
  x@parameters@data
  current <- data.frame(x@exprs)
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
  
  somePDFPath = paste(output_dir,"SCIENCE_FLOW_STEP1_FSC-A_VS_SSC-A_GATING_",sample_names[k],"_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=8, height=6,pointsize=10)
  for(i in 1:length(gs_data)){
    local_min <- density(fs_data[[i]]@exprs[,"FSC-A"])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,"FSC-A"])$y)))==+2)+1]
    current_ff <- gh_pop_get_data(gs_data[[i]])
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

