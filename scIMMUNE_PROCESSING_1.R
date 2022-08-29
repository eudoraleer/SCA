###########################################################################################################################
# SCIENCE FIGURES: scIMMUNE
###########################################################################################################################

# Remove P214 Spleen
setwd("~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/")
project_name <- "scImmune_All_Tissues"
input_dir <- "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/"
phenodata_dir <- "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/RAW_TCRBCR_ALL_TISSUES.csv"
output_dir = "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/"
dir.create(output_dir, recursive = T)

Packages <- c("Seurat","monocle","SingleCellExperiment","ggthemes","plot3D","ggalluvial","dplyr","alluvial", "tcR","scRepertoire","immunarch","powerTCR","circlize","scales")
suppressMessages(lapply(Packages, library, character.only = TRUE))

source('~/Dropbox/KI/Studies/LXXLZ/Website/Web_Scripts/SCA_Analysis_Functions_V1.0.0.R')
# source('/mnt/SCA/PIPELINES/SCA_Analysis_Functions_V1.0.0.R')
# source('/proj/store_eudoraleer/PROJECTS/scRNASEQ/SCRIPTS/SCA_Analysis_Functions_V1.0.0.R')

color_conditions <- color_ini()

#####################################################################################################
contig_monarch_bak <- contig_monarch

# head(contig_monarch)
x$SIDCELLID <- paste(x$SID, x$CELL_ID, sep = "_")
x$SIDCELLID <- gsub("-|\\.","_",x$SIDCELLID)

for(i in 1:length(contig_monarch$data)){
  c <- contig_monarch$data[[i]]
  c$SID <- meta[match(names(contig_monarch$data)[i], meta$Sample),"SID"]
    c$SIDCELLID <- paste(c$SID, gsub("(.*?);.*","\\1",c$Barcode), sep = "_")
    c$SIDCELLID <- gsub("-|\\.","_",c$SIDCELLID)
    c$SUBTISSUE <- meta[match(names(contig_monarch$data)[i], meta$Sample),"SUBTISSUE"]
    c$TISSUE <- meta[match(names(contig_monarch$data)[i], meta$Sample),"TISSUE"]
    c$CELL_SUBTYPE <- x@meta.data[match(c$SIDCELLID, x$SIDCELLID),"CELL_SUBTYPE"]
    c$CELL_SUBTYPE[which(is.na(c$CELL_SUBTYPE))] <- ifelse(x@meta.data[match(names(contig_monarch$data)[i], x$Sample),"CELL_SUBTYPE"] == "BCR","B-CELLS","T-CELLS")
    contig_monarch$data[[i]] <- c
    
}

current <- meta_explore[meta_explore$CELL_TYPE == "TCR",]
current_meta <- x@meta.data[which(x@meta.data$Sample %in% current$Sample),]
current$SUBTISSUE <- paste(x@meta.data[match(current$Sample, x@meta.data$Sample),"TISSUE"], ":", current$SUBTISSUE, sep = "")
current_meta$SUBTISSUE <- paste(x@meta.data[match(current_meta$Sample, x@meta.data$Sample),"TISSUE"], ":", current_meta$SUBTISSUE, sep = "")

################################################################################################################
############################ By BCR / TCR separately ###########################################################
for(i in 1:length(contig_table)){
  
  current_type <- contig_types[i]
  current_single <- lengthContig(contig_table[[i]], cloneCall="aa", chains = "single", exportTable = TRUE,group = "CELL_SUBTYPE")
  write.csv(current_single, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_AA_SEPARATE_BY_CELL_SUBTYPES_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
  
  current_combined <- lengthContig(contig_table[[i]], cloneCall="aa", chains = "combined", exportTable = TRUE,group = "CELL_SUBTYPE")
  write.csv(current_combined, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_AA_COMBINED_BY_CELL_SUBTYPES_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
  
  current_single <- lengthContig(contig_table[[i]], cloneCall="nt", chains = "single", exportTable = TRUE,group = "CELL_SUBTYPE")
  write.csv(current_single, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_NUCLEOTIDES_SEPARATE_BY_CELL_SUBTYPES_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
  
  current_combined <- lengthContig(contig_table[[i]], cloneCall="nt", chains = "combined", exportTable = TRUE,group = "CELL_SUBTYPE")
  write.csv(current_combined, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_NUCLEOTIDES_COMBINED_BY_CELL_SUBTYPES_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
}

######### SUBTISSUES CLONOTYPE FREQUENCY #########################################################

clonotypes_freq <- NULL
current_files <- list.files(path = output_dir, pattern = "^04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_.*BY_CELL_SUBTYPES_.*csv", full.names = T)
tissues_kmer_summary <- NULL

for(i in 1:length(current_files)){
  current_celltype <- NULL
  current_datatype <- NULL
  combined <- NULL
  combined <- ifelse(gsub(".*(COMBINED).*","\\1",current_files[i]) == "COMBINED", T, F)
  current <- readfile(current_files[i])
  if(length(grep("BCR", current_files[i])) > 0){
    current_celltype <- "BCR"
  }else{
    current_celltype <- "TCR"
  }
  
  if(combined == T){
    colnames(current) <- c("LENGTH","SEQUENCE","CELL_SUBTYPE","SUBTISSUE")
  }else{
    colnames(current) <- c("LENGTH","SEQUENCE","SUBTISSUE","CHAIN","CELL_SUBTYPE")
  }
  
  current_datatype <- gsub(".*TABLE_CLONOTYPE_ABUNDANCE_(.*?)_.*","\\1",current_files[i])
  current$LENGTH <- nchar(gsub("_","",current$SEQUENCE))
  # temp <- data.frame(table(current[,c("LENGTH","SUBTISSUE")]))
  # temp <- temp[!temp$Freq == 0,]
  current <- split(current, current$SUBTISSUE)
  current <- lapply(current, function(x){x <- split(x,x$CELL_SUBTYPE)})
  for(k in 1:length(current)){
    c <- current[[k]]
  for(j in 1:length(c)){
    temp <- c[[j]]
    current_median <- NULL
    if(combined == T){
      current_median <- as.numeric(as.character(summary(temp$LENGTH)[grep("Median",names(summary(temp$LENGTH)))]))
      current_frequency <- table(temp$LENGTH)/nrow(temp)
      temp$FREQUENCY <- as.numeric(as.character(current_frequency[match(temp$LENGTH,names(current_frequency))]))
    }else{
      temp <- split(temp, temp$CHAIN)
      temp <- lapply(temp, function(x){
        current_frequency <- NULL
        current_frequency <- table(x$LENGTH)/nrow(x)
        x <- cbind(x, MEDIAN_LENGTH = as.numeric(summary(x$LENGTH)[grep("Median",names(summary(x$LENGTH)))]),
                   FREQUENCY = as.numeric(as.character(current_frequency[match(x$LENGTH,names(current_frequency))])))})
      temp <- do.call(rbind.data.frame, temp)
    }
    
    if(combined == T){
      temp2 <- data.frame(
        SUBTISSUE = unique(temp$SUBTISSUE),
        CELL_SUBTYPE = unique(temp$CELL_SUBTYPE),
        DATA_TYPE = paste(current_celltype, "_", current_datatype, sep = ""),
        CHAIN = "COMBINE",
        STATUS = ifelse(combined == T, "COMBINED","SINGLE"),
        MEDIAN_LENGTH = current_median)
    }else{
      temp2 <- cbind(
        SUBTISSUE = unique(temp$SUBTISSUE),
        CELL_SUBTYPE = unique(temp$CELL_SUBTYPE),
        DATA_TYPE = paste(current_celltype, "_", current_datatype, sep = ""),
        STATUS = ifelse(combined == T, "COMBINED","SINGLE"),
        data.frame(unique(temp[,c("CHAIN","MEDIAN_LENGTH")])
        ))
    }
    tissues_kmer_summary <- rbind(tissues_kmer_summary,temp2)
    c[[j]] <- temp
  }
  current[[k]] <- do.call(rbind.data.frame, c)
  }
  current <- do.call(rbind.data.frame, current)
  
  dir.create(paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/", sep = ""))
  write.csv(current, paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/WITH_FREQUENCY_INFO_CELL_SUBTYPES_",gsub(".*\\/(.*)","\\1",current_files[i]), sep = ""), quote=F, row.names = F)
  write.csv(tissues_kmer_summary, paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/SUBTISSUES_KMER_MEDIAN_LENGTH_SUMMARY_CELL_SUBTYPES_",gsub(".*\\/(.*)","\\1",current_files[i]), sep = ""), quote=F, row.names = F)
  
  # current$SUBTISSUE <- factor(current$SUBTISSUE, levels = sort(unique(current$SUBTISSUE)))
  tissues <- sort(unique(current$SUBTISSUE))
  for(l in 1:length(tissues)){
  current_file_name <- gsub(".*\\/(.*).csv",paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/",tissues[l],"_\\1.png",sep = ""),current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=4000, units = "px", res = 300)
  p1 <- ggplot(current[which(current$SUBTISSUE == tissues[l]),],aes(LENGTH,FREQUENCY,fill=CELL_SUBTYPE))+ 
    geom_bar(stat = "identity", position = "dodge")+ #, color = "black", size = 1
    theme_classic()+
    # facet_wrap(~CELL_SUBTYPE)+
    # scale_x_continuous(breaks = seq(0, max(current$LENGTH)+2, by = 5))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 0, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          # legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          strip.text.x = element_text(size = 25),
          strip.text.y = element_text(size = 25),
          plot.title = element_text(size = 20, face=1, hjust = 0.5))+
    scale_fill_manual(values = gen_colors(color_conditions$bright, length(unique(current$SUBTISSUE))))+
    guides(color=guide_legend(title="SUBTISSUES", ncol = 3),fill = guide_legend(title="SUBTISSUES", ncol = 3))
  if(combined == F){
    p1 <- p1 + facet_wrap(~CHAIN, ncol=1)
  }
  print(p1+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), tissues[l]," ",current_celltype, ":CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
  # current$SUBTISSUE <- factor(current$SUBTISSUE, levels = rev(sort(unique(current$SUBTISSUE))))
  current_file_name <- gsub(".*\\/(.*).csv",paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/",tissues[l],"_\\1_BOXPLOT.png",sep = ""),current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=5000, units = "px", res = 300)
  p2 <- ggplot(current[which(current$SUBTISSUE == tissues[l]),],aes(LENGTH,CELL_SUBTYPE,fill=CELL_SUBTYPE))+ geom_boxplot()+ theme_classic()+
    # scale_x_continuous(breaks = seq(0, max(current$LENGTH)+2, by = 5))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 0, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          # legend.key.size = unit(1, "cm"),
          legend.position = "none",
          strip.text.x = element_text(size = 25),
          strip.text.y = element_text(size = 25),
          plot.title = element_text(size = 20, face=1, hjust = 0.5))+
    scale_fill_manual(values = gen_colors(color_conditions$bright, length(unique(current$SUBTISSUE))))+
    guides(color=guide_legend(title="SUBTISSUES", ncol = 3),fill = guide_legend(title="SUBTISSUES", ncol = 3))
  if(combined == F){
    p2 <- p2 + facet_wrap(~CHAIN, ncol=1)
  }
  print(p2+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), tissues[l]," ",current_celltype, ":CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
  library(ggridges)
  # current$SUBTISSUE <- factor(current$SUBTISSUE, levels = sort(unique(current$SUBTISSUE)))
  current_file_name <- gsub(".*\\/(.*).csv",paste(output_dir,"/05.0CLONOTYPE_ABUNDANCE_CELL_SUBTYPES/",tissues[l],"_\\1_RIDGE.png",sep = ""),current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=5000, units = "px", res = 300)
  p2 <- ggplot(current[which(current$SUBTISSUE == tissues[l]),],aes(LENGTH,CELL_SUBTYPE,fill=CELL_SUBTYPE))+ geom_density_ridges()+ theme_classic()+
    # scale_x_continuous(breaks = seq(0, max(current$LENGTH)+2, by = 5))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 0, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          # legend.key.size = unit(1, "cm"),
          legend.position = "none",
          strip.text.x = element_text(size = 25),
          strip.text.y = element_text(size = 25),
          plot.title = element_text(size = 20, face=1, hjust = 0.5))+
    scale_fill_manual(values = gen_colors(color_conditions$bright, length(unique(current$SUBTISSUE))))+
    guides(color=guide_legend(title="SUBTISSUES", ncol = 3),fill = guide_legend(title="SUBTISSUES", ncol = 3))
  if(combined == F){
    p2 <- p2 + facet_wrap(~CHAIN, ncol=1)
  }
  print(p2+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), tissues[l]," ",current_celltype, ":CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  }
}


#### ALL WRONG!!!!! scRepertoire Package Bug
# abundance_contig <- NULL
# abundance_names <- NULL
# clonecall_types
# k <- 1
# for(i in 1:length(clonecall_types)){
#   for(j in 1:length(cell_types)){
#     abundance_contig[[k]] <- abundanceContig(contig_table[[j]], cloneCall = clonecall_types[i], group = c("CELL_SUBTYPE","TISSUE"),exportTable = T)
#     abundance_contig[[k]] <- abundance_contig[[k]][order(abundance_contig[[k]]$Abundance, decreasing = T),]
#     write.csv(abundance_contig[[k]], paste(output_dir,"scIMMUNE_CLONOTYPES_ABUNDANCE_TISSUE_",cell_types[j],"_",toupper(clonecall_types[i]),"_",current_type,".csv", sep = ""), quote = F, row.names = F)
#     abundance_names <- c(abundance_names, paste(cell_types[j],"_",toupper(clonecall_types[i]),"_",current_type, sep = ""))
#     k <- k+1
#   }
# }

abundance_contig <- NULL
abundance_names <- NULL
clonecall_cols <- c("CTgene","CTnt","CTstrict","CTaa")
clonecall_cols
clonecall_types

m <- 1
for(i in 1:length(clonecall_types)){
  for(j in 1:length(cell_types)){
    current <- contig_table[[j]]
    for(k in 1:length(current)){
      current[[k]] <- data.frame(table(current[[k]][,c(clonecall_cols[i],"SUBTISSUE","CELL_SUBTYPE")]))
    }
    current <- do.call(rbind.data.frame, current)
    abundance_contig[[m]] <- current
    abundance_names <- c(abundance_names, paste(cell_types[j],"_",toupper(clonecall_types[i]),sep = ""))
    write.csv(abundance_contig[[m]], paste(output_dir,"scIMMUNE_CLONOTYPES_ABUNDANCE_TISSUE_",cell_types[j],"_",toupper(clonecall_types[i]),".csv", sep = ""), quote = F, row.names = F)
    
    m <- m+1
    
  }
}
saveRDS(abundance_contig,"~/Desktop/abundance_contig_5.2SCA_PLOT_ALLUVIAL.RDS")

dir.create("~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/5.2SCA_PLOT_ALLUVIAL/")
abundance_contig <- readRDS("~/Desktop/abundance_contig_5.2SCA_PLOT_ALLUVIAL.RDS")
n <- 10
for(k in 1:length(abundance_contig)){
  current <- split(abundance_contig[[k]], abundance_contig[[k]]$SUBTISSUE)
  for(j in 1:length(current)){
    temp <- current[[j]]
    temp <- temp[temp[grep("freq", colnames(temp), ignore.case = T)] > 0,]
    temp$PROPORTION <- temp[,grep("freq", colnames(temp), ignore.case = T)]/sum(temp[,grep("freq", colnames(temp), ignore.case = T)])
    temp <- split(temp, temp$CELL_SUBTYPE)
    temp <- temp[unlist(lapply(temp, nrow)) > 0]
    temp <- lapply(temp, function(y){
      y <- y[order(y$PROPORTION, decreasing = T),]
      y <- y[1:ifelse(nrow(y)<n, nrow(y), n),]
    })
    temp <- do.call(rbind.data.frame, temp)
    current[[j]] <- temp
    temp <- NULL
  }
  current <- do.call(rbind.data.frame,current)
  abundance_contig[[k]] <- current
  current <- split(current, current$CELL_SUBTYPE)
  current <- current[unlist(lapply(current, nrow)) > 0]
  for(l in 1:length(current)){
    plotx <- current[[l]]
    cname <- colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)]
    cct <- as.character(unique(plotx$CELL_SUBTYPE))
    colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] <- "CLONE_TYPE"
    library(reshape2)
    plotx <- dcast(plotx, CLONE_TYPE~SUBTISSUE, value.var = "PROPORTION")
    plotx[is.na(plotx)] <- 0
    write.table(plotx,paste(output_dir,"5.2SCA_PLOT_ALLUVIAL/5.2SCA_PLOT_ALLUVIAL_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    class(plotx) <- c("immunr_dynamics","data.table","data.frame")
    p1 <- vis(plotx)+ggtitle(paste(cct, ": ",gsub("(.*?)_(.*)","\\1",abundance_names[k]), " TOP ",n," CLONOTYPES (",cname,")", sep = ""))+
      scale_fill_manual(values = gen_colors(color_conditions$colorful, nrow(plotx)))+
      theme(plot.margin = unit(c(2,2,2,2), "cm"),
            axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "cm"),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 20, face=2),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "none")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+xlab("SUB-TISSUES")+ylab("PROPORTION")
    somePNGPath = paste(output_dir,"5.2SCA_PLOT_ALLUVIAL/5.2SCA_PLOT_ALLUVIAL_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".png", sep = "")
    png(somePNGPath, width=5000, height=3000, units = "px", res = 300)
    print(p1)
    dev.off()
  }
}

saveRDS(abundance_contig, "~/Desktop/abundance_contig_subcelltypes_5.2SCA_PLOT_ALLUVIAL.RDS")

abundance_contig <- readRDS("~/Desktop/abundance_contig_5.2SCA_PLOT_ALLUVIAL.RDS")
clonecall_cols

n <- 10
for(k in 1:length(abundance_contig)){
  current <- abundance_contig[[k]]
  current$CELL_SUBTYPE <- gsub(".*B-CELL.*","B-CELL",current$CELL_SUBTYPE, ignore.case = T)
  current$CELL_SUBTYPE <- gsub(".*T-CELL.*","T-CELL",current$CELL_SUBTYPE, ignore.case = T)
  current <- split(current, current$SUBTISSUE)
  for(j in 1:length(current)){
    temp <- current[[j]]
    temp <- temp[temp[grep("freq", colnames(temp), ignore.case = T)] > 0,]
    temp$PROPORTION <- temp[,grep("freq", colnames(temp), ignore.case = T)]/sum(temp[,grep("freq", colnames(temp), ignore.case = T)])
    temp <- split(temp, temp$CELL_SUBTYPE)
    temp <- temp[unlist(lapply(temp, nrow)) > 0]
    temp <- lapply(temp, function(y){
      y <- y[order(y$PROPORTION, decreasing = T),]
      y <- y[1:ifelse(nrow(y)<n, nrow(y), n),]
    })
    temp <- do.call(rbind.data.frame, temp)
    current_clones <- as.character(unique(temp[,colnames(temp)[grep("SUBTISSUE|CELL_SUBTYPE|Freq|PROPORTION",colnames(temp), ignore.case = T, invert = T)]]))
    temp2 <- NULL
    for(i in 1:length(current_clones)){
      currentcurrent <- temp[which(temp[,colnames(temp)[grep("SUBTISSUE|CELL_SUBTYPE|Freq|PROPORTION",colnames(temp), ignore.case = T, invert = T)]] == current_clones[i]),]
      current_name <- colnames(currentcurrent)[grep("SUBTISSUE|CELL_SUBTYPE|Freq|PROPORTION",colnames(currentcurrent), ignore.case = T, invert = T)]
      currentcurrent <- data.frame(unique(currentcurrent[which(currentcurrent[,current_name] == current_clones[i]),c(current_name,"SUBTISSUE","CELL_SUBTYPE")]), Freq = sum(currentcurrent$Freq), PROPORTION = sum(currentcurrent$PROPORTION))
      temp2 <- rbind(temp2,currentcurrent) 
    }
    temp <- temp2
    rm(temp2)
    current[[j]] <- temp
    temp <- NULL
  }
  current <- do.call(rbind.data.frame,current)
  abundance_contig[[k]] <- current
  current <- split(current, current$CELL_SUBTYPE)
  current <- current[unlist(lapply(current, nrow)) > 0]
  for(l in 1:length(current)){
    plotx <- current[[l]]
    cname <- colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)]
    cct <- as.character(unique(plotx$CELL_SUBTYPE))
    colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] <- "CLONE_TYPE"
    library(reshape2)
    plotx <- dcast(plotx, CLONE_TYPE~SUBTISSUE, value.var = "PROPORTION")
    plotx[is.na(plotx)] <- 0
    write.table(plotx,paste(output_dir,"5.2SCA_PLOT_ALLUVIAL/5.2SCA_PLOT_ALLUVIAL_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    class(plotx) <- c("immunr_dynamics","data.table","data.frame")
    p1 <- vis(plotx)+ggtitle(paste(cct, ": ",gsub("(.*?)_(.*)","\\1",abundance_names[k]), " TOP ",n," CLONOTYPES (",cname,")", sep = ""))+
      scale_fill_manual(values = gen_colors(color_conditions$colorful, nrow(plotx)))+
      theme(plot.margin = unit(c(2,2,2,2), "cm"),
            axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "cm"),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 20, face=2),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "none")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+xlab("SUB-TISSUES")+ylab("PROPORTION")
    somePNGPath = paste(output_dir,"5.2SCA_PLOT_ALLUVIAL/5.2SCA_PLOT_ALLUVIAL_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".png", sep = "")
    png(somePNGPath, width=5000, height=3000, units = "px", res = 300)
    print(p1)
    dev.off()
  }
  
  
}

saveRDS(abundance_contig, "~/Desktop/abundance_contig_TBcells_5.2SCA_PLOT_ALLUVIAL.RDS")

# for(j in 1:length(current)){
  #   temp <- split(current[[j]], current[[j]]$CELL_SUBTYPE)
  #   temp <- temp[unlist(lapply(temp, nrow)) > 0]
  #   temp <- lapply(temp, function(y){
  #     y <- data.frame(y, PROPORTION = y[,grep("freq", colnames(y), ignore.case = T)]/sum(y[,grep("freq", colnames(y), ignore.case = T)]))
  #   })
  #   temp <- do.call(rbind.data.frame, temp)
  # }
#   temp <- temp[temp$PROPORTION > 0,]
# }

library(reshape2)
dir.create("~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE")
abundance_contig <- readRDS("~/Desktop/abundance_contig_subcelltypes_5.2SCA_PLOT_ALLUVIAL.RDS")

n <- 10
for(k in 1:length(abundance_contig)){
  current <- abundance_contig[[k]]
  cell_types <- as.character(sort(unique(current$CELL_SUBTYPE)))
  for(j in 1:length(cell_types)){
    plotx <- current[current$CELL_SUBTYPE == cell_types[j],]
    plotx <- plotx[order(plotx$Freq, decreasing = T),]
    ctop <- unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])[1:ifelse(length(unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])) < n, length(unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])), n)]
    plotx <- plotx[which(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] %in% ctop),]
    cname <- colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)]
    cct <- as.character(unique(plotx$CELL_SUBTYPE))
    colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] <- "CLONE_TYPE"
    plotx <- dcast(plotx, CLONE_TYPE~SUBTISSUE, value.var = "PROPORTION")
    plotx[is.na(plotx)] <- 0
    plotx$CLONE_TYPE <- gsub("^NA;|;NA$","",plotx$CLONE_TYPE)
    plotx$CLONE_TYPE <- gsub(";NA;",";",plotx$CLONE_TYPE)
    plotx$CLONE_TYPE <- gsub("^None;|;None$","",plotx$CLONE_TYPE, ignore.case = T)
    plotx$CLONE_TYPE <- gsub(";NONE;",";",plotx$CLONE_TYPE, ignore.case = T)
    
    write.table(plotx,paste(output_dir,"5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE/5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    class(plotx) <- c("immunr_dynamics","data.table","data.frame")
    p1 <- vis(plotx)+ggtitle(paste(cct, " ONLY: ",gsub("(.*?)_(.*)","\\1",abundance_names[k]), " TOP ",n," CLONOTYPES (",cname,")", sep = ""))+
      scale_fill_manual(values = gen_colors(color_conditions$colorful, nrow(plotx)))+
      theme(plot.margin = unit(c(2,2,2,2), "cm"),
            axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 10,14)),
            legend.key.size = unit(1, "cm"),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 30, face=2),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "bottom")+
      guides(color=guide_legend(title="SUBTISSUE", ncol = 2),fill = guide_legend(ncol = ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 1,2)))+xlab("SUB-TISSUES")+
      ylab("PROPORTION")
    somePNGPath = paste(output_dir,"5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE/5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".png", sep = "")
    png(somePNGPath, width=6000, height=ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 6000,4000), units = "px", res = 300)
    print(p1)
    dev.off()
  }
}

abundance_contig <- readRDS("~/Desktop/abundance_contig_TBcells_5.2SCA_PLOT_ALLUVIAL.RDS")

n <- 10
for(k in 1:length(abundance_contig)){
  current <- abundance_contig[[k]]
  cell_types <- as.character(sort(unique(current$CELL_SUBTYPE)))
  for(j in 1:length(cell_types)){
    plotx <- current[current$CELL_SUBTYPE == cell_types[j],]
    plotx <- plotx[order(plotx$Freq, decreasing = T),]
    ctop <- unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])[1:ifelse(length(unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])) < n, length(unique(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)])), n)]
    plotx <- plotx[which(plotx[,grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] %in% ctop),]
    cname <- colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)]
    cct <- as.character(unique(plotx$CELL_SUBTYPE))
    colnames(plotx)[grep(paste(clonecall_cols,collapse = "|"), colnames(plotx), ignore.case = T)] <- "CLONE_TYPE"
    plotx <- dcast(plotx, CLONE_TYPE~SUBTISSUE, value.var = "PROPORTION")
    plotx[is.na(plotx)] <- 0
    plotx$CLONE_TYPE <- gsub("^NA;|;NA$","",plotx$CLONE_TYPE)
    plotx$CLONE_TYPE <- gsub(";NA;",";",plotx$CLONE_TYPE)
    plotx$CLONE_TYPE <- gsub("^None;|;None$","",plotx$CLONE_TYPE, ignore.case = T)
    plotx$CLONE_TYPE <- gsub(";NONE;",";",plotx$CLONE_TYPE, ignore.case = T)
    
    write.table(plotx,paste(output_dir,"5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE/5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    class(plotx) <- c("immunr_dynamics","data.table","data.frame")
    p1 <- vis(plotx)+ggtitle(paste(gsub("CELL$","CELLS",cct, ignore.case = T), " ONLY: ",gsub("(.*?)_(.*)","\\1",abundance_names[k]), " TOP ",n," CLONOTYPES (",cname,")", sep = ""))+
      scale_fill_manual(values = gen_colors(color_conditions$colorful, nrow(plotx)))+
      theme(plot.margin = unit(c(2,2,2,2), "cm"),
            axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 10,14)),
            legend.key.size = unit(1, "cm"),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 30, face=2),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "bottom")+
      guides(color=guide_legend(title="SUBTISSUE", ncol = 2),fill = guide_legend(ncol = ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 1,2)))+xlab("SUB-TISSUES")+
      ylab("PROPORTION")
    somePNGPath = paste(output_dir,"5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE/5.5SCA_PLOT_ALLUVIAL_TOP10_EACH_CELL_SUBTYPE_",gsub("\\s+|\\/","_",cct), "_",abundance_names[k],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_",cname,".png", sep = "")
    png(somePNGPath, width=6000, height=ifelse(length(grep("CTStrict|Ctnt", cname, ignore.case = T)) > 0, 6000,4000), units = "px", res = 300)
    print(p1)
    dev.off()
  }
}

# for(k in 1:length(abundance_contig)){
#   current <- abundance_contig[[k]]
#   cell_types <- unique(current$CELL_SUBTYPE)
#   for(j in 1:length(cell_types)){
#     plotx <- current[which(current$CELL_SUBTYPE == cell_types[j]),]
#     plotx <- dcast(plotx, CTaa ~SUBTISSUE, value.var = "PROPORTION")
#     plotx[is.na(plotx)] <- 0
#     row.names(plotx) <- plotx$CTaa
#     plotx <- plotx[,grep("CTaa",colnames(plotx),invert = T)]
#     plotx <- scale(plotx)
#     plotx <- t(scale(t(plotx)))
#     heatmap.2(as.matrix(t(plotx)), trace="none",key=T, keysize=1,
#               col=jet2.col (n = 100, alpha = 1),srtCol=45, scale="none", Colv = T,
#               Rowv = T,
#               density.info="none", cexCol=0.8,cexRow=0.8)
#     
#   }
#   
#   
# }


######## Find top genes and do pathway analysis
rnaseq_coords$CELL_ID
table(x$TCRBCR)
head(rnaseq_coords)
dim(rnaseq_coords)

which(paste(rnaseq_coords$SAMPLE_ID,rnaseq_coords$CELL_ID, sep = "") %in% 
      paste(x$SAMPLE_ID,x$CELL_ID, sep = ""))

mapp_ref <- data.frame(SID = paste(rnaseq_coords$SAMPLE_ID,rnaseq_coords$CELL_ID, sep = "_"),
                       CELL_SUBTYPE = x@meta.data[match(paste(rnaseq_coords$SAMPLE_ID,rnaseq_coords$CELL_ID, sep = ""), paste(x$SAMPLE_ID,x$CELL_ID, sep = "")),"CELL_SUBTYPE"])

mapp_ref$CELL_SUBTYPE[which(is.na(mapp_ref$CELL_SUBTYPE))] <- "OTHERS"
saveRDS(mapp_ref, "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/Original_cell_order_final_tcrbcr_subtypes_20210725.RDS")

row.names(x@meta.data) <- colnames(x)
# x2 <- subset(x, subset = TCRBCR != "B-CELLS")
# x2 <- subset(x2, subset = TCRBCR != "OTHER CELLS")
rnaseq_coords$CELL_TYPE <- x@meta.data[match(paste(rnaseq_coords$CELL_ID,rnaseq_coords$SAMPLE_ID),paste(x$CELL_ID,x$SAMPLE_ID)),"CELL_SUBTYPE"]
rnaseq_coords[which(is.na(rnaseq_coords$CELL_TYPE)),"CELL_TYPE"] <- "TO_REMOVE"
saveRDS(rnaseq_coords,"~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/RDATA_RDS/FINAL_TCRBCR_CELL_SUBTYPES_RNASEQ_COORDS.RDS")

# data_integrated <- subset(data_integrated, subset = CELL_SUBTYPE != "OTHERS")
# data_integrated <- subset(data_integrated, subset = CELL_TYPE != "TO_REMOVE")

#### Top markers: Cell Types for Each Tissue
celltypes_tissues <- data.frame(table(x@meta.data[,c("TISSUE","CELL_SUBTYPE")]))
plotx <- reshape2::dcast(celltypes_tissues, TISSUE ~ CELL_SUBTYPE, value.var = "Freq")
write.table(plotx, "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/FINAL_CELL_TYPE_CELL_NUM_TISSUE_WISE.txt", quote = F, row.names = F, sep = "\t")

celltypes_tissues <- celltypes_tissues[celltypes_tissues$Freq > 0,]
celltypes_tissues$ID <- paste(celltypes_tissues$TISSUE, celltypes_tissues$CELL_SUBTYPE, sep = ":")
to_remove <- celltypes_tissues[which(celltypes_tissues$Freq <= 100),"ID"]
# plotx <- reshape2::dcast(celltypes_tissues, TISSUE ~ CELL_SUBTYPE, value.var = "Freq")
# row.names(plotx) <- plotx$TISSUE
# plotx <- plotx[,grep("^TISSUE$", colnames(plotx), ignore.case = T, invert = T)]
#     heatmap.2(as.matrix(t(plotx)), trace="none",key=T, keysize=1,
#               col=jet2.col (n = 100, alpha = 1),srtCol=45, scale="none", Colv = T,
#               Rowv = T,
#               density.info="none", cexCol=0.8,cexRow=0.8)

files <- list.files(path = "~/Dropbox/KI/Studies/LXXLZ/Website/10X/IMMUNE_PROFILING/Data/SAMPLE_DATA/OUTPUT/RDATA_RDS/",
                    pattern = "_FindAllMarkers_Cell_Type_Result_RNA_Assay_Integrated_Data.RDS", full.names = T)
files <- files[grep("TCRBCR_FindAllMarkers_Cell_Type_Result_RNA_Assay_Integrated_Data", files, ignore.case = T, invert = T)]

for(j in 1:length(files)){
  cname <- gsub(".*\\/(.*)_FindAllMarkers_Cell_Type_Result_RNA_Assay_Integrated_Data.RDS","\\1",files[j])
  current_markers <- readRDS(files[j])
  current_markers <- current_markers[current_markers$p_val_adj < 0.05,]
  current_markers <- current_markers[order(current_markers$avg_log2FC, decreasing = T),]
  current_markers <- current_markers[which(!paste(cname,current_markers$cluster, sep = ":") %in% to_remove),]
  write.table(current_markers, paste(output_dir,"FILTERED_P005_LOW_CELL_NUM_100_TISSUE_",cname,"_FindAllMarkers.txt", sep = ""), quote = F, row.names = F, sep = "\t")
  n <- 50
  topn <- split(current_markers, current_markers$cluster)
  topn <- lapply(topn, function(x){
    x <- x[order(x$avg_log2FC, decreasing = T),]
    x <- x[1:n,]
  })
  topn <- do.call(rbind.data.frame, topn)
  
  n <- 200
  p_threshold <- 0.01
  current_markers <- current_markers[which(current_markers$p_val_adj < p_threshold),]
  clusters <- unique(current_markers$cluster)
  current_go <- paste(output_dir, "/TOP_",n,"_GENES_GO_TERMS_TISSUE_CELLTYPES/", sep = "")
  dir.create(current_go)
  current_kegg <- paste(output_dir, "/TOP_",n,"_GENES_KEGG_PATHWAYS_TISSUE_CELLTYPES/", sep = "")
  dir.create(current_kegg)
  for(i in 1:length(clusters)){
    current <- current_markers[which(current_markers$cluster == clusters[i]),]
    current <- current[order(current$avg_log2FC, decreasing = T),]
    current <- data.frame(SYMBOL = unique(current[,"gene"])[1:ifelse(nrow(current)<n, nrow(current),n)])
    gene_col_provided <- "SYMBOL"
    entrez_out <- select(org.Hs.eg.db, keys = current$SYMBOL, keytype = 'SYMBOL', columns = 'ENTREZID')
    current$ENTREZID <- entrez_out[match(current$SYMBOL, entrez_out$SYMBOL),"ENTREZID"]
    current <- current[which(!is.na(current$ENTREZID)),]
    go <- goana(current$ENTREZID, species="Hs")
    go <- go[go$P.DE < p_threshold,]
    go <- go[order(go$P.DE, decreasing = F),]
    write.csv(go, paste(current_go,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_GO_TERMS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)
    
    n <- 50
    TOP_GO <- topGO(go, ontology =c("BP","CC","MF"), n=n)
    TOP_GO <- TOP_GO[order(TOP_GO$P.DE, decreasing = F),]
    TOP_GO$Term <- factor(TOP_GO$Term, levels = rev(unique(TOP_GO$Term)))
    TOP_GO$NegLogP <- -log(TOP_GO$P.DE)
    
    write.csv(TOP_GO, paste(current_go,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_GO_BPCCMF_TERMS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)
    
    somePNGPath <- paste(current_go,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_GO_BPCCMF_TERMS_P",p_threshold,".png",sep = "")
    png(somePNGPath, width = 4000, height =2000, units = "px", res = 300)
    print(ggplot(TOP_GO, aes(NegLogP, Term, color = Term))+
            geom_point(size = 3, alpha = 0.8, stroke = 1)+
            theme_classic()+theme(legend.position = "none",
                                  title = element_text(face = 2, size = 12)) + 
            ylab("GO TERMS")+
            ggtitle(paste(gsub("_"," ",cname), " - TOP GO TERMS: ",gsub("_"," ",clusters[i]),sep="")))
    dev.off()
    
    keg <- kegga(current$ENTREZID, species="Hs")
    keg <- keg[keg$P.DE < p_threshold,]
    keg <- keg[order(keg$P.DE, decreasing = F),]
    write.csv(keg, paste(current_kegg,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)
    
    n <- 50
    TOP_KEGG <- topKEGG(keg, n=n)
    write.csv(TOP_KEGG, paste(current_kegg,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)
    
    TOP_KEGG <- TOP_KEGG[order(TOP_KEGG$P.DE, decreasing = F),]
    TOP_KEGG$Pathway <- factor(TOP_KEGG$Pathway, levels = rev(unique(TOP_KEGG$Pathway)))
    TOP_KEGG$NegLogP <- -log(TOP_KEGG$P.DE)
    
    somePNGPath <- paste(current_kegg,"/",cname,"_",gsub("\\/","_",clusters[i]),"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".png",sep = "")
    png(somePNGPath, width = 4000, height =2000, units = "px", res = 300)
    print(ggplot(TOP_KEGG, aes(NegLogP, Pathway, color = Pathway))+
            geom_point(size = 3, alpha = 0.8, stroke = 1)+
            theme_classic()+theme(legend.position = "none",
                                  title = element_text(face = 2, size = 12)) +
            ggtitle(paste(gsub("_"," ",cname), " - TOP KEGG PATHWAYS: ",gsub("_"," ",clusters[i]),sep="")))
    dev.off()
    
    
  }
}
  
  save.image(paste(output_dir,"/",project_name,"_20210626.RData", sep = ""))
  
#   #### Top markers: clusters
  current_markers <- readRDS(paste(output_dir,"RDATA_RDS/FindAllMarkers_Result_RNA_Assay_Integrated_Data.RDS", sep = ""))
  # current_markers <- NULL
  # files <- list.files("OUTPUT/BEFORE_FINAL/ALL_SIG_GENES_TCELL_CLUSTERS/", full.names = T, pattern = "SIG_GENES_T_CELL_CLUSTER")
  # for(i in 1:length(files)){
  #   current <- read.table(files[i], header=T)
  #   current_markers <- rbind(current_markers, current)
  # }

  current_markers <- current_markers[order(current_markers$avg_log2FC, decreasing = T),]

  n <- 200
  p_threshold <- 0.01
  current_markers <- current_markers[which(current_markers$p_val_adj < p_threshold),]
  current_dir <- paste(output_dir, "/TOP_",n,"_GENES_GO_TERMS_INTEGRATIVE_DATA_RNA_ASSAY_CLUSTERS/", sep = "")
  dir.create(current_dir)
  clusters <- unique(current_markers$cluster)
  for(i in 1:length(clusters)){
    current <- current_markers[which(current_markers$cluster == clusters[i]),]
    current <- current[order(current$avg_log2FC, decreasing = T),]
    current <- data.frame(SYMBOL = unique(current[,"gene"])[1:ifelse(nrow(current)<n, nrow(current),n)])
    gene_col_provided <- "SYMBOL"
    entrez_out <- select(org.Hs.eg.db, keys = current$SYMBOL, keytype = 'SYMBOL', columns = 'ENTREZID')
    current$ENTREZID <- entrez_out[match(current$SYMBOL, entrez_out$SYMBOL),"ENTREZID"]
    current <- current[which(!is.na(current$ENTREZID)),]
    go <- goana(current$ENTREZID, species="Hs")
    go <- go[go$P.DE < p_threshold,]
    go <- go[order(go$P.DE, decreasing = F),]
    write.csv(go, paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_GO_TERMS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)

    n <- 50
    TOP_GO <- topGO(go, ontology =c("BP","CC","MF"), n=n)
    TOP_GO <- TOP_GO[order(TOP_GO$P.DE, decreasing = F),]
    TOP_GO$Term <- factor(TOP_GO$Term, levels = rev(unique(TOP_GO$Term)))
    TOP_GO$NegLogP <- -log(TOP_GO$P.DE)

    write.csv(TOP_GO, paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_GO_BPCCMF_TERMS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)

    somePNGPath <- paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_GO_BPCCMF_TERMS_P",p_threshold,".png",sep = "")
    png(somePNGPath, width = 3500, height =2000, units = "px", res = 300)
    print(ggplot(TOP_GO, aes(NegLogP, Term, color = Term))+
            geom_point(size = 3, alpha = 0.8, stroke = 1)+
            theme_classic()+theme(legend.position = "none",
                                  title = element_text(face = 2, size = 12)) +
            ylab("GO TERMS")+
            ggtitle(paste("TOP GO TERMS: CLUSTER ",clusters[i],sep="")))
    dev.off()

    keg <- kegga(current$ENTREZID, species="Hs")
    keg <- keg[keg$P.DE < p_threshold,]
    keg <- keg[order(keg$P.DE, decreasing = F),]
    write.csv(keg, paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)

    n <- 50
    TOP_KEGG <- topKEGG(keg, n=n)
    write.csv(TOP_KEGG, paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".csv",sep = ""), quote = F, row.names = T)

    TOP_KEGG <- TOP_KEGG[order(TOP_KEGG$P.DE, decreasing = F),]
    TOP_KEGG$Pathway <- factor(TOP_KEGG$Pathway, levels = rev(unique(TOP_KEGG$Pathway)))
    TOP_KEGG$NegLogP <- -log(TOP_KEGG$P.DE)

    somePNGPath <- paste(current_dir,"/CLUSTER_",clusters[i],"_TOP_",n,"_GENES_KEGG_PATHWAYS_P",p_threshold,".png",sep = "")
    png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
    print(ggplot(TOP_KEGG, aes(NegLogP, Pathway, color = Pathway))+
            geom_point(size = 3, alpha = 0.8, stroke = 1)+
            theme_classic()+theme(legend.position = "none",
                                  title = element_text(face = 2, size = 12)) +
            ggtitle(paste("TOP KEGG PATHWAYS: CLUSTER ",clusters[i],sep="")))
    dev.off()


  }



# library("gplots")
# library("plot3D")
# library("org.Hs.eg.db")
# library("limma")
# somePNGPath = paste("OUTPUT/TOP_",n,"_GENES_TCRBCR_CELL_SUBTYPES.png", sep = "")
# png(somePNGPath, width=2000, height=n*220, units = "px", res = 150)
# heatmap.2(plot_median_CELL_SUBTYPE,margin=c(30, 15), trace="none",key=T, keysize=1,
#           col=jet2.col (n = 100, alpha = 1),main = project_name,
#           srtCol=45,
#           scale="none", Colv = T, Rowv = T,
#           density.info="none", cexCol=2,cexRow=2)
# dev.off()









