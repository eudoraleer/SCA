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

pheno_data <- pheno_ini(phenodata_dir, pipeline = "SCIENCE")
#####################################################################################################
# library(stringr)
# multitissue_files <- pheno_data[pheno_data$TISSUE == "MULTI_TISSUES",]
# pheno_data <- pheno_data[!pheno_data$TISSUE == "MULTI_TISSUES",]
# 
# for(i in 1:nrow(multitissue_files)){
#   x <- readfile(paste("CONTIGS/",multitissue_files[i,"FILE"], sep = ""))
#   current_tissues <- unique(x$TISSUE)
#   for(j in 1:length(current_tissues)){
#     current <- x[x$TISSUE == current_tissues[j],]
#     current_meta <- multitissue_files[i,]
#     current_meta$SAMPLE_ID <- paste(multitissue_files[i,"SAMPLE_ID"],"_",current_tissues[j], sep = "")
#     current_meta$FILE <- paste(multitissue_files[i,"SAMPLE_ID"],"_",current_tissues[j], ".csv",sep = "")
#     current_meta$TISSUE <- current_tissues[j]
#     current_meta$SUBTISSUE <- current_tissues[j]
#     # paste(str_to_title(gsub("_"," ",current_tissues[j])),collapse = "_")
#     pheno_data <- rbind(pheno_data, current_meta)
#     write.csv(current,paste("CONTIGS/",multitissue_files[i,"SAMPLE_ID"],"_",current_tissues[j],".csv",sep=""),sep = ",",quote = F,row.names=F)
#   rm(current)
#   rm(current_meta)
#     }
#   system(paste("mv ", paste("CONTIGS/",multitissue_files[i,"FILE"], sep = "")," MULTITISSUES/", sep = ""))
# rm(x)
# rm(current_tissues)
# }

# write.csv(pheno_data, "TCRBCR_ALL_TISSUES.csv", row.names = F, quote=F)

pheno_data$TISSUE <- gsub("PBMC","BLOOD",pheno_data$TISSUE)
pheno_data$TISSUE <- gsub("CORD_BLOOD","BLOOD",pheno_data$TISSUE)
write.csv(pheno_data, "TCRBCR_ALL_TISSUES_EXCLUDING_NUM_94.csv", row.names = F, quote=F)

x <- readfile("CONTIGS/NUM_94_GSE149356_Processed_TCR.Final.csv")
current <- data.frame(table(x[,c("donor","TISSUE")]))
current <- current[which(!current$Freq == 0),]

summary <- data.frame(TISSUE = ifelse(current$TISSUE %in% c("CORD_BLOOD","PBMC"),"BLOOD",current$TISSUE),
                      SUBTISSUE = current$TISSUE,
                      CELL_COUNT = current$Freq,
                      FILE = "NUM_94_GSE149356_Processed_TCR.Final.csv")
k <- NULL

for(i in 1:nrow(pheno_data)){
  x <- readfile(paste("CONTIGS/",pheno_data$FILE[i], sep = ""))
  print(unique(x$TISSUE))
  if(!is.null(unique(x$TISSUE))){
  if(length(unique(x$TISSUE)) > 1){
    k <- c(k, pheno_data$FILE[i])
  }
  }
  summary <- rbind(summary,data.frame(TISSUE = ifelse(pheno_data$TISSUE[i] %in% c("CORD_BLOOD","PBMC"),"BLOOD",pheno_data$TISSUE[i]),
                                      SUBTISSUE = pheno_data$SUBTISSUE[i],
                                      CELL_COUNT = length(unique(x$barcode)),
                                      FILE = pheno_data$FILE[i]))
  rm(x)
}

write.csv(summary,"OUTPUT/scIMMUNE_SUMMARY_CELL_COUNTS_RAW.csv", quote=F, row.names = F)

current <- data.frame(table(pheno_data[,c("TISSUE","SUBTISSUE")]))
current <- current[which(!current$Freq == 0),]
write.csv(current, "OUTPUT/scIMMUNE_SUBTISSUE_SAMPLE_COUNTS_RAW.csv", row.names = F, quote = F)
write.csv(data.frame(table(pheno_data$TISSUE)), "OUTPUT/scIMMUNE_TISSUE_SAMPLE_COUNTS_RAW.csv", row.names = F, quote = F)

contig_files <- paste(input_dir,"CONTIGS/",pheno_data$FILE,sep = "")
k <- 0

not_formated <- NULL

for(i in 1:length(contig_files)){
  current <- readfile(contig_files[i])
  if(length(grep("is_cell",colnames(current), ignore.case = T)) == 0){
    print(colnames(current))
    print(paste(i,":",contig_files[i]))
    k <- k+1
    if(length(grep("NUM_94_",contig_files[i],ignore.case = T)) > 0){
      not_formated <- rbind(not_formated, data.frame(FILE = contig_files[i],
                       SAMPLE_ID = gsub(".*\\/(.*)\\.csv","\\1",contig_files[i]),
                       PROJECT = "P210", TISSUE = "BLOOD", SUBTISSUE = "PBMC|CORD_BLOOD",CELL_TYPE = "TCR"))
    }else{
    not_formated <- rbind(not_formated, pheno_data[match(gsub(".*\\/(.*)","\\1",contig_files[i]),pheno_data$FILE),])
    }
    }
}
current <- data.frame(table(not_formated[,c("TISSUE","SUBTISSUE","CELL_TYPE")]))
current <- current[which(!current$Freq == 0),]

write.csv(not_formated,"OUTPUT/NOT_FORMATED.csv", row.names = F, quote=F)

####################################################################################################
contig_path_monarch <- unique(gsub("(.*\\/).*","\\1",contig_files))
contig_monarch <- repLoad(contig_path_monarch,.mode = "paired")
# repLoad("~/Desktop/test_contigs/")

filtered_pheno <- pheno_data[match(names(contig_monarch$data), gsub("(.*)\\.csv","\\1",pheno_data$FILE, ignore.case = T)),]
write.csv(filtered_pheno, paste("FILTERED_TCRBCR_ALL_TISSUES.csv"), quote = F, row.names = F)

not_formated <- pheno_data[which(!pheno_data$SAMPLE_ID %in% filtered_pheno$SAMPLE_ID),]
write.csv(not_formated, paste("OUTPUT/NOT_FORMATED.csv"), quote = F, row.names = F)

names(contig_monarch$data) <- filtered_pheno$SAMPLE_ID
contig_monarch$meta <- filtered_pheno[match(names(contig_monarch$data), filtered_pheno$SAMPLE_ID),]

if(median(nchar(names(contig_monarch$data))) > 10){
  contig_monarch$meta$Sample <- paste(contig_monarch$meta$SUBTISSUE,"_S",1:nrow(contig_monarch$meta),"_",contig_monarch$meta$CELL_TYPE, sep = "")
  names(contig_monarch$data) <- contig_monarch$meta$Sample
}else{
  contig_monarch$meta$Sample <- pheno_data[match(names(contig_monarch$data), pheno_data$SAMPLE_ID),"SAMPLE_ID"]
  names(contig_monarch$data) <- contig_monarch$meta$Sample
}

write.csv(contig_monarch$meta, paste(output_dir,"0.0SCA_TABLE_META_ALL_SAMPLES_",project_name,".csv", sep = ""), quote = F, row.names = F)

filtered_pheno$Sample <- contig_monarch$meta[match(filtered_pheno$SAMPLE_ID,contig_monarch$meta$SAMPLE_ID),"Sample"]
write.csv(filtered_pheno, paste("FILTERED_TCRBCR_ALL_TISSUES.csv"), quote = F, row.names = F)

meta_explore <- repExplore(contig_monarch$data, .method = "volume")
write.csv(meta_explore, paste(output_dir,"0.0FILTERED_SAMPLES_UNIQUE_CLONOTYPES.csv"), quote = F, row.names = F)

meta <- contig_monarch$meta
meta_explore$CELL_TYPE <- meta[match(meta_explore$Sample,meta$Sample),"CELL_TYPE"]
meta_explore$TISSUE <- meta[match(meta_explore$Sample,meta$Sample),"TISSUE"]
meta_explore$SUBTISSUE <- meta[match(meta_explore$Sample,meta$Sample),"SUBTISSUE"]
meta_explore$SUBTISSUE <- paste(meta_explore$TISSUE, ":", meta_explore$SUBTISSUE, sep = "")

current <- meta_explore[meta_explore$CELL_TYPE == "TCR",]
# current_meta <- meta[which(meta$Sample %in% current$Sample),]
current$Sample <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", current$Sample, sep = "")
p1 <- vis(current, .test = F)+
  # scale_fill_manual(values = gen_colors(color_conditions$ggplot, length(unique(meta$Sample))))+
  facet_wrap(~CELL_TYPE) +
  theme(axis.text.x = element_text(size = 7),
                                 axis.text.y = element_text(size = 15),
                                 axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 
                                 legend.text = element_text(size = 15),
                                 strip.text.x = element_text(size = 20),
                                 legend.position = "none")

current <- meta_explore[meta_explore$CELL_TYPE == "BCR",]
# current_meta <- meta[which(meta$Sample %in% current$Sample),]
current$Sample <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", current$Sample, sep = "")
p2 <- vis(current, .test = F)+
  # scale_fill_manual(values = gen_colors(color_conditions$ggplot, length(unique(meta$Sample))))+
  facet_wrap(~CELL_TYPE) +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

library(patchwork)
somePNGPath = paste(output_dir,"01.0SCA_PLOT_NUMBER_UNIQUE_CLONOTYPE_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
print(p1+p2)
dev.off()

current <- meta_explore[meta_explore$CELL_TYPE == "TCR",]
current_meta <- meta[which(meta$Sample %in% current$Sample),]
# current$Sample <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", names(current), sep = "")
p1 <- vis(current, .by = c("TISSUE"), .meta = current_meta, .test = F)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$TISSUE))))+
  facet_wrap(~CELL_TYPE) + theme(axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 15),
                                 axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 legend.text = element_text(size = 15),
                                 strip.text.x = element_text(size = 20),
                                 legend.position = "none")

current <- meta_explore[meta_explore$CELL_TYPE == "BCR",]
current_meta <- meta[which(meta$Sample %in% current$Sample),]
# current$Sample <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", names(current), sep = "")
p2 <- vis(current, .by = c("TISSUE"), .meta = current_meta, .test = F)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$TISSUE))))+
  facet_wrap(~CELL_TYPE) + theme(axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 15),
                                 axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 legend.text = element_text(size = 15),
                                 strip.text.x = element_text(size = 20),
                                 legend.position = "none")

somePNGPath = paste(output_dir,"01.0SCA_PLOT_NUMBER_UNIQUE_CLONOTYPE_TISSUES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=5500, height=2500, units = "px", res = 300)
print(p1+p2)
dev.off()

current <- meta_explore[meta_explore$CELL_TYPE == "TCR",]
current_meta <- meta[which(meta$Sample %in% current$Sample),]
current$SUBTISSUE <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", current$SUBTISSUE, sep = "")
current_meta$SUBTISSUE <- paste(meta[match(current_meta$Sample, meta$Sample),"TISSUE"], ":", current_meta$SUBTISSUE, sep = "")

p1 <- vis(current, .by = c("SUBTISSUE"), .meta = current_meta, .test = F)+
  scale_fill_manual(values = gen_colors(color_conditions$tableau20, length(unique(meta$SUBTISSUE))))+
  facet_wrap(~CELL_TYPE)+  theme(axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 15),
                                 axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 
                                 legend.text = element_text(size = 15),
                                 strip.text.x = element_text(size = 20),
                                 legend.position = "none")

current <- meta_explore[meta_explore$CELL_TYPE == "BCR",]
current_meta <- meta[which(meta$Sample %in% current$Sample),]
current$SUBTISSUE <- paste(meta[match(current$Sample, meta$Sample),"TISSUE"], ":", current$SUBTISSUE, sep = "")
current_meta$SUBTISSUE <- paste(meta[match(current_meta$Sample, meta$Sample),"TISSUE"], ":", current_meta$SUBTISSUE, sep = "")

p2 <- vis(current, .by = c("SUBTISSUE"), .meta = current_meta, .test = F)+
  scale_fill_manual(values = gen_colors(color_conditions$tableau20, length(unique(meta$SUBTISSUE))))+
  facet_wrap(~CELL_TYPE)+  theme(axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 15),
                                 axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 
                                 legend.text = element_text(size = 15),
                                 strip.text.x = element_text(size = 20),
                                 legend.position = "none")

somePNGPath = paste(output_dir,"01.0SCA_PLOT_NUMBER_UNIQUE_CLONOTYPE_SUBTISSUES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=5500, height=2500, units = "px", res = 300)
print(p1+p2)
dev.off()

tissues <- unique(meta$TISSUE)
current_colors <- gen_colors(color_conditions$tableau20, length(tissues))
for(i in 1:length(tissues)){
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$TISSUE == tissues[i]),"Sample"])]
current_meta <- meta[which(meta$TISSUE == tissues[i]),]
if(length(unique(current_meta$SUBTISSUE)) > 1){
  cc_temp <- gen_colors(current_colors[i], length(unique(current_meta$SUBTISSUE)))
}else{
  cc_temp <- current_colors[i]
}

exp_len_aa <- repExplore(current, .method = "len", .col = "aa")
exp_len_aa$SUBTISSUE <- meta[match(exp_len_aa$Sample, meta$Sample),"SUBTISSUE"]
exp_len_aa$CELL_TYPE <- meta[match(exp_len_aa$Sample, meta$Sample),"CELL_TYPE"]

exp_len_nt <- repExplore(current, .method = "len", .col = "nt")
exp_len_nt$SUBTISSUE <- meta[match(exp_len_nt$Sample, meta$Sample),"SUBTISSUE"]
exp_len_nt$CELL_TYPE <- meta[match(exp_len_nt$Sample, meta$Sample),"CELL_TYPE"]

current <- NULL
current <- exp_len_aa[which(exp_len_aa$CELL_TYPE == "TCR"),]
temp_meta <- current_meta[which(current_meta$CELL_TYPE == "TCR"),]
pt <- vis(current, .by = c("SUBTISSUE"), .meta = temp_meta, .test = F) + 
  facet_wrap(~CELL_TYPE)+ ggtitle("TCR: AMINO ACID")+
  scale_fill_manual(values = cc_temp) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")+
facet_wrap(~SUBTISSUE)
p1 <- pt
rm(pt)

if(length(which(exp_len_aa$CELL_TYPE == "BCR")) > 0){
current <- exp_len_aa[which(exp_len_aa$CELL_TYPE == "BCR"),]
temp_meta <- current_meta[which(current_meta$CELL_TYPE == "BCR"),]
pb <- vis(current, .by = c("SUBTISSUE"), .meta = temp_meta, .test = F) + 
facet_wrap(~CELL_TYPE)+ggtitle("BCR: AMINO ACID")+
scale_fill_manual(values = cc_temp) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")+
  facet_wrap(~SUBTISSUE)

p1 <- p1 + pb
rm(pb)
}

# if(length(unique(exp_len_aa$CELL_TYPE)) == 1){
#   p1 <- p1+facet_wrap(~SUBTISSUE)+ggtitle(unique(exp_len_aa$CELL_TYPE))+
#     ggtitle(paste("DISTRIBUTION OF CDR3 LENGTH: AMINO ACID\n\n",unique(exp_len_aa$CELL_TYPE), sep = ""))+
#     theme(plot.title = element_text(hjust = 0.5))
# }else{
  p1 <- p1+
    # facet_wrap(~SUBTISSUE)+
    plot_annotation(title = "DISTRIBUTION OF CDR3 LENGTH: AMINO ACID", 
    theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
# }


  current <- NULL
  current <- exp_len_nt[which(exp_len_nt$CELL_TYPE == "TCR"),]
  temp_meta <- current_meta[which(current_meta$CELL_TYPE == "TCR"),]
  pt <- vis(current, .by = c("SUBTISSUE"), .meta = temp_meta, .test = F) + 
    facet_wrap(~CELL_TYPE)+ggtitle("TCR: NUCLEOTIDE")+
    scale_fill_manual(values = cc_temp) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")+
    facet_wrap(~SUBTISSUE)
  p2 <- pt
  rm(pt)
  
  if(length(which(exp_len_nt$CELL_TYPE == "BCR")) > 0){
    current <- exp_len_nt[which(exp_len_nt$CELL_TYPE == "BCR"),]
    temp_meta <- current_meta[which(current_meta$CELL_TYPE == "BCR"),]
    pb <- vis(current, .by = c("SUBTISSUE"), .meta = temp_meta, .test = F) + 
      facet_wrap(~CELL_TYPE)+ggtitle("BCR: NUCLEOTIDE")+
      scale_fill_manual(values = cc_temp) +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
            legend.text = element_text(size = 15),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(hjust = 0.5),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "right")+
      facet_wrap(~SUBTISSUE)
    
    p2 <- p2 + pb
    rm(pb)
  }
  
  # if(length(unique(exp_len_nt$CELL_TYPE)) == 1){
  #   p2 <- p2+facet_wrap(~SUBTISSUE)+ggtitle(unique(exp_len_nt$CELL_TYPE))+
  #     ggtitle(paste("DISTRIBUTION OF CDR3 LENGTH: NUCLEOTIDE\n\n",unique(exp_len_nt$CELL_TYPE), sep = ""))+
  #     theme(plot.title = element_text(hjust = 0.5))
  # }else{
  p2 <- p2+
    # facet_wrap(~SUBTISSUE)+
    plot_annotation(title = "DISTRIBUTION OF CDR3 LENGTH: NUCLEOTIDE", 
                    theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
  # }
  
somePNGPath = paste(output_dir,"01.1SCA_PLOT_DISTR_CDR3_LENGTHS_",tissues[i],"_",project_name,".png", sep = "")
png(somePNGPath, width=5000, height=2800, units = "px", res = 300)
print(p1/p2+plot_annotation(title = paste("Distribution of CDR3 Lengths - ", tissues[i], sep = ""),
theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
dev.off()

rm(p1)
rm(p2)

# imm_top <- data.frame(repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000)))
# imm_top$SUBTISSUE <- meta[match(row.names(imm_top), meta$Sample),"SUBTISSUE"]
# imm_top$CELL_TYPE <- meta[match(row.names(imm_top), meta$Sample),"CELL_TYPE"]
# imm_top <- as.matrix(imm_top)
# class(imm_top) <- c("immunr_top_prop","matrix","array")
}

exp_cnt <- repExplore(contig_monarch$data, .method = "count")
exp_cnt$CELL_TYPE <- meta[match(exp_cnt$Sample, meta$Sample),"CELL_TYPE"]

# exp_vol <- repExplore(contig_monarch$data, .method = "volume")

p <- vis(exp_cnt, .by = c("SUBTISSUE"), .meta = contig_monarch$meta, .test = F) + 
      guides(color = guide_legend(ncol = 1)) + facet_wrap(~CELL_TYPE)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")+
      scale_color_manual(values = gen_colors(color_conditions$colorful, length(unique(meta$SUBTISSUE))))
# p3 <- vis(exp_vol, .by = c("SUBTISSUE"), .meta = contig_monarch$meta, .test = F)

somePNGPath = paste(output_dir,"01.2SCA_PLOT_CLONOTYPE_ABUNDANCE_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=5000, height=2000, units = "px", res = 300)
print(p)
dev.off()

# p4 <- vis(exp_len_aa, .by = "CELL_TYPE", .meta = contig_monarch$meta)
# p5 <- vis(exp_cnt, .by = "GROUP", .meta = contig_monarch$meta)
# p6 <- vis(exp_vol, .by = c("GROUP", "CELL_TYPE"), .meta = contig_monarch$meta)
# proportion of repertoire occupied by the pools of cell clones

imm_pr <- repClonality(contig_monarch$data, .method = "clonal.prop")
write.csv(imm_pr, paste(output_dir,"01.3SCA_TABLE_CLONE_PROPORTIONS_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)

# Top most abundance cell types
imm_top <- repClonality(contig_monarch$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
write.csv(imm_top, paste(output_dir,"01.4SCA_TABLE_TOP_N_MOST_ABUNDANT_CLONOTYPES_PORPORTION_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)

imm_rare <- repClonality(contig_monarch$data, .method = "rare")
write.csv(imm_rare, paste(output_dir,"01.5SCA_TABLE_RELATIVE_ABUNDANCE_RARE_CLONOTYPES_PORPORTION_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)

# clonal space homeostasis
# the proportion of the repertoire occupied by the clones of a given size
imm_hom <- repClonality(contig_monarch$data, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
write.csv(imm_hom, paste(output_dir,"01.6SCA_TABLE_CLONALSPACE_HOMEOSTASIS_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)

current <- contig_monarch$data
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))

somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=6500, height=3000, units = "px", res = 300)
p <- vis(imm_top)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - TOP CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))

somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_top)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - TCR TOP CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))

somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_top)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - BCR TOP CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_GROUPED_BY_TOP_CLONOGROUPS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_top, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - TCR TOP CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_GROUPED_BY_TOP_CLONOGROUPS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_top, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - BCR TOP CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data
imm_top <- repClonality(current, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
somePNGPath = paste(output_dir,"01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_GROUPED_BY_TOP_CLONOGROUPS_ALL_SAMPLES_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_top, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - ALL TOP CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_rare <- repClonality(current, .method = "rare")

somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_rare)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - RARE CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_rare <- repClonality(current, .method = "rare")

somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_rare)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - TCR RARE CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_rare <- repClonality(current, .method = "rare")

somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_rare)+scale_fill_manual(values = color_conditions$bright) +
  ggtitle("scIMMUNE - BCR RARE CLONAL PROPORTIONS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
imm_rare <- repClonality(current, .method = "rare")
somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_GROUPED_BY_RARE_CLONOGROUPS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_rare, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - TCR RARE CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
imm_rare <- repClonality(current, .method = "rare")
somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_GROUPED_BY_RARE_CLONOGROUPS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_rare, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - BCR RARE CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data
imm_rare <- repClonality(current, .method = "rare")
somePNGPath = paste(output_dir,"01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_GROUPED_BY_RARE_CLONOGROUPS_ALL_SAMPLES_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_rare, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - ALL RARE CLONAL PROPORTIONS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current_colors <- c("#ff6666","#ccff66","#5D2E8C","#2EC4B6")

current <- contig_monarch$data
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))

somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_hom)+scale_fill_manual(values = current_colors) +
  ggtitle("scIMMUNE - CLONALSPACE HOMEOSTASIS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))

somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_hom)+scale_fill_manual(values = current_colors) +
  ggtitle("scIMMUNE - TCR CLONALSPACE HOMEOSTASIS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
names(current) <- paste(meta[match(names(current), meta$Sample),"TISSUE"], ":", names(current), sep = "")
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))

somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
p <- vis(imm_hom)+scale_fill_manual(values = current_colors) +
  ggtitle("scIMMUNE - BCR CLONALSPACE HOMEOSTASIS")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_GROUPED_BY_RARE_CLONOGROUPS_TCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_hom, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - TCR CLONALSPACE HOMEOSTASIS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_GROUPED_BY_RARE_CLONOGROUPS_BCR_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_hom, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - BCR CLONALSPACE HOMEOSTASIS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

current <- contig_monarch$data
imm_hom <- repClonality(current, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
somePNGPath = paste(output_dir,"01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_GROUPED_BY_RARE_CLONOGROUPS_ALL_SAMPLES_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=2000, units = "px", res = 200)
p <- vis(imm_hom, .by = c("SUBTISSUE"), .meta = contig_monarch$meta)+
  scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
  ggtitle("scIMMUNE - CLONALSPACE HOMEOSTASIS")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
print(p)
dev.off()

# Repertoire overlap and public (shared) clonotypes
imm_ov1 <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
names(imm_ov1) <- paste(meta[match(names(imm_ov1), meta$Sample),"TISSUE"], ":", names(imm_ov1), sep = "")
imm_ov1 <- repOverlap(imm_ov1, .method = "public", .verbose = F)
p <- vis(imm_ov1) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("CLONOTYPES SHARED BETWEEN TCR REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 50, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 50, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 30),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 50, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")

somePNGPath = paste(output_dir,"02.0SCA_PLOT_TCR_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=6000, units = "px", res = 200)
print(p)
dev.off()

imm_ov1 <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
names(imm_ov1) <- paste(meta[match(names(imm_ov1), meta$Sample),"TISSUE"], ":", names(imm_ov1), sep = "")
imm_ov1 <- repOverlap(imm_ov1, .method = "public", .verbose = F)
p <- vis(imm_ov1) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("CLONOTYPES SHARED BETWEEN BCR REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 40, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 40, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 30, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")

somePNGPath = paste(output_dir,"02.0SCA_PLOT_BCR_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=6000, units = "px", res = 300)
print(p)
dev.off()

imm_ov1 <- contig_monarch$data
names(imm_ov1) <- paste(meta[match(names(imm_ov1), meta$Sample),"TISSUE"], ":", names(imm_ov1), sep = "")
imm_ov1 <- repOverlap(imm_ov1, .method = "public", .verbose = F)
p <- vis(imm_ov1) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("CLONOTYPES SHARED BETWEEN ALL REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 40, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 40, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 40, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
# p$layers[[2]] <- NULL
somePNGPath = paste(output_dir,"02.0SCA_PLOT_ALL_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=6000, units = "px", res = 200)
print(p)
dev.off()

# Repertoire overlap for morisita (shared) clonotypes
imm_ov2 <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
names(imm_ov2) <- paste(meta[match(names(imm_ov2), meta$Sample),"TISSUE"], ":", names(imm_ov2), sep = "")
imm_ov2 <- repOverlap(imm_ov2, .method = "morisita", .verbose = F)
p <- vis(imm_ov2) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("MORISITA OVERLAP INDEX - TCR REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 50, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 50, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 50, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
# guides(fill=guide_legend(keywidth=1,keyheight=1,default.unit="inch", title="Overlap values"))
p$layers[[2]] <- NULL

somePNGPath = paste(output_dir,"02.1SCA_PLOT_TCR_MORISITA_OVERLAP_INDEX_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=8000, height=8000, units = "px", res = 300)
print(p)
dev.off()

imm_ov2 <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
names(imm_ov2) <- paste(meta[match(names(imm_ov2), meta$Sample),"TISSUE"], ":", names(imm_ov2), sep = "")
imm_ov2 <- repOverlap(imm_ov2, .method = "morisita", .verbose = F)
p <- vis(imm_ov2) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("MORISITA OVERLAP INDEX - BCR REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 40, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 40, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 30, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")

p$layers[[2]] <- NULL

somePNGPath = paste(output_dir,"02.1SCA_PLOT_BCR_MORISITA_OVERLAP_INDEX_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=6000, units = "px", res = 300)
print(p)
dev.off()

imm_ov2 <- contig_monarch$data
names(imm_ov2) <- paste(meta[match(names(imm_ov2), meta$Sample),"TISSUE"], ":", names(imm_ov2), sep = "")
imm_ov2 <- repOverlap(imm_ov2, .method = "morisita", .verbose = F)
p <- vis(imm_ov2) + 
  scale_fill_continuous_tableau(palette = "Blue-Teal")+
  # scale_fill_gradientn(colours = jet.col (n = 100, alpha = 1))+
  ggtitle("MORISITA OVERLAP INDEX - ALL REPERTOIRES") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 40, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 40, margin=margin(0,10,0,0)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 40, face = 2, hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "right")
p$layers[[2]] <- NULL
somePNGPath = paste(output_dir,"02.1SCA_PLOT_ALL_MORISITA_OVERLAP_INDEX_REPERTOIRE_OVERLAPS_",project_name,".png", sep = "")
png(somePNGPath, width=6000, height=6000, units = "px", res = 200)
print(p)
dev.off()

# Build public repertoire table using CDR3 nucleotide sequences
pr.nt <- pubRep(contig_monarch$data, "nt", .verbose = F)
write.csv(pr.nt, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_ALL_CDR3_NUCLEOTIDES_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

# Public repertoire table: CDR3 aminoacid sequences and V alleles
pr.aav <- pubRep(contig_monarch$data, "aa+v", .verbose = F)
write.csv(pr.aav, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_ALL_CDR3_AA_V_ALLELES_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

pr.aa <- pubRep(contig_monarch$data, "aa", .verbose = F)
write.csv(pr.aa, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_ALL_CDR3_AA_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

# Coding sequences
pr.aav.cod <- pubRep(contig_monarch$data, "aa+v", .coding = T)
write.csv(pr.aav.cod, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_CODING_CDR3_AA_V_ALLELES_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

# Coding sequences
pr.nt.cod <- pubRep(contig_monarch$data, "nt", .coding = T)
write.csv(pr.nt.cod, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_CODING_CDR3_NUCLEOTIDES_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

# Coding sequences
pr.aa.cod <- pubRep(contig_monarch$data, "aa", .coding = T)
write.csv(pr.aa.cod, paste(output_dir,"02.2SCA_TABLE_ALL_CLONOTYPES_BY_CODING_CDR3_AA_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

# Gene usage computation for known gene segments
gene_stats <- gene_stats()
write.csv(gene_stats, paste(output_dir,"02.3SCA_TABLE_KNOWN_GENE_SEGMENTS_BY_SPECIES.csv", sep = ""), quote = F, row.names = F)

all_genes <- colnames(gene_stats)[(!grepl("alias|species", colnames(gene_stats), ignore.case = T)) &
                                    (gene_stats[grep("hs", gene_stats$alias,ignore.case = T),] != 0)]
all_gene_usage <- NULL
for(i in 1:length(all_genes)){
  all_gene_usage <- rbind(all_gene_usage, data.frame(Gene_Type = all_genes[i], geneUsage(contig_monarch$data, .gene = paste("hs.",all_genes[i], sep = ""))))
}

write.csv(all_gene_usage, paste(output_dir,"02.4SCA_TABLE_GENES_USAGE_BY_GENE_TYPE_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

contig_monarch$meta$SUBTISSUE <- paste(contig_monarch$meta$TISSUE,":", contig_monarch$meta$SUBTISSUE, sep = "")
meta$SUBTISSUE <- paste(meta$TISSUE,":", meta$SUBTISSUE, sep = "")
groups <- unique(contig_monarch$meta$CELL_TYPE)

n <- 20
for(i in 1:length(groups)){
  for(j in 1:length(all_genes)){
    somePNGPath = paste(output_dir,"02.5SCA_PLOT_TOP",n,"_",toupper(all_genes[j]),"_GENE_USAGE_",groups[i],"_REPERTOIRE_",project_name,".png", sep = "")
    png(somePNGPath, width=5000, height=2500, units = "px", res = 300)
    current_gu <- geneUsage(contig_monarch$data[which(toupper(names(contig_monarch$data)) %in% toupper(unlist(contig_monarch$meta[contig_monarch$meta$CELL_TYPE == groups[i],"Sample"])))], .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
    current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
    current_gu[is.na(current_gu)] <- 0
    current_gu <- current_gu[grep("^None;.*None$|^None$", current_gu$Names, ignore.case = T, invert = T),]
    
    if(nrow(current_gu) > n){
      current_gu <- current_gu[order(rowMeans(current_gu[,grep("Name", colnames(current_gu), ignore.case = T, invert = T)]), decreasing = T),]
      current_gu <- current_gu[1:n,]
      current_gu <- current_gu[!is.na(current_gu$Names),]
    }else{
      n <- nrow(current_gu)
    }
    
    current_gu$Names <- gsub(";NA|;None|None;|NA;","",current_gu$Names, ignore.case = T)
    # vis(current_gu, .by = "GROUP", .meta = contig_monarch$meta, .plot = "box", .grid = T)
    p1 <- vis(current_gu, .by = "SUBTISSUE", .meta = contig_monarch$meta, .test=F) +
      scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 25, face = 2, hjust = 0.5),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "right")+ guides(fill = guide_legend(ncol = 1))+
      ggtitle(paste("TOP ",n," ",toupper(all_genes[j]), " GENE USAGE: ", groups[i], " REPERTOIRES", sep = ""))
    print(p1)
    n <- 20
    dev.off()
    
  }
}

n <- 20
for(i in 1:length(groups)){
  for(j in 1:length(all_genes)){
    for(k in 1:length(tissues)){
      if(length(which((toupper(names(contig_monarch$data)) %in% toupper(unlist(contig_monarch$meta[contig_monarch$meta$CELL_TYPE == groups[i],"Sample"]))) & names(contig_monarch$data) %in% meta[which(meta$TISSUE == tissues[k]),"Sample"])) > 0){
        current_gu <- geneUsage(contig_monarch$data[which((toupper(names(contig_monarch$data)) %in% toupper(unlist(contig_monarch$meta[contig_monarch$meta$CELL_TYPE == groups[i],"Sample"]))) & names(contig_monarch$data) %in% meta[which(meta$TISSUE == tissues[k]),"Sample"])], .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
        current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
        current_gu[is.na(current_gu)] <- 0
        current_gu <- current_gu[grep("^None;.*None$|^None$", current_gu$Names, ignore.case = T, invert = T),]
        
        if(nrow(current_gu) > n){
          current_gu <- current_gu[order(rowMeans(current_gu[,grep("Name", colnames(current_gu), ignore.case = T, invert = T)]), decreasing = T),]
          current_gu <- current_gu[1:n,]
          current_gu <- current_gu[!is.na(current_gu$Names),]
        }else{
          n <- nrow(current_gu)
        }
        
        current_gu$Names <- gsub(";NA|;None|None;|NA;","",current_gu$Names, ignore.case = T)
        # vis(current_gu, .by = "GROUP", .meta = contig_monarch$meta, .plot = "box", .grid = T)
        if(ncol(current_gu) > 2){
          p1 <- vis(current_gu, .by = "SUBTISSUE", .meta = contig_monarch$meta, .test=F) +
            theme(axis.text.x = element_text(size = 15),
                  axis.text.y = element_text(size = 15),
                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 10),
                  strip.text.x = element_text(size = 15),
                  plot.title = element_text(size = 25, face = 2, hjust = 0.5),
                  strip.background = element_rect(colour="white", fill="white"),
                  legend.position = "right")+ guides(fill = guide_legend(ncol = 1))+
            ggtitle(paste(tissues[k], " TOP ",n," ",toupper(all_genes[j]), " GENE USAGE: ", groups[i], " REPERTOIRES", sep = "")) +
            scale_fill_manual(values = gen_colors(color_conditions$bright, length(unique(meta$SUBTISSUE))))
        }else{
          p1 <- ggplot(current_gu, aes(x=Names, y=Clones, fill=Names)) + 
            geom_bar(stat="identity", alpha=1)+
            xlab("Gene") + ylab("Frequency") + theme_classic() +
            theme(plot.margin = unit(c(1,1,1,1), "cm"),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 15),
                  axis.text.y = element_text(size = 15),
                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 10),
                  strip.text.x = element_text(size = 15),
                  plot.title = element_text(size = 25, face = 2, hjust = 0.5),
                  strip.background = element_rect(colour="white", fill="white"),
                  legend.position = "right")+ guides(fill = guide_legend(ncol = 1))+
            ggtitle(paste(tissues[k], " TOP ",n," ",toupper(all_genes[j]), " GENE USAGE: ", groups[i], " REPERTOIRES", sep = "")) +
            scale_fill_manual(values = gen_colors(color_conditions$colorful, n))
        }
        somePNGPath = paste(output_dir,"02.6GENE_USAGE_TISSUE/02.6SCA_PLOT_",tissues[k],"TOP_",n,"_",toupper(all_genes[j]),"_GENE_USAGE_",groups[i],"_REPERTOIRE_",project_name,".png", sep = "")
        png(somePNGPath, width=5000, height=2500, units = "px", res = 300)
        print(p1)
        n <- 20
        dev.off()
      }
    }
  }
}

dir.create("OUTPUT/02.7GENE_USAGE_CORRELATIONS/")
for(j in 1:length(all_genes)){
  current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
  current_gu <- geneUsage(current, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
  current_gu <- current_gu[!is.na(current_gu$Names),]
  current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
  imm_gu_js <- geneUsageAnalysis(current_gu, .method = "js", .verbose = F)
  colnames(imm_gu_js) <- paste(meta[match(colnames(imm_gu_js), meta$Sample),"SUBTISSUE"],colnames(imm_gu_js), sep = "|")
  row.names(imm_gu_js) <- paste(meta[match(row.names(imm_gu_js), meta$Sample),"SUBTISSUE"],row.names(imm_gu_js), sep = "|")
  row.names(imm_gu_js) <- gsub("\\.",":",row.names(imm_gu_js))
  imm_gu_cor <- geneUsageAnalysis(current_gu, .method = "cor", .verbose = F)
  colnames(imm_gu_cor) <- paste(meta[match(colnames(imm_gu_cor), meta$Sample),"SUBTISSUE"],colnames(imm_gu_cor), sep = "|")
  row.names(imm_gu_cor) <- paste(meta[match(row.names(imm_gu_cor), meta$Sample),"SUBTISSUE"],row.names(imm_gu_cor), sep = "|")
  row.names(imm_gu_cor) <- gsub("\\.",":",row.names(imm_gu_cor))
  
  p1 <- vis(imm_gu_js, .title = paste("Gene usage JS-divergence: ", toupper(all_genes[j]), sep = ""), .leg.title = "JS", .text.size = 1.5)+
    # scale_fill_gradientn(colours = jet2.col(n = 100, alpha = 1))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
  p2 <- vis(imm_gu_cor, .title = paste("Gene usage correlation: ", toupper(all_genes[j]), sep = ""), .leg.title = "Cor", .text.size = 1.5)+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
  p1$layers[[2]] <- NULL
  p2$layers[[2]] <- NULL
  
  somePNGPath = paste(output_dir,"02.7GENE_USAGE_CORRELATIONS/02.7SCA_PLOT_",toupper(all_genes[j]),"_GENE_USAGE_JS-DIVERGENCE_CORRELATION_BCR_",project_name,".png", sep = "")
  png(somePNGPath, width=7000, height=3500, units = "px", res = 300)
  print(p1 + p2 + plot_annotation(title = paste(toupper(all_genes[j]), ": CORRELATION ANALYSIS FOR BCR REPERTOIRES", sep = ""), 
  theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
  
  current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
  current_gu <- geneUsage(current, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
  current_gu <- current_gu[!is.na(current_gu$Names),]
  current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
  imm_gu_js <- geneUsageAnalysis(current_gu, .method = "js", .verbose = F)
  colnames(imm_gu_js) <- paste(meta[match(colnames(imm_gu_js), meta$Sample),"SUBTISSUE"],colnames(imm_gu_js), sep = "|")
  row.names(imm_gu_js) <- paste(meta[match(row.names(imm_gu_js), meta$Sample),"SUBTISSUE"],row.names(imm_gu_js), sep = "|")
  row.names(imm_gu_js) <- gsub("\\.",":",row.names(imm_gu_js))
  imm_gu_cor <- geneUsageAnalysis(current_gu, .method = "cor", .verbose = F)
  colnames(imm_gu_cor) <- paste(meta[match(colnames(imm_gu_cor), meta$Sample),"SUBTISSUE"],colnames(imm_gu_cor), sep = "|")
  row.names(imm_gu_cor) <- paste(meta[match(row.names(imm_gu_cor), meta$Sample),"SUBTISSUE"],row.names(imm_gu_cor), sep = "|")
  row.names(imm_gu_cor) <- gsub("\\.",":",row.names(imm_gu_cor))
  
  p1 <- vis(imm_gu_js, .title = paste("Gene usage JS-divergence: ", toupper(all_genes[j]), sep = ""), .leg.title = "JS", .text.size = 1.5)+
    # scale_fill_gradientn(colours = jet2.col(n = 100, alpha = 1))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
  p2 <- vis(imm_gu_cor, .title = paste("Gene usage correlation: ", toupper(all_genes[j]), sep = ""), .leg.title = "Cor", .text.size = 1.5)+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
  p1$layers[[2]] <- NULL
  p2$layers[[2]] <- NULL
  
somePNGPath = paste(output_dir,"02.7GENE_USAGE_CORRELATIONS/02.7SCA_PLOT_",toupper(all_genes[j]),"_GENE_USAGE_JS-DIVERGENCE_CORRELATION_TCR_",project_name,".png", sep = "")
  png(somePNGPath, width=8500, height=4200, units = "px", res = 300)
  print(p1 + p2 + plot_annotation(title = paste(toupper(all_genes[j]), ": CORRELATION ANALYSIS FOR TCR REPERTOIRES", sep = ""),theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
  
}

# somePDFPath = paste(output_dir,"02.8SCA_PLOT_HCLUSTERING_K-MEANS_ALL_SAMPLES_ANALYSIS_",project_name,".pdf", sep = "")
# pdf(file=somePDFPath, width=14, height=6,pointsize=10)
# # for(i in 1:length(cell_types)){
# # for(j in 1:length(all_genes)){
# current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[i]),"Sample"])]
#   current_gu <- geneUsage(current, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
#   current_gu <- current_gu[!is.na(current_gu$Names),]
#   current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
#   colnames(current_gu) <- gsub("(.*?_S.*)_.*","\\1",colnames(current_gu))
#   print(vis(geneUsageAnalysis(current_gu, "cosine+hclust", .verbose = F))+ggtitle(paste("Optimal number of clusters\n(Clustered based on: ", toupper(all_genes[j]), ")", sep = "")))
# # }
# # }
# dev.off()

dir.create("OUTPUT/02.8SAMPLE_DISTANCE/")
cell_types <- sort(unique(meta$CELL_TYPE))
for(i in 1:length(cell_types)){
  # somePDFPath = paste(output_dir,"02.9SAMPLE_DISTANCE/02.9SCA_PLOT_GENE_WISE_SAMPLE_DISTANCE_PCA_MDS_",cell_types[i],"_REPERTOIRES_",project_name,".pdf", sep = "")
  # pdf(file=somePDFPath, width=30, height=10,pointsize=10)
for(j in 1:length(all_genes)){
  somePNGPath = paste(output_dir,"02.8SAMPLE_DISTANCE/02.8SCA_PLOT_GENE_",all_genes[j],"_SAMPLE_DISTANCE_PCA_MDS_",cell_types[i],"_REPERTOIRES_",project_name,".png", sep = "")
  png(somePNGPath, width=6000, height=2000, units = "px", res = 300)
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[i]),"Sample"])]
    current_gu <- geneUsage(current, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
    current_gu <- current_gu[!is.na(current_gu$Names),]
    current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
    imm_cl_pca <- geneUsageAnalysis(current_gu, "js+pca+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
    imm_cl_mds <- geneUsageAnalysis(current_gu, "js+mds+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
    imm_cl_tsne <- geneUsageAnalysis(current_gu, "js+tsne+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .perp = .01, .verbose = F)
    p1 <- vis(imm_cl_pca, .plot = "clust") + ggtitle(paste("PCA: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+MDS+K-MEANS)", sep = "")) + scale_color_manual(values = color_conditions$cold) +
      scale_fill_manual(values = color_conditions$cold) + xlab("PC1")+ylab("PC2")+guides(color = guide_legend(title="CLUSTERS"))
    p2 <- vis(imm_cl_mds, .plot = "clust") + ggtitle(paste("MDS: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+MDS+K-MEANS)", sep = "")) + scale_color_manual(values = color_conditions$cold) +
      scale_fill_manual(values = color_conditions$cold) + xlab("MDS1")+ylab("MDS2")+guides(color = guide_legend(title="CLUSTERS"))
    p3 <- vis(imm_cl_tsne, .plot = "clust") + ggtitle(paste("tSNE: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+tSNE+K-MEANS)", sep = "")) + scale_color_manual(values = color_conditions$cold) +
      scale_fill_manual(values = color_conditions$cold) + xlab("tSNE1")+ylab("tSNE2")+guides(color = guide_legend(title="CLUSTERS"))
    print(p1+p2+p3+plot_annotation(title = paste(toupper(all_genes[j]), ": SAMPLE DISTANCE FOR ",cell_types[i]," REPERTOIRES", sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))

  dev.off()
}
}
  
for(i in 1:length(cell_types)){
  somePNGPath = paste(output_dir,"02.9SCA_PLOT_ALL_GENE_SAMPLE_DISTANCE_PCA_MDS_",cell_types[i],"_REPERTOIRES_",project_name,".png", sep = "")
  png(somePNGPath, width=6000, height=2000, units = "px", res = 300)
  current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[i]),"Sample"])]
  current_gu <- geneUsage(current, .norm = T)
  current_gu <- current_gu[!is.na(current_gu$Names),]
  current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
  imm_cl_pca <- geneUsageAnalysis(current_gu, "js+pca+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
  imm_cl_mds <- geneUsageAnalysis(current_gu, "js+mds+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
  imm_cl_tsne <- geneUsageAnalysis(current_gu, "js+tsne+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .perp = .01, .verbose = F)
  p1 <- vis(imm_cl_pca, .plot = "clust") + ggtitle(paste("ALL GENES PCA: JS-DIVERGENCE+MDS+K-MEANS", sep = "")) + scale_color_manual(values = color_conditions$cold) +
    scale_fill_manual(values = color_conditions$cold) + xlab("PC1")+ylab("PC2")+guides(color = guide_legend(title="CLUSTERS"))
  p2 <- vis(imm_cl_mds, .plot = "clust") + ggtitle(paste("ALL GENES MDS: JS-DIVERGENCE+MDS+K-MEANS", sep = "")) + scale_color_manual(values = color_conditions$cold) +
    scale_fill_manual(values = color_conditions$cold) + xlab("MDS1")+ylab("MDS2")+guides(color = guide_legend(title="CLUSTERS"))
  p3 <- vis(imm_cl_tsne, .plot = "clust") + ggtitle(paste("ALL GENES tSNE: JS-DIVERGENCE+tSNE+K-MEANS", sep = "")) + scale_color_manual(values = color_conditions$cold) +
    scale_fill_manual(values = color_conditions$cold) + xlab("tSNE1")+ylab("tSNE2")+guides(color = guide_legend(title="CLUSTERS"))
  print(p1+p2+p3+plot_annotation(title = paste(toupper(all_genes[j]), ": SAMPLE DISTANCE FOR ",cell_types[i]," REPERTOIRES", sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  
  dev.off()
}

dir.create("OUTPUT/03.0SPECTRATYPE_SUBTISSUES/")
subtissues <- sort(unique(meta$SUBTISSUE))
for(i in 1:length(subtissues)){
  for(j in 1:length(cell_types)){
    if(length(which((names(contig_monarch$data) %in% meta[which(meta$SUBTISSUE == subtissues[i]),"Sample"]) & (names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[j]),"Sample"]))) > 0){
    somePNGPath = paste(output_dir,"/03.0SPECTRATYPE_SUBTISSUES/03.0SCA_PLOT_",cell_types[j],"_SPECTRATYPING_",subtissues[i],"_",project_name,".png", sep = "")
    png(somePNGPath, width=5500, height=3000, units = "px", res = 300)
  current <- contig_monarch$data[which((names(contig_monarch$data) %in% meta[which(meta$SUBTISSUE == subtissues[i]),"Sample"]) & (names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[j]),"Sample"]))]
  current <- do.call(rbind.data.frame, current)
  p1 <- vis(spectratype(current, .quant = "id", .col = "nt"))+
    xlab('CDR3 Nucleotide Length')+ ggtitle("")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(vjust = 0.5, hjust=1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
  
  if(nrow(current) > 10){
  p2 <- vis(spectratype(current,.col = "aa+v")) + 
    xlab('CDR3 AA+V Length')+ ggtitle("")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(vjust = 0.5, hjust=1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face = 2, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")
p <- p1+p2
rm(p1)
rm(p2)

  }else{
    p <- p1
    rm(p1)
    
  }
  
  print(p + plot_annotation(title = paste(cell_types[j], " REPERTOIRE SPECTRATYPING\n\n", subtissues[i], sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  dev.off()
  rm(p)
    }
}
}

# Diversity estimation
diversity_estimates <- NULL
diversity_methods <- c("chao1", "hill", "div", "gini.simp", "inv.simp", "d50","raref")

# BCR
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
diversity_estimates_bcr <- NULL
for(i in 1:length(diversity_methods)){
  if(diversity_methods[i] == "raref"){
    temp <- current[which(as.numeric(unlist(lapply(current,function(x){nrow(x)}))) > 50)]
  }else{
    temp <- current
  }
  diversity_estimates_bcr[[i]] <- repDiversity(temp, .method = diversity_methods[i])
}

for(i in 1:length(diversity_estimates_bcr)){
  write.csv(diversity_estimates_bcr[[i]], paste(output_dir,"/BCR_",diversity_methods[i],"_SAMPLEWISE_RESULTS.csv", sep = ""), quote=F, row.names = F)
}

# TCR
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
diversity_estimates_tcr <- NULL
for(i in 1:length(diversity_methods)){
  if(diversity_methods[i] == "raref"){
    temp <- current[which(as.numeric(unlist(lapply(current,function(x){nrow(x)}))) > 50)]
  }else{
    temp <- current
  }
  diversity_estimates_tcr[[i]] <- repDiversity(temp, .method = diversity_methods[i])
}

for(i in 1:length(diversity_estimates_tcr)){
  write.csv(diversity_estimates_tcr[[i]], paste(output_dir,"/TCR_",diversity_methods[i],"_SAMPLEWISE_RESULTS.csv", sep = ""), quote=F, row.names = F)
}


# BCR - TISSUEWISE
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}

names(current_merged) <- current_tissues


diversity_estimates_bcr <- NULL
for(i in 1:length(diversity_methods)){
  if(diversity_methods[i] == "raref"){
    temp <- current_merged[which(as.numeric(unlist(lapply(current_merged,function(x){nrow(x)}))) > 50)]
  }else{
    temp <- current_merged
  }
  diversity_estimates_bcr[[i]] <- repDiversity(temp, .method = diversity_methods[i])
}

for(i in 1:length(diversity_estimates_bcr)){
  write.csv(diversity_estimates_bcr[[i]], paste(output_dir,"/BCR_",diversity_methods[i],"_TISSUEWISE_RESULTS.csv", sep = ""), quote=F, row.names = F)
}

# TCR - TISSUEWISE
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}

names(current_merged) <- current_tissues

diversity_estimates_tcr <- NULL
for(i in 1:length(diversity_methods)){
  if(diversity_methods[i] == "raref"){
    temp <- current_merged[which(as.numeric(unlist(lapply(current_merged,function(x){nrow(x)}))) > 50)]
  }else{
    temp <- current_merged
  }
  diversity_estimates_tcr[[i]] <- repDiversity(temp, .method = diversity_methods[i])
}

for(i in 1:(length(diversity_estimates_tcr)-1)){
  write.csv(diversity_estimates_tcr[[i]], paste(output_dir,"/TCR_",diversity_methods[i],"_TISSUEWISE_RESULTS.csv", sep = ""), quote=F, row.names = F)
}


dir.create("OUTPUT/03.1DIVERSITY_ESTIMATES/")

for(i in 1:length(diversity_estimates_bcr)){
  somePNGPath = paste(output_dir,"/03.1DIVERSITY_ESTIMATES/03.1SCA_PLOT_BCR_REPETOIRE_DIVERSITY_ESTIMATES_",toupper(diversity_methods[i]),"_",project_name,".png", sep = "")
  png(somePNGPath, width=5000, height=3000, units = "px", res = 250)
  p1 <- vis(diversity_estimates_bcr[[i]], .by = c("SUBTISSUE"), .meta = contig_monarch$meta, .test=F)+
    scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) + 
    scale_color_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 25),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1))
  print(p1+plot_annotation(title = "BCR REPETROIRE DIVERSITY ESTIMATION", 
                              theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
}

for(i in 1:length(diversity_estimates_tcr)){
  somePNGPath = paste(output_dir,"/03.1DIVERSITY_ESTIMATES/03.1SCA_PLOT_TCR_REPETOIRE_DIVERSITY_ESTIMATES_",toupper(diversity_methods[i]),"_",project_name,".png", sep = "")
  png(somePNGPath, width=5000, height=3000, units = "px", res = 250)
  p1 <- vis(diversity_estimates_tcr[[i]], .by = c("SUBTISSUE"), .meta = contig_monarch$meta, .test=F)+
    scale_fill_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) + 
    scale_color_manual(values = gen_colors(color_conditions$general, length(unique(meta$SUBTISSUE)))) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 25),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "right")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))
  print(p1+plot_annotation(title = "TCR REPETROIRE DIVERSITY ESTIMATION", 
                           theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
}


dir.create("OUTPUT/03.2SCA_PLOT_ALLUVIAL/")
# TCR - TISSUEWISE
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}

names(current_merged) <- current_tissues

n <- 10
for(i in 1:length(current_merged)){
  tc1 <- trackClonotypes(current_merged, list(i, n), .col = "nt")
  tc2 <- trackClonotypes(current_merged, list(i, n), .col = "aa+v")
  p1 <- vis(tc1)+ggtitle(paste("TCR TOP ",n," CLONOTYPES (CDR3 NUCLEOTIDE) IN ", names(current_merged)[i], sep = ""))+
    scale_fill_manual(values = gen_colors(color_conditions$tenx, n))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, "cm"),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 30, face=2),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+xlab("SUB-TISSUES")+ylab("PROPORTION")
  
  p2 <- vis(tc2)+ggtitle(paste("TCR TOP ",n," CLONOTYPES (CDR3 AA+V) IN ", names(current_merged)[i], sep = ""))+
    scale_fill_manual(values = gen_colors(color_conditions$tenx, n))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, "cm"),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 30, face=2),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 2))+xlab("SUB-TISSUES")+ylab("PROPORTION")

  somePNGPath = paste(output_dir,"03.2SCA_PLOT_ALLUVIAL/03.2SCA_PLOT_TCR_ALLUVIAL_",names(current_merged)[i],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_NUCLEOTIDES_",project_name,".png", sep = "")
  png(somePNGPath, width=6500, height=5000, units = "px", res = 300)
  print(p1)
  dev.off()
  
  somePNGPath = paste(output_dir,"03.2SCA_PLOT_ALLUVIAL/03.2SCA_PLOT_TCR_ALLUVIAL_",names(current_merged)[i],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_CDR3_AA+V_",project_name,".png", sep = "")
  png(somePNGPath, width=5500, height=3000, units = "px", res = 200)
  print(p2)
  dev.off()
}


# BCR - TISSUEWISE
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "BCR"),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}

names(current_merged) <- current_tissues

n <- 10
for(i in 1:length(current_merged)){
  tc1 <- trackClonotypes(current_merged, list(i, n), .col = "nt")
  tc2 <- trackClonotypes(current_merged, list(i, n), .col = "aa+v")
  p1 <- vis(tc1)+ggtitle(paste("BCR TOP ",n," CLONOTYPES (CDR3 NUCLEOTIDE) IN ", names(current_merged)[i], sep = ""))+
    scale_fill_manual(values = gen_colors(color_conditions$tenx, n))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, "cm"),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 30, face=2),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+xlab("SUB-TISSUES")+ylab("PROPORTION")
  
  p2 <- vis(tc2)+ggtitle(paste("BCR TOP ",n," CLONOTYPES (CDR3 AA+V) IN ", names(current_merged)[i], sep = ""))+
    scale_fill_manual(values = gen_colors(color_conditions$tenx, n))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, "cm"),
          strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 30, face=2),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom")+guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 2))+xlab("SUB-TISSUES")+ylab("PROPORTION")
  
  somePNGPath = paste(output_dir,"03.2SCA_PLOT_ALLUVIAL/03.2SCA_PLOT_BCR_ALLUVIAL_",names(current_merged)[i],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_NUCLEOTIDES_",project_name,".png", sep = "")
  png(somePNGPath, width=6500, height=5000, units = "px", res = 300)
  print(p1)
  dev.off()
  
  somePNGPath = paste(output_dir,"03.2SCA_PLOT_ALLUVIAL/03.2SCA_PLOT_BCR_ALLUVIAL_",names(current_merged)[i],"_TOP_",n,"_MOST_ABUNDANT_CLONOTYPES_CDR3_AA+V_",project_name,".png", sep = "")
  png(somePNGPath, width=5500, height=3000, units = "px", res = 200)
  print(p2)
  dev.off()
}


# DB
# vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")
# vdjdb
# vdjdb_st = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/SearchTable-2019-10-17%2012_36_11.989.tsv.gz", "vdjdb-search", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
# vdjdb_st
# tbadb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "CMV")
# tbadb

vdjdb <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = "HomoSapiens")
# mcpas <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human")

result_vdjdb <- dbAnnotate(contig_monarch$data, vdjdb, "CDR3.aa", "cdr3")
write.csv(result_vdjdb, paste(output_dir,"03.4SCA_TABLE_MATCHED_WITH_VDJDB_CDR3AA_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = F)

plotx <- result_vdjdb
plotx <- data.frame(plotx)
plotx <- reshape2::melt(plotx)
plotx <- plotx[which(!plotx$value == 0),]
plotx <- plotx[which(!plotx$variable == "Samples"),]
colnames(plotx) <- c("CDR3.aa","SAMPLES","COUNT")
plotx$CELL_TYPE <- gsub(".*_(.*)","\\1",plotx$SAMPLES)
plotx$SUBTISSUE <- meta[match(plotx$SAMPLES, meta$Sample),"SUBTISSUE"]
plotx <- plyr::count(plotx, vars = c("SUBTISSUE","CDR3.aa"), wt_var = "COUNT")
colnames(plotx) <- c("SUBTISSUE","CDR3.aa","COUNT")
plotx <- cbind(plotx, vdjdb[match(plotx$CDR3.aa,vdjdb$cdr3),])
plotx <- plotx[order(plotx$COUNT, decreasing = T),]
write.csv(plotx, paste(output_dir,"03.4SCA_TABLE_SUBTISSUES_VDJDB_CDR3AA_TCR_",project_name,".csv", sep = ""), quote = F, row.names = F)

write.table(plotx[which(plotx$COUNT > 5),], paste(output_dir,"03.4SCA_TABLE_FILTERED_SUBTISSUES_VDJDB_CDR3AA_TCR_",project_name,".csv", sep = ""), quote = F, row.names = F, sep = ";")

plotx <- plotx[which(plotx$COUNT > 5),]
ggplot(plotx, aes(CDR3.aa,COUNT, fill = Chain))+
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~Pathology+SUBTISSUE)+
  scale_fill_manual(values=gen_colors(color_conditions$tenx, length(unique(plotx$CDR3.aa))))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
 


# library(gplots)
# 
# somePNGPath = paste(output_dir,"scIMMUNE_VDJDB_COUNT_HEATMAP.png", sep = "")
# png(somePNGPath, width=4000, height=3100, units = "px", res = 300)
# heatmap.2(as.matrix(t(plotx)),margin=c(15, 15), trace="none",key=T, keysize=1,
#           col=jet2.col (n = 100, alpha = 1),srtCol=45, scale="none", Colv = T,
#           Rowv = T,
#           density.info="none", cexCol=0.8,cexRow=0.8)
# # corrplot(corr_medians, col = rev(col2(100)), order = "hclust", tl.col = "black", tl.cex = 2, type="upper",tl.pos="tp")
# dev.off()


# resulplotxt_mcpas <- dbAnnotate(contig_monarch$data, mcpas, c("CDR3.aa", "V.name"), c("CDR3.beta.aa", "TRBV"))
# write.csv(number_contig, paste(output_dir,"01.0SCA_TABLE_PERCENT_UNIQUE_CLONOTYPE_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)

# K-mer statistics
kmer_length <- 5
n <- 10

kmers <- getKmers(contig_monarch$data, kmer_length)
kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
write.csv(kmers, paste(output_dir,"03.5SCA_TABLE_ALL_AA_KMER-",kmer_length,"_ALL_SAMPLES_",project_name,".csv", sep = ""), quote = F, row.names = F)


#### TCR SUBTISSUES ####
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == "TCR"),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}
names(current_merged) <- current_tissues

kmers_tcr <- getKmers(current_merged, kmer_length)
kmers_tcr <- kmers[grep("\\;", kmers_tcr$Kmer, ignore.case = T, invert = T),]

# somePNGPath = paste(output_dir,"03.5SCA_PLOT_TOP_",n,"_AA_KMER-",kmer_length,"_ALL_SAMPLES_ANALYSIS_",project_name,".png", sep = "")
# png(somePNGPath, width=5500, height=3000, units = "px", res = 300)
# 
# p1 <- vis(kmers_tcr, .position = "stack", .head = n)
# p2 <- vis(kmers_tcr, .head = n, .position = "dodge")
# print(p1/p2)
# 
# dev.off()

# Sequence motif analysis
dir.create("OUTPUT/03.6_KMER_MOTIFS/")
kmer_length <- 9

for(k in 1:length(cell_types)){
current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[k]),"Sample"])]
current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
current_tissues <- unique(current_labels)

current_merged <- NULL
for(i in 1:length(current_tissues)){
  temp <- current[which(current_labels == current_tissues[i])]
  temp <- do.call(rbind.data.frame, temp)
  current_merged[[i]] <- temp
}
names(current_merged) <- current_tissues

for(i in 1:length(current_merged)){
  
  kmers <- getKmers(current_merged[[i]], kmer_length)
  kmers_aa_stats <- kmer_profile(kmers)
  colnames(kmers_aa_stats) <- paste("POS_", 1:ncol(kmers_aa_stats), sep = "")
  write.csv(kmers_aa_stats, paste(output_dir,"03.6_KMER_MOTIFS/03.6SCA_TABLE_",cell_types[k],"_",names(current_merged)[i],"_AA_SEQUENCE_KMER-",kmer_length,"_MOTIFS_COUNTS_",project_name,".csv", sep = ""), quote = F, row.names = T)
  
  somePNGPath = paste(output_dir,"03.6_KMER_MOTIFS/03.7SCA_PLOT_",cell_types[k],"_",names(current_merged)[i],"_TOP_",n,"_KMER-",kmer_length,"_SEQUENCE_MOTIFS_",project_name,".png", sep = "")
  png(somePNGPath, width=6000, height=3000, units = "px", res = 300)
  p1 <- vis(kmers_aa_stats) +
    scale_color_manual(values = gen_colors(color_conditions$general,nrow(kmers_aa_stats)))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 10),
          # legend.key.size = unit(1, "cm"),
          # strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face=1, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"))+
    guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+
    xlab("KMER POSITION")+ylab("PROPORTION")
  
  p2 <- vis(kmers_aa_stats, .plot = "seq")+scale_fill_manual(values = gen_colors(color_conditions$general,nrow(kmers_aa_stats)))+
    theme(plot.margin = unit(c(2,2,2,2), "cm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          # legend.key.size = unit(1, "cm"),
          # strip.text.x = element_text(size = 15),
          plot.title = element_text(size = 20, face=1, hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"))+
    guides(color=guide_legend(title="CHEMISTRY", ncol = 5),fill = guide_legend(title="CHEMISTRY", ncol = 5))+xlab("KMER POSITIONS")+ylab("PROPORTION")
  
  print(p1+p2+plot_annotation(title = paste(cell_types[k], " K-MER PROFILE (K = ",kmer_length,") - ", names(current_merged)[i], sep = ""),theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
  dev.off()
}
}

for(k in 1:length(cell_types)){
  current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[k]),"Sample"])]
  current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
  current_tissues <- unique(current_labels)
  
  current_merged <- NULL
  for(i in 1:length(current_tissues)){
    temp <- current[which(current_labels == current_tissues[i])]
    temp <- do.call(rbind.data.frame, temp)
    current_merged[[i]] <- temp
  }
  names(current_merged) <- current_tissues
  
ov <- repOverlap(current_merged)
grid.col = gen_colors(color_conditions$tableau20, length(current_merged))
names(grid.col) <- gsub(".*\\:(.*)","\\1",names(current_merged))
colnames(ov) <- gsub(".*\\:(.*)","\\1",colnames(ov))
row.names(ov) <- gsub(".*\\:(.*)","\\1",row.names(ov))
somePNGPath = paste(output_dir, "03.7SCA_PLOT_CIRCOS_DIAGRAM_REPERTOIRE_OVERLAPS_",cell_types[k],"_",project_name,".png", sep = "")
png(somePNGPath, width=5500, height=5000, units = "px", res = 300)
vis(ov, .plot = "circos", annotationTrack = c("grid", "axis"),preAllocateTracks = 1, grid.col = grid.col)
# title(paste(cell_types[k]," REPERTOIRE OVERLAPS: SUBTISSUE-LEVEL", sep = ""), cex = 4)
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 1.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
dev.off()

somePNGPath = paste(output_dir, "03.8SCA_PLOT_HEATMAP_REPERTOIRE_OVERLAPS_",cell_types[k],"_",project_name,".png", sep = "")
png(somePNGPath, width=4000, height=3000, units = "px", res = 300)
heatmap.2(as.matrix(t(ov)),
          margin=c(15, 15), main = paste(cell_types[k], " REPERTOIRE OVERLAPS"),
          trace="none",key=T, keysize=1,col=jet2.col (n = 100, alpha = 1),srtCol=45, scale="none", Colv = T,Rowv = T,density.info="none", cexCol=1.0,cexRow=1.0)
dev.off()

}

################################################################################################################
############################ By BCR / TCR separately ###########################################################
# contig_files <- contig_files[grep("NUM_94_", contig_files, ignore.case = T, invert = T)]
# pheno_data <- pheno_data[grep("NUM_94_", pheno_data$FILE, ignore.case = T, invert = T),]
contig_files <- paste("CONTIGS/",meta$FILE, sep = "")
contig_list <- lapply(contig_files, function(x){x <- readfile(x)})
names(contig_list) <- meta$Sample

contig_underscore <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))}))
contig_hyphen <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("-", x$barcode))}))

if(length(which(contig_underscore == 0))>1 & length(which(contig_hyphen > 0)) == length(contig_list)){
  contig_list <- contig_list
} else{
  
  for (i in seq_along(contig_list)) {
    contig_list[[i]] <- stripBarcode(contig_list[[i]], column = 1, connector = "_",
                                     num_connects = max(as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))})))+1)
  }
  
}

summary <- lapply(contig_list, function(x){rbind(nrow(x))})

bcr_list <- meta[meta$CELL_TYPE == "BCR","Sample"]
tcr_list <- meta[meta$CELL_TYPE == "TCR","Sample"]

for(i in 1:length(contig_list)){
  if(!is.logical(contig_list[[i]]$is_cell)){
    contig_list[[i]]$is_cell <- ifelse(toupper(contig_list[[i]]$is_cell) == "TRUE", TRUE, FALSE)
  }
  
  if(!is.logical(contig_list[[i]]$high_confidence)){
    contig_list[[i]]$high_confidence <- ifelse(toupper(contig_list[[i]]$high_confidence) == "TRUE", TRUE, FALSE)
  }
  
  if(!is.logical(contig_list[[i]]$full_length)){
    contig_list[[i]]$full_length <- ifelse(toupper(contig_list[[i]]$full_length) == "TRUE", TRUE, FALSE)
  }
  
  if(!is.logical(contig_list[[i]]$productive)){
    contig_list[[i]]$productive <- ifelse(toupper(contig_list[[i]]$productive) == "TRUE", TRUE, FALSE)
  }
}

# contig_bcr <- combineBCR(contig_list[c(15,31)], 
#                          samples = meta[match(bcr_list,meta$Sample),"Sample"][c(15,31)],
#                          ID = meta[match(bcr_list,meta$Sample),"SUBTISSUE"][c(15,31)])

contig_bcr <- combineBCR(contig_list[which(names(contig_list) %in% bcr_list)], 
                           samples = meta[match(bcr_list,meta$Sample),"Sample"],
                           removeNA = TRUE,
                           ID = meta[match(bcr_list,meta$Sample),"SUBTISSUE"])

tab_presence <- NULL
tgd_presence <- NULL
contig_tcr_tab <- NULL
contig_tcr_tgd <- NULL

  contig_tcr <- contig_list[which(names(contig_list) %in% tcr_list)]
  tab_presence <- as.character(unlist(lapply(contig_tcr, function(x){x <- ifelse(length(grep("TRA|TRB",unique(x$chain), ignore.case = T)) > 0, "TRUE")})))
  if(length(grep("TRUE", tab_presence, ignore.case = T)) > 0){
    contig_tcr_tab <- combineTCR(contig_tcr[which(tab_presence == "TRUE")], 
                                 samples = meta[match(names(contig_tcr[which(tab_presence == "TRUE")]),meta$Sample),"Sample"],
                                 removeNA = TRUE,
                                 ID = meta[match(names(contig_tcr[which(tab_presence == "TRUE")]),meta$Sample),"SUBTISSUE"],
                                 cells = "T-AB")
  }
  

bcr_subtissues <- sort(unique(meta[meta$CELL_TYPE == "BCR","SUBTISSUE"]))
subtissues_contig_bcr <- NULL
names(contig_bcr) <- gsub("(.*BCR)_.*","\\1",names(contig_bcr))
for(j in 1:length(bcr_subtissues)){
  current <- contig_bcr[which(names(contig_bcr) %in% meta[which(meta$CELL_TYPE == "BCR" & meta$SUBTISSUE == bcr_subtissues[j]),"Sample"])]
  current <- do.call(rbind.data.frame, current)
  subtissues_contig_bcr[[j]] <- current
}

names(subtissues_contig_bcr) <- bcr_subtissues


tcr_subtissues <- sort(unique(meta[meta$CELL_TYPE == "TCR","SUBTISSUE"]))
subtissues_contig_tcr <- NULL
names(contig_tcr_tab) <- gsub("(.*TCR)_.*","\\1",names(contig_tcr_tab))
for(j in 1:length(tcr_subtissues)){
  current <- contig_tcr_tab[which(names(contig_tcr_tab) %in% meta[which(meta$CELL_TYPE == "TCR" & meta$SUBTISSUE == tcr_subtissues[j]),"Sample"])]
  current <- do.call(rbind.data.frame, current)
  subtissues_contig_tcr[[j]] <- current
}

names(subtissues_contig_tcr) <- tcr_subtissues

contig_pheno <- meta
contig_types <- c("BCR","TCR_AB")
contig_table <- NULL
contig_table[[1]] <- subtissues_contig_bcr
contig_table[[2]] <- subtissues_contig_tcr

for(i in 1:length(contig_table)){
  
  current_type <- contig_types[i]
  
  # if(length(contig_table[[i]]) > 0){
  #   contig_table[[i]] <- lapply(contig_table[[i]], function(x){
  #     data.frame(x,
  #                SAMPLE_ID = contig_pheno[match(unique(x$sample), contig_pheno$Sample),"Sample"],
  #                TISSUE = contig_pheno[match(unique(x$sample), contig_pheno$Sample),"TISSUE"],
  #                SUBTISSUE = contig_pheno[match(unique(x$sample), contig_pheno$Sample),"SUBTISSUE"],
  #                CELL_TYPE = contig_pheno[match(unique(x$sample), contig_pheno$Sample),"CELL_TYPE"])
  #   })
    
    # number_contig <- quantContig(contig_table[[i]], cloneCall="gene+nt", scale = T, exportTable = T)
    # write.csv(number_contig, paste(output_dir,"03.8SCA_TABLE_PERCENT_UNIQUE_CLONOTYPE_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    
   #  somePNGPath = paste(output_dir,"03.9SCA_PLOT_CLONOTYPE_ABUNDANCE_",current_type,"_",project_name,".png", sep = "")
   # png(somePNGPath, width=5500, height=5000, units = "px", res = 300)
   #  print(abundanceContig(contig_table[[i]], cloneCall = "gene", scale = F) + 
   #          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom", 
   #                plot.title = element_text(size = 30, face = "bold")) +
   #          guides(fill = guide_legend(ncol = 3)))
   #  dev.off()
    
    # abundance_contig <- abundanceContig(contig_table[[i]], cloneCall = "gene", exportTable = T)
    # write.csv(abundance_contig, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_GENE_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    
    # abundance_contig <- abundanceContig(contig_table[[i]], cloneCall = "aa", exportTable = T)
    # write.csv(abundance_contig, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_AA_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    current_single <- lengthContig(contig_table[[i]], cloneCall="aa", chains = "single", exportTable = TRUE)
    write.csv(current_single, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_AA_SEPARATE_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    
    current_combined <- lengthContig(contig_table[[i]], cloneCall="aa", chains = "combined", exportTable = TRUE)
    write.csv(current_combined, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_AA_COMBINED_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    
    current_single <- lengthContig(contig_table[[i]], cloneCall="nt", chains = "single", exportTable = TRUE)
    write.csv(current_single, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_NUCLEOTIDES_SEPARATE_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
    
    current_combined <- lengthContig(contig_table[[i]], cloneCall="nt", chains = "combined", exportTable = TRUE)
    write.csv(current_combined, paste(output_dir,"04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_NUCLEOTIDES_COMBINED_",current_type,"_",project_name,".csv", sep = ""), quote = F, row.names = F)
}

######### SUBTISSUES CLONOTYPE FREQUENCY #########################################################

clonotypes_freq <- NULL
current_files <- list.files(path = output_dir, pattern = "^04.0SCA_TABLE_CLONOTYPE_ABUNDANCE_.*csv", full.names = T)
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
    colnames(current) <- c("LENGTH","SEQUENCE","SUBTISSUE")
  }else{
    colnames(current) <- c("LENGTH","SEQUENCE","SUBTISSUE","CHAIN")
  }
  
  current_datatype <- gsub(".*TABLE_CLONOTYPE_ABUNDANCE_(.*?)_.*","\\1",current_files[i])
  current$LENGTH <- nchar(gsub("_","",current$SEQUENCE))
  # temp <- data.frame(table(current[,c("LENGTH","SUBTISSUE")]))
  # temp <- temp[!temp$Freq == 0,]
  current <- split(current, current$SUBTISSUE)
  for(j in 1:length(current)){
    temp <- current[[j]]
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
        CELL_TYPE = current_celltype,
        DATA_TYPE = current_datatype,
        CHAIN = "COMBINE",
        STATUS = ifelse(combined == T, "COMBINED","SINGLE"),
        MEDIAN_LENGTH = current_median)
    }else{
      temp2 <- cbind(
        SUBTISSUE = unique(temp$SUBTISSUE),
        CELL_TYPE = current_celltype,
        DATA_TYPE = current_datatype,
        STATUS = ifelse(combined == T, "COMBINED","SINGLE"),
        data.frame(unique(temp[,c("CHAIN","MEDIAN_LENGTH")])
        ))
    }
    tissues_kmer_summary <- rbind(tissues_kmer_summary,temp2)
    current[[j]] <- temp
  }
  current <- do.call(rbind.data.frame, current)
  write.csv(current, paste(output_dir,"/WITH_FREQUENCY_INFO_",gsub(".*\\/(.*)","\\1",current_files[i]), sep = ""), quote=F, row.names = F)
  write.csv(tissues_kmer_summary, paste(output_dir,"/SUBTISSUES_KMER_MEDIAN_LENGTH_SUMMARY.csv", sep = ""), quote=F, row.names = F)
  
  current$SUBTISSUE <- factor(current$SUBTISSUE, levels = sort(unique(current$SUBTISSUE)))
  current_file_name <- gsub(".csv",".png",current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=4000, units = "px", res = 300)
  p1 <- ggplot(current,aes(LENGTH,FREQUENCY,fill=SUBTISSUE))+ geom_bar(stat = "identity", position = "dodge")+ theme_classic()+ 
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
    print(p1+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), current_celltype, ": SUBTISSUE SPECIFIC CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
  current$SUBTISSUE <- factor(current$SUBTISSUE, levels = rev(sort(unique(current$SUBTISSUE))))
  current_file_name <- gsub(".csv","_BOXPLOT.png",current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=5000, units = "px", res = 300)
  p2 <- ggplot(current,aes(LENGTH,SUBTISSUE,fill=SUBTISSUE))+ geom_boxplot()+ theme_classic()+
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
  print(p2+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), current_celltype, ": SUBTISSUE SPECIFIC CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
  library(ggridges)
  current$SUBTISSUE <- factor(current$SUBTISSUE, levels = sort(unique(current$SUBTISSUE)))
  current_file_name <- gsub(".csv","_RIDGE.png",current_files[i])
  somePNGPath = current_file_name
  png(somePNGPath, width=5000, height=5000, units = "px", res = 300)
  p2 <- ggplot(current,aes(LENGTH,SUBTISSUE,fill=SUBTISSUE))+ geom_density_ridges()+ theme_classic()+
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
  print(p2+plot_annotation(title = paste(ifelse(combined==T,"COMBINED ",""), current_celltype, ": SUBTISSUE SPECIFIC CDR3 ",current_datatype," DISTRIBUTION", sep = ""), theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
  dev.off()
  
}

head(tissues_kmer_summary)
# Sequence motif analysis
dir.create("OUTPUT/04.1_KMER_MOTIFS_MEDIAN_LENGTH_SUBTISSUE_SPECIFIC/")
chains <- c("IGH","IGL","TRA","TRB")
for(k in 1:length(cell_types)){
  current <- contig_monarch$data[which(names(contig_monarch$data) %in% meta[which(meta$CELL_TYPE == cell_types[k]),"Sample"])]
  current_labels <- meta[match(names(current), meta$Sample),"SUBTISSUE"]
  current_tissues <- unique(current_labels)
  
  current_merged <- NULL
  for(i in 1:length(current_tissues)){
    temp <- current[which(current_labels == current_tissues[i])]
    temp <- do.call(rbind.data.frame, temp)
    current_merged[[i]] <- temp
  }
  names(current_merged) <- current_tissues
  
  for(i in 1:length(current_merged)){
    for(m in 1:length(chains)){
      kmer_length <- NULL
      kmer_length <- ceiling(tissues_kmer_summary[which(tissues_kmer_summary$SUBTISSUE == names(current_merged)[i] & tissues_kmer_summary$DATA_TYPE == "AA" & tissues_kmer_summary$CELL_TYPE == current_celltype & tissues_kmer_summary$CHAIN == chains[m]),"MEDIAN_LENGTH"])
      
      if(length(kmer_length) > 0){
        
        kmers <- getKmers(current_merged[[i]], kmer_length)
        kmers <- kmers[which(nchar(gsub(";","",kmers$Kmer)) == kmer_length),]
        for(j in 1:5){
          temp <- NULL
          temp <- getKmers(current_merged[[i]], kmer_length+j)
          temp <- temp[which(nchar(gsub(";","",temp$Kmer)) == kmer_length),]
          if(nrow(temp) > 0){
            kmers <- rbind(kmers, temp)
          }
        }
        
        kmers_aa_stats <- kmer_profile(kmers)
        colnames(kmers_aa_stats) <- paste("POS_", 1:ncol(kmers_aa_stats), sep = "")
        write.csv(kmers_aa_stats, paste(output_dir,"04.1_KMER_MOTIFS_MEDIAN_LENGTH_SUBTISSUE_SPECIFIC/04.1SCA_TABLE_CHAIN_",chains[m], "_", cell_types[k],"_",names(current_merged)[i],"_AA_SEQUENCE_KMER-",kmer_length,"_MOTIFS_COUNTS_",project_name,".csv", sep = ""), quote = F, row.names = T)
        
        colnames(kmers_aa_stats) <- gsub("POS_","", colnames(kmers_aa_stats))
        somePNGPath = paste(output_dir,"04.1_KMER_MOTIFS_MEDIAN_LENGTH_SUBTISSUE_SPECIFIC/04.1SCA_PLOT_",chains[m], "_",cell_types[k],"_",names(current_merged)[i],"_TOP_",n,"_KMER-",kmer_length,"_SEQUENCE_MOTIFS_",project_name,".png", sep = "")
        png(somePNGPath, width=6500, height=3000, units = "px", res = 300)
        p1 <- vis(kmers_aa_stats) +
          scale_color_manual(values = gen_colors(color_conditions$general,nrow(kmers_aa_stats)))+
          theme(plot.margin = unit(c(2,2,2,2), "cm"),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 10),
                # legend.key.size = unit(1, "cm"),
                # strip.text.x = element_text(size = 15),
                plot.title = element_text(size = 20, face=1, hjust = 0.5),
                strip.background = element_rect(colour="white", fill="white"))+
          guides(color=guide_legend(title="SUBTISSUE", ncol = 1),fill = guide_legend(ncol = 1))+
          xlab("KMER POSITION")+ylab("PROPORTION")
        
        p2 <- vis(kmers_aa_stats, .plot = "seq")+scale_fill_manual(values = gen_colors(color_conditions$general,nrow(kmers_aa_stats)))+
          theme(plot.margin = unit(c(2,2,2,2), "cm"),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                legend.title = element_text(size = 15),
                legend.text = element_text(size = 15),
                # legend.key.size = unit(1, "cm"),
                # strip.text.x = element_text(size = 15),
                plot.title = element_text(size = 20, face=1, hjust = 0.5),
                strip.background = element_rect(colour="white", fill="white"))+
          guides(color=guide_legend(title="CHEMISTRY", ncol = 5),fill = guide_legend(title="CHEMISTRY", ncol = 5))+xlab("KMER POSITIONS")+ylab("PROPORTION")
        
        print(p1+p2+plot_annotation(title = paste(cell_types[k], " ",chains[m]," K-MER PROFILE (K = ",kmer_length,") - ", names(current_merged)[i], sep = ""),theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
        dev.off()
      }
    }
  }
}





save.image(paste(output_dir,"/",project_name,".RData", sep = ""))
