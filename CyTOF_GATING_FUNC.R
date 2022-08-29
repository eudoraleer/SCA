###################################################################
##  Do not change anything unless you know the meaning!!!
##  Gating_fun.R Version 1.1.1
##  New Features in V1.1.1:
##  1. Acclimatizing the code to fit the output from Helios Machine
###################################################################

Gating_all <- function(input_fcs,barcode_file,barcode_name,cytof_machine,
                       Event_length_1=NULL,Event_length_2=NULL,
                       DNA_1=NULL,DNA_2=NULL)
{
  library('flowCore')
  library('quantmod')
  library('MASS')
  library('RColorBrewer')
  data=read.FCS(input_fcs,transformation=FALSE)
  file_list=NULL
  barcode_table=read.csv(barcode_file)
  marker_num=length(barcode_name)
  marker_name=barcode_name
  for(i in 1:(2^marker_num)) 
    file_list[i]=paste(barcode_table[i,dim(barcode_table)[2]],'_',barcode_table[i,1],'.fcs',sep='')
  x_name=parameters(data)$name
  if (length(grep("Helios",cytof_machine, ignore.case = T)) > 0){
    if (length(grep("_", parameters(data)$desc, invert = T)) > 0){
      Nas <- grep("_", parameters(data)$desc, invert = T)
      to_Nas <- Nas[!is.na(parameters(data)$desc[Nas])]
      parameters(data)$desc[to_Nas] <- NA
      
      for (k in 1:length(to_Nas)){
        desc_to_Nas <- which(description(data) %in% parameters(data)$desc[to_Nas])
        key_to_Nas <- which(keyword(data) %in% parameters(data)$desc[to_Nas])
        description(data)[desc_to_Nas[k]] <- NA
        keyword(data)[key_to_Nas[k]] <- NA
      }
    }
    parameters(data)$desc <- gsub("\\d+\\w+?_(.*)+_*.*","\\1",parameters(data)$desc)
    parameters(data)$desc <- gsub("_.*","",parameters(data)$desc)
    if(length(which(parameters(data)$desc == "Live")) > 0){
      parameters(data)$desc <- gsub("Live","DeadLive",parameters(data)$desc)
    }
    x_desc <- parameters(data)$desc
    
    for (i in grep("P[0-9]+",names(description(data)))){
      if (description(data)[i] != "Event_length"){
        description(data)[i] <- gsub("\\d+\\w+?_(.*)+_*.*","\\1",description(data)[i])
        description(data)[i] <- gsub("_.*","",description(data)[i])
        if (description(data)[i] == "Live"){
          description(data)[i] <- "DeadLive"
        }
      }
    }
    
    for (i in grep("P[0-9]+",names(keyword(data)))){
      if (keyword(data)[i] != "Event_length"){
        keyword(data)[i] <- gsub("\\d+\\w+?_(.*)+_*.*","\\1",keyword(data)[i])
        keyword(data)[i] <- gsub("_.*","",keyword(data)[i])
        if (keyword(data)[i] == "Live"){
          keyword(data)[i] <- "DeadLive"
        }
      }
    }
    
  } else if (length(grep("CyTOF",cytof_machine, ignore.case = T)) > 0){
    x_desc=parameters(data)$desc
  }
  marker_id=NULL
  for(i in 1:length(x_name))
  {
    if (x_name[i]=='Event_length') Event_length=i
    if (is.na(x_desc[i]) == F)
    {
      if (x_desc[i]=='DNA') DNA=i
      if (x_desc[i]=='DeadLive') DeadLive=i
      for (j in 1:marker_num) if (x_desc[i]==marker_name[j]) marker_id[j]=i
    }
  }
  
  ss=strsplit(input_fcs,'/')[[1]]
  plot_name=gsub('.fcs','',ss[length(ss)])
  pdf(paste(plot_name,'_gating_plot.pdf',sep=''),pointsize=30,width=50, height=20)
  par(mfrow=c(2,5))
  data_tmp=data
  xi=(exprs(data_tmp)[exprs(data_tmp)[,Event_length]>1 & exprs(data_tmp)[,Event_length]<150 &
                        !is.na(exprs(data_tmp)[,Event_length])
                      ,Event_length])
  #xi=xi[!is.na(xi)]
  plot(density(xi),main=parameters(data)$name[Event_length])
  cutoff=Find_one_pop(xi)
  if (is.null(Event_length_1)) Event_length_1=cutoff[1]
  if (is.null(Event_length_2)) Event_length_2=100
  abline(v=Find_one_pop(xi))
  abline(v=c(Event_length_1,Event_length_2),col='red')
  xi=log2(exprs(data_tmp)[exprs(data_tmp)[,DNA]>1,DNA])
  xi=xi[!is.na(xi)]
  plot(density(xi),main=parameters(data)$desc[DNA])
  abline(v=Find_one_pop(xi))
  cutoff=Find_one_pop(xi)
  if (is.null(DNA_1)) DNA_1=cutoff[1]
  if (is.null(DNA_2)) DNA_2=cutoff[3]
  abline(v=Find_one_pop(xi))
  abline(v=c(DNA_1,DNA_2),col='red')
  
  i=Event_length;j=DNA;
  xi=    (exprs(data_tmp)[exprs(data_tmp)[,i]>1 & exprs(data_tmp)[,j]>1 & exprs(data_tmp)[,i]<150 &
                            !is.na(exprs(data_tmp)[,j]) & !is.na(exprs(data_tmp)[,i]),i]); 
  xj=log2(exprs(data_tmp)[exprs(data_tmp)[,i]>1 & exprs(data_tmp)[,j]>1 & exprs(data_tmp)[,i]<150 &
                            !is.na(exprs(data_tmp)[,j]) & !is.na(exprs(data_tmp)[,i]),j]); #even length vs DNA
  smoothScatter((xi),(xj),xlab=parameters(data)$name[i],ylab=parameters(data)$desc[j],
                xlim=c(sort(xi)[0.01*length(xi)]-1,1.25*sort(xi)[0.99*length(xi)]),
                ylim=c(sort(xj)[0.01*length(xj)]-1,1+sort(xj)[0.99*length(xj)]),
                nrpoints = 10000)
  nk <- 11
  my.cols <- rev(brewer.pal(nk, "RdYlBu"))
  z <- kde2d(xi,xj,n=50)
  contour(z, drawlabels=FALSE, nlevels=nk, col=my.cols, add=TRUE, lwd=2)
  
  title(paste('data(',as.character(dim(exprs(data_tmp))[1]),')',sep=''))
  cutoff=Find_one_pop(xi[!is.na(xi)])
  if (is.null(Event_length_1)) Event_length_1=cutoff[1]
  i1=Event_length_1
  if (is.null(Event_length_2)) Event_length_2=cutoff[3]
  i2=Event_length_2
  cutoff=Find_one_pop(xj[!is.na(xj)])
  if (is.null(DNA_1)) DNA_1=(cutoff[1])
  j1=DNA_1
  if (is.null(DNA_2)) DNA_2=(cutoff[3])
  j2=DNA_2
  sqrcut <- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
  lines((rbind(sqrcut,sqrcut[1,])),col='red')
  sqrcut <- matrix(c(i1,i1,i2,i2,2^j1,2^j2,2^j2,2^j1),ncol=2,nrow=4)
  colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
  data_0<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= sqrcut)))
  text(i2,(j1+j2)/2,as.character(dim(exprs(data_0))[1]),col='red',cex=2)
  
  
  xi=log2(exprs(data_tmp)[exprs(data_tmp)[,DeadLive]>1 & !is.na(exprs(data_tmp)[,DeadLive]), DeadLive])
  plot(density(xi),main=parameters(data)$desc[DeadLive])
  abline(v=Find_two_pop(xi))
  abline(v=Find_two_pop(xi)[c(1)],col='red')
  
  
  data_tmp=data_0
  i=DNA;j=DeadLive; 
  xi=log2(exprs(data_tmp)[exprs(data_tmp)[,i]>1 & exprs(data_tmp)[,j]>1 & 
                            !is.na(exprs(data_tmp)[,j]) & !is.na(exprs(data_tmp)[,i]),i]); 
  xj=log2(exprs(data_tmp)[exprs(data_tmp)[,i]>1 & exprs(data_tmp)[,j]>1 & 
                            !is.na(exprs(data_tmp)[,j]) & !is.na(exprs(data_tmp)[,i]),j]); #DNA vs live/Dead
  smoothScatter((xi),(xj),xlab=parameters(data)$desc[i],ylab=parameters(data)$desc[j],
                xlim=c(sort(xi)[0.01*length(xi)]-1,1+sort(xi)[1*length(xi)]),
                ylim=c(sort(xj)[0.01*length(xj)]-1,1+sort(xj)[1*length(xj)]),
                nrpoints = 10000)
  my.cols <- rev(brewer.pal(nk, "RdYlBu"))
  z <- kde2d(xi,xj,n=50)
  contour(z, drawlabels=FALSE, nlevels=nk, col=my.cols, add=TRUE, lwd=2)
  title(paste('data_0(',as.character(length(xi)),')',sep=''))
  cutoff=Find_one_pop(xi)
  i1=sort(xi)[1]
  i2=sort(xi)[length(xi)]
  cutoff=Find_two_pop(xj)
  j1=sort(xj)[1]
  j2=cutoff[1]
  
  #	j2= 10
  
  sqrcut <- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
  lines((rbind(sqrcut,sqrcut[1,])),col='red')
  colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
  data_1<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= 2^sqrcut)))
  text(i2,(j1+j2)/2,as.character(dim(exprs(data_1))[1]),col='red',cex=2)
  
  
  marker_filter=NULL
  
  for(i in 1:marker_num)
  {
    print(marker_name[i])
    data_tmp=data_1
    xi=log2(exprs(data_tmp)[exprs(data_tmp)[,marker_id[i]]>1,marker_id[i]])
    plot(density(xi),main=parameters(data)$desc[marker_id[i]])
    if (length(Find_two_pop(xi)) == 2){
      print(paste0("Marker ",marker_name[i]," does not have a cutting point."))
      abline(v=Find_two_pop(xi))
    }else{
      marker_filter[i]=2^Find_two_pop(xi)[c(1)]
      abline(v=Find_two_pop(xi))
      abline(v=log2(marker_filter[i]),col='red',cex=2)
    }
    #if(marker_name[i]=='CD45B') marker_filter[i]=2^6
    #if(marker_name[i]=='CD45C') marker_filter[i]=2^7
  }
  #marker_filter
  
  
  data_tmp=data_1
  i=DNA;j=marker_id[1]; 
  xi=log2(exprs(data_tmp)[exprs(data_tmp)[,i]>0 & exprs(data_tmp)[,j]>0,i]); 
  xj=log2(exprs(data_tmp)[exprs(data_tmp)[,j]>0 & exprs(data_tmp)[,i]>0,j]); #DNA vs CD45A
  smoothScatter((xi),(xj),xlab=parameters(data)$desc[i],ylab=parameters(data)$desc[j],
                xlim=c(sort(xi)[1]-1,1+sort(xi)[1*length(xi)]),
                ylim=c(sort(xj)[1]-1,1+sort(xj)[1*length(xj)]),
                nrpoints = 10000)
  my.cols <- rev(brewer.pal(nk, "RdYlBu"))
  z <- kde2d(xi,xj,n=50)
  contour(z, drawlabels=FALSE, nlevels=nk, col=my.cols, add=TRUE, lwd=2)
  title(paste('data_1(',as.character(length(xi)),')',sep=''))
  #abline(v=j1,col='red')
  i1=0
  i2=2^sort(xi)[length(xi)]
  #j1=2^Find_two_pop(log2(exprs(data_tmp)[exprs(data_tmp)[,j]>1,j]))[1]
  j1=marker_filter[1]
  j2=2^sort(xj)[length(xj)]
  abline(h=log2(j1),col='red')
  sqrcut <- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
  colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
  data_A<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= sqrcut)))
  text(sort(xi)[1]+1,1+sort(xj)[1*length(xj)],as.character(dim(exprs(data_A))[1]),col='red',cex=2)
  i1=0
  i2=2^sort(xi)[length(xi)]
  j1=0
  #j2=2^Find_two_pop(log2(exprs(data_tmp)[exprs(data_tmp)[,j]>1,j]))[1]
  j2=marker_filter[1]
  sqrcut<- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
  colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
  data_a<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= sqrcut)))
  text(sort(xi)[1]+1,sort(xj)[1]-1,as.character(dim(exprs(data_a))[1]),col='green',cex=2)
  dev.off()
  
  pdf(paste(plot_name,'_plot_all_data.pdf',sep=''),pointsize=30,width=60, height=40)
  par(mfrow=c(4,6))
  
  for(k in 1:nrow(barcode_table))
    #k=1
  {
    title_name=paste(barcode_table[k,dim(barcode_table)[2]],'_',barcode_table[k,2],sep='')
    #print(title_name)
    if (barcode_table[k,2]=='A') data_tmp=data_A
    if (barcode_table[k,2]=='a') data_tmp=data_a
    
    i = marker_id[1]
    x_lab = parameters(data)$desc[i]
    
    for (n in 2:marker_num)
    {
      title_name=paste(title_name,barcode_table[k,1+n],sep='')
      #print(title_name)
      i=marker_id[1]; 
      j=marker_id[n];
      
      if (n != 2) {
        
        x_lab = paste(x_lab, parameters(data)$desc[marker_id[n-1]], sep ="_and_")
      }
      
      xi=log2(exprs(data_tmp)[exprs(data_tmp)[,i]>0 & exprs(data_tmp)[,j]>0,i]); 
      xj=log2(exprs(data_tmp)[exprs(data_tmp)[,j]>0 & exprs(data_tmp)[,i]>0,j]); #DNA vs CD45A
      smoothScatter((xi),(xj),xlab=x_lab,ylab=parameters(data)$desc[j],
                    xlim=c(sort(xi)[1]-1,1+sort(xi)[1*length(xi)]),
                    ylim=c(sort(xj)[1]-1,1+sort(xj)[1*length(xj)]),
                    nrpoints = 10000)
      my.cols <- rev(brewer.pal(nk, "RdYlBu"))
      z <- kde2d(xi,xj,n=50)
      contour(z, drawlabels=FALSE, nlevels=nk, col=my.cols, add=TRUE, lwd=2)
      title(paste(title_name,'(',as.character(length(xi)),')',sep=''))
      i1=0
      i2=2^sort(xi)[length(xi)]
      if (barcode_table[k,1+n]==c('A','B','C','D','E')[n])
      {
        #j1=2^Find_two_pop(log2(exprs(data_tmp)[exprs(data_tmp)[,j]>1,j]))[1]
        j1=marker_filter[n]
        j2=2^sort(xj)[length(xj)]
        abline(h=log2(j1),col='red')
        sqrcut <- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
        colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
        data_1<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= sqrcut)))
        text(sort(xi)[1]+1,1+sort(xj)[1*length(xj)],as.character(dim(exprs(data_1))[1]),col='red',cex=2)
      }
      if (barcode_table[k,1+n]==c('a','b','c','d','e')[n])
      {
        j1=0
        #j2=2^Find_two_pop(log2(exprs(data_tmp)[exprs(data_tmp)[,j]>1,j]))[1]
        j2=marker_filter[n]
        abline(h=log2(j2),col='red')
        sqrcut <- matrix(c(i1,i1,i2,i2,j1,j2,j2,j1),ncol=2,nrow=4)
        colnames(sqrcut) <- c(parameters(data)$name[i],parameters(data)$name[j])
        data_1<-Subset(data_tmp, filter(data_tmp, polygonGate(filterId="test", .gate= sqrcut)))
        text(sort(xi)[1]+1,sort(xj)[1]-1,as.character(dim(exprs(data_1))[1]),col='green',cex=2)
      }
      data_tmp=data_1
    }
    suppressWarnings(write.FCS(data_tmp, file_list[k]))
    print(c(file_list[k],as.character(dim(exprs(data_1))[1])))
  }
  dev.off()
  
  
}


Find_one_pop<-function(x)
{
  x_den=density(x)
  x_peak_pos=findPeaks(x_den$y) - 1
  x_val_pos=c(1,findPeaks(-x_den$y)-1,length(x_den$x))
  max_pos=x_peak_pos[1]
  for (i in 1:length(x_peak_pos))
  {
    if (x_den$y[x_peak_pos[i]]>x_den$y[max_pos]) max_pos=x_peak_pos[i]
  }
  peak_val=x_den$x[max_pos]
  if (length(x_val_pos)>0)
  {
    left_pos=x_val_pos[1]
    right_pos=x_val_pos[1]
    for (i in 1:length(x_val_pos))
    {
      if (x_den$x[x_val_pos[i]] < peak_val)  
      {
        if (peak_val - x_den$x[x_val_pos[i]] < abs(peak_val-x_den$x[left_pos])) left_pos=x_val_pos[i]
      }
      if (x_den$x[x_val_pos[i]] > peak_val) 
      {
        if (x_den$x[right_pos] - peak_val<0) right_pos=x_val_pos[i]
        if (x_den$x[x_val_pos[i]] - peak_val < abs(x_den$x[right_pos] - peak_val)) right_pos=x_val_pos[i]
      }
    }
    left_val=x_den$x[left_pos]
    right_val=x_den$x[right_pos]
  }
  #    else
  #    {
  #        left_val=x_den$x[1]
  #        right_val=x_den$x[length(x_den$x)]
  #    }
  return(c(left_val,peak_val,right_val))
}



Find_two_pop<-function(x)
{
  x_den=density(x)
  x_peak_pos=findPeaks(x_den$y)-1
  x_val_pos=findPeaks(-x_den$y)-1
  peak1=sort(x_den$y[findPeaks(x_den$y)-1])[length(x_peak_pos)]
  if (length(sort(x_den$y[findPeaks(x_den$y)-1])) >1){
    
    peak2=sort(x_den$y[findPeaks(x_den$y)-1])[length(x_peak_pos)-1]
    
  } else {
    
    peak2=sort(x_den$y[findPeaks(x_den$y)-1])[length(x_peak_pos)]
  }
  for (i in 1:length(x_peak_pos))
  {
    if (x_den$y[x_peak_pos[i]]==peak1) max_pos1=x_peak_pos[i]
    if (x_den$y[x_peak_pos[i]]==peak2) max_pos2=x_peak_pos[i]
    
  }
  peak_val1=x_den$x[max_pos1]
  peak_val2=x_den$x[max_pos2]
  if (length(x_val_pos)>0)
  {
    sep_pos=-1
    
    for (i in 1:length(x_val_pos))
    {
      
      if ((x_den$x[x_val_pos[i]] - peak_val1)*(x_den$x[x_val_pos[i]] - peak_val2)<0)  
      {
        if (sep_pos==-1) sep_pos=x_val_pos[i]
        if (x_den$y[x_val_pos[i]] < x_den$y[sep_pos]) sep_pos=x_val_pos[i]
      }
    }
    sep_val=x_den$x[sep_pos]	
  }
  
  if (exists("sep_val")){
    return(c(sep_val,peak_val1,peak_val2))
  } else {
    return(c(peak_val1,peak_val2))
  }
}