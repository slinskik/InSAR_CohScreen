###########################
# Check coherance with time
#
# created 2017-10-10 k Slinski
###########################
library(ncdf4)
library(raster)
library(rasterVis)
library(viridis)
library(scales)
library(lubridate)
library(plyr)
library(ggplot2)
library(reshape)

rm(list=ls())
setwd('K:\\phyrev\\_WA_INSAR\\GMTSAR\\Ara\\20171012')

# #**** order matters!! load dependencies before main package 
# library(sp, lib.loc="/home/slinskik/rPkg")
# library(viridisLite, lib.loc="/home/slinskik/rPkg")
# library(methods, lib.loc="/home/slinskik/rPkg")
# library(plyr, lib.loc="/home/slinskik/rPkg")
# library(raster, lib.loc="/home/slinskik/rPkg")
# library(rasterVis, lib.loc="/home/slinskik/rPkg")
# library(viridis, lib.loc="/home/slinskik/rPkg")
# library(scales, lib.loc="/home/slinskik/rPkg")
# library(ncdf4, lib.loc="/home/slinskik/rPkg")
# library(lubridate, lib.loc="/home/slinskik/rPkg")
# library(ggplot2, lib.loc="/home/slinskik/rPkg")
# library(reshape, lib.loc="/home/slinskik/rPkg")
# library(lattice, lib.loc="/home/slinskik/rPkg")
# library(latticeExtra, lib.loc="/home/slinskik/rPkg")
# library(RColorBrewer, lib.loc="/home/slinskik/rPkg")

# ------------------------------------------------
#   Average Monthly Coherence Function
# ------------------------------------------------
coh_mo=function(pairs, s, ext){
  list=c()
  list3=c()
  months=c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
  mean=c()
  n=c()
  for (m in 1:12){
    pairsU=pairs[pairs$month==m|pairs$month2==m,]
    pairsU=pairsU$dir
    if (length(pairsU)==1){
      file=paste0('intf_all/',pairsU,'/corr_ll.grd')
      ras=extend(raster(file),ext)
    } else{
      file=paste0('intf_all/',pairsU[1],'/corr_ll.grd')
      ras=extend(raster(file),ext)
      for(i in 2:length(pairsU)){
        file=paste0('intf_all/',pairsU[i],'/corr_ll.grd')
        ras=stack(ras, extend(raster(file),ext))
      }
    }
    ras=mean(ras, na.rm=T)
    n=rbind(n, length(pairsU))
    name=paste0('month_', as.character(m))
    names(ras)=months[m]
    assign(name, ras)
    list=rbind(list, name)
  }
  list=as.vector(list)
  stack=get(list[1])
    for (i in 2:length(list)){
    ras=get(list[i])
    stack=stack(stack,ras)
  }
  mean=cellStats(stack, 'mean', na.rm=T)
  for (i in 1:5){
    stackU=reclassify(stack, c(-Inf,.1*i,NA))
    names(stackU)=months
    stackT=reclassify(stackU, c(-Inf,Inf,1))
    sum=as.vector(cellStats(stackT, 'sum'))
    stackT=reclassify(stack, c(-Inf,Inf,1))
    sum2=as.vector(cellStats(stackT, 'sum'))
    frac=sum/sum2
    assign(paste0('stack', as.character(i)), stackU)
    assign(paste0('frac', as.character(i)), frac)
  }
  out=list(stack, data.frame(months, mean, n, frac1, frac2, frac3, frac4, frac5), stack1, stack2, stack3, stack4, stack5)
  names(out)=c('stack', 'df', 'stack1', 'stack2', 'stack3', 'stack4', 'stack5')
  return(out)
}

# ------------------------------------------------
#   Get pairs
# ------------------------------------------------
pairs=list.dirs('intf_all', full.names = F, recursive = F)
pairs=pairs[-1*grep('_$', pairs)]
dates=matrix(unlist(strsplit(pairs, '_')), ncol=2, byrow=T)
pairs=data.frame(dir=pairs, d1=dates[,1], d2=dates[,2], stringsAsFactors = F)
pairs$d1=as.Date(pairs$d1, '%Y%j')+1
pairs$d2=as.Date(pairs$d2, '%Y%j')+1
pairs$delta=as.numeric(pairs$d2-pairs$d1)
pairs$month=month(pairs$d1)
pairs$month2=month(pairs$d2)
ddply(pairs, 'delta', summarize, count=length(month))
#pairs=subset(pairs, delta<=48)
sets=sort(unique(pairs$delta))
rm(dates)
# ------------------------------------------------
#   Get extent
# ------------------------------------------------
exts=c()
for(i in pairs$dir){
  #print(i)
  file=paste0('intf_all/',i,'/corr_ll.grd')
  ext=extent(raster(file))
  x=ext[1:2]
  y=ext[3:4]
  exts=rbind(exts,as.vector(ext))
}
ext=extent(min(exts[,1]), max(exts[,2]), min(exts[,3]), max(exts[,4]))
rm(exts, i, file, x, y)

# ------------------------------------------------
#   Average Coherance Across All Pairs
# ------------------------------------------------
file=paste0('intf_all/',pairs$dir[1],'/corr_ll.grd')
ras=extend(raster(file),ext)
for(i in 2:length(pairs$dir)){
  #print(i)
  file=paste0('intf_all/',pairs$dir[i],'/corr_ll.grd')
  ras=stack(ras, extend(raster(file),ext))
}
Coh_mean=mean(ras, na.rm=T)
mean=cellStats(Coh_mean, 'mean', na.rm=T)

tiff('coh_area.tiff', width=10, height=10, units='in', res=300)
levelplot(Coh_mean, margin=F, xlab='Longitude', ylab='Latitude',par.settings = viridisTheme, at=unique(c(seq(0, 1, length=20))),
          main = paste0('Average Coherence \n','(mean=',as.character(round(mean, 2)), ', n=',as.character(length(pairs$dir)), ')'))
dev.off()
rm(Coh_mean, i, file, ras)

# ------------------------------------------------
#   Coherence by timestep
# ------------------------------------------------
means=c()
months=c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
s=365
out=coh_mo(pairs, s, ext)
file=paste0('coh_',as.character(s), 'days.tiff')
tiff()
  gplot(out$stack, maxpixels=50000)+geom_tile(aes(fill = value))+facet_wrap(~ variable, nrow=2)+
    scale_fill_viridis(option="viridis", na.value="transparent", limits=c(0,1), name='Coherence')+
    xlab('Longitude')+ylab('Latitude')+ggtitle(paste0('Average Coherence of ',as.character(s), ' Day Pairs'))+
    theme_bw()+coord_equal(expand=F)+
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=10),  
          axis.title.x=element_text(size=12), axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          plot.title=element_text(size=16, hjust=.5),
          legend.text=element_text(size=10), legend.title=element_blank(), legend.key = element_blank(),
          panel.border = element_rect(linetype = 1, size=1,colour = "black"),
          strip.background = element_rect(colour=NA, fill=NA),strip.text.x=element_text(size=12))
  ggsave(filename=file, width=16, height=9, dpi = 300)
dev.off()
  
stack=stack(out$stack, out$stack1, out$stack2, out$stack3, out$stack4, out$stack5)
file=paste0('frac_coh_',as.character(s), 'days.tiff')
tiff()
  gplot(stack, maxpixels=30000)+geom_tile(aes(fill = value))+facet_wrap(~ variable, ncol=12)+
    scale_fill_viridis(option="viridis", na.value="transparent", limits=c(0,1), name='Coherence')+
    xlab('Longitude')+ylab('Latitude')+ggtitle(paste0('Average Coherence of ',as.character(s), ' Day Pairs'))+
    theme_bw()+coord_equal(expand=F)+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),  
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          line=element_blank(),
          plot.title=element_text(size=16, hjust=.5),
          legend.text=element_text(size=8), legend.title=element_blank(), legend.key = element_blank(),
          panel.border = element_rect(linetype = 1,colour = "black"),
          strip.background = element_rect(colour=NA, fill=NA),strip.text.x=element_blank())
  ggsave(filename=file, width=16, height=9, dpi = 300)
dev.off()
  
