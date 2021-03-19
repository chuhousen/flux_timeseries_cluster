rm(list=ls()) # clean working memory

require("pdc")
require("proxy")
require("dtw")
require("graphics")
require("rgdal")
require("rworldmap")
require("ggplot2")   # for general plotting
require("ggmap")    # for fortifying shapefiles
require("raster")
require("RColorBrewer")
require("ggrepel")
require("ape") ##  for displaying more appealing trees
require("dtwclust")
#require("TSclust")  lost dependency

source("math_utility.R")
#source("BIF_working.R")
source("get_ecomap.R")

ver<-"20201117"

path1<-"D:\\Housen\\Flux\\Data-exploring\\01_AMF_Processing\\02_local_working_hchu\\"
path2<-"D:\\AmeriFlux-Data\\"
path3<-"D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\"
base.in<-paste(path2,"00_olaf_data_ameriflux\\BASE-BADM\\",sep="")
map.in<-paste(path3,"Ecoregion\\",sep="")
path.in<-paste(path3,"diurnal-seasonality\\",ver,"\\",sep="")
path.out<-paste(path3,"cluster-output\\",ver,"\\",sep="")

## extract badm
full.ls<-read.csv(paste(path.in,"ALL_BASE_site_short_list.csv",sep=""),
                  na="-9999",
                  header=T,
                  stringsAsFactors=F)
full.ls<-full.ls[!is.na(full.ls$SITE_ID),]
rownames(full.ls)<-paste(full.ls$SITE_ID)

#### downselecting 
full.ls<-full.ls[((full.ls$LE_F==1&full.ls$H_F==1)|
                    (full.ls$NEE_F==1&full.ls$H_F==1)|
                    (full.ls$LE_F==1&full.ls$NEE_F==1))&
                   full.ls$USTAR_F==1&
                   full.ls$TA_F==1&
                   full.ls$VPD_F==1,]

## original colnames from diurnal-seasonal files, revising the colname
sel.var<-c("FC","SC","NEE","LE","H","USTAR","NETRAD",
           "TA","VPD","WS",
           "G","SWC",#"TS",
           "SW_IN","LW_IN","LW_OUT","SW_OUT")
drop.ls<-c("FC","SC","WS","G","LW_IN","LW_OUT","SW_OUT")

for(l in 1:length(sel.var)){
  if(sel.var[l] %in% drop.ls){
    full.ls<-full.ls[,-which(colnames(full.ls)==sel.var[l])]
    full.ls<-full.ls[,-which(colnames(full.ls)==paste0(sel.var[l],"_F"))]
  }else{
    colnames(full.ls)[which(colnames(full.ls)==sel.var[l])]<-paste0(sel.var[l],"_availability")
    colnames(full.ls)[which(colnames(full.ls)==paste0(sel.var[l],"_F"))]<-paste0(sel.var[l],"_availability_after_filling")  
  }
}
#colnames(full.ls)

## find a list of preprocessed diurnal-seasonal files
src.list.in<-list.files(path.in)
src.list.in<-src.list.in[which(substr(src.list.in,
                                      start=nchar(src.list.in)-3,
                                      stop=nchar(src.list.in))==".csv")]
src.list.in<-src.list.in[which(substr(src.list.in,
                                      start=nchar(src.list.in)-9,
                                      stop=nchar(src.list.in)-4)=="MEDIAN")]

src.list.in<-src.list.in[which(substr(src.list.in,
                                      start=1,
                                      stop=6) %in% full.ls$SITE_ID)]

### Drop problematic cases 
## drop.ls<-c("BR-CST","CA-SCB","US-Elm","US-Esm","US-Hn2","US-NC2",) 

## 
n.grp<-c(10:30)  # target range of tree branches to keep

target.var.ls<-c("NEE","LE","H",
                 "USTAR","NETRAD",
                 "TA","VPD",#"WS",
                 "SWC","SW_IN")

### read in Ecoregion map
ecomap1<-readOGR(paste(map.in,"wwf_terr_ecos_oRn",sep=""),
                 layer=c("wwf_terr_ecos_oRn"),verbose=F)
ecomap3<-readOGR(paste(map.in,"NA_CEC_Eco_Level3",sep=""),
                 layer=c("NA_CEC_Eco_Level3"),verbose=F)

#proj4string(ecomap1)
ecomap1 <- spTransform(ecomap1, 
                       CRS("+proj=longlat +datum=WGS84"))
#proj4string(ecomap3)
ecomap3 <- spTransform(ecomap3, 
                       CRS("+proj=longlat +datum=WGS84"))

### work on extraction of ecoregion for each site
pts<-SpatialPoints(full.ls[,c("GRP_LOCATION.LOCATION_LONG",
                              "GRP_LOCATION.LOCATION_LAT")],
                   proj4string = CRS("+proj=longlat +datum=WGS84"))
full.ls<-cbind.data.frame(full.ls,
                          over(pts,ecomap1)[,c(4,5,6,7)],
                          over(pts,ecomap3)[,c(1:6)])

## Work by variables
for(l1 in 1:length(target.var.ls)){
  
  target.var<-target.var.ls[l1]
  
  src.ls<-sub("_MEDIAN.csv","",src.list.in)
  src.ls<-src.ls[!duplicated(src.ls)]
  
  site.ls<-substr(src.ls,start=1,stop=6)  
  res.ls<-substr(src.ls,start=8,stop=9)
  ver.ls<-substr(src.ls,start=11,stop=nchar(src.ls))
  
  print(paste("################################################"))
  print(paste("######  ",target.var,"  ##################"))
  
  #########################################################################################################
  #### work on prepaare site list
  
  ## read in pre-processed diurnal-seasonal files
  drop.case<-NULL
  for(i0 in 1:length(src.ls)){
    
    data.tmp<-read.csv(paste(path.in,src.ls[i0],"_MEDIAN.csv",sep=""),header=T,
                       na.strings = "-9999",stringsAsFactors=F)
    data.tmp<-data.tmp[,target.var]
    
    ## check whether missing values
    if(sum(is.na(data.tmp))!=0){
      drop.case<-c(drop.case,site.ls[i0])
    }
  }
  
  print(paste("[Info]",length(src.ls)-length(drop.case),"has sufficient data"))
  
  if(length(drop.case)>0){
    for(i1 in 1:length(drop.case)){
      src.ls<-src.ls[site.ls!=drop.case[i1]]
      ver.ls<-ver.ls[site.ls!=drop.case[i1]]
      res.ls<-res.ls[site.ls!=drop.case[i1]]
      site.ls<-site.ls[site.ls!=drop.case[i1]]
    }
    print(paste("[Info]","drop",length(drop.case),"sites because of gaps"))
    print(paste("[Info]","drop",paste(drop.case,collapse=" ")))
  }
  lat.ls<-as.numeric(full.ls[site.ls,"GRP_LOCATION.LOCATION_LAT"])
  lon.ls<-as.numeric(full.ls[site.ls,"GRP_LOCATION.LOCATION_LONG"])
  
  ################################################################################################################
  #### prepare full diurnal-seasonal time series 
  
  num.ts <- length(src.ls) # number of time series
  num.dim <- length(target.var) # number of dimensions
  len.ts <- 24*24 # length of time series
  
  data.pre <- NULL #array(dim = c(len.ts, num.ts, num.dim))
  
  for(i2 in 1:length(src.ls)){
    data.tmp<-read.csv(paste(path.in,src.ls[i2],"_MEDIAN.csv",sep=""),
                       header=T,
                       na.strings = "-9999",
                       stringsAsFactors=F)
    data.tmp<-data.tmp[,target.var]
    data.pre<-cbind(data.pre,data.tmp)
    colnames(data.pre)[which(colnames(data.pre)=="data.tmp")]<-paste(site.ls[i2])
  }
  
  write.csv(data.pre,
            paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,"-compiled_ts.csv",sep=""),
            row.names = F)
  
  #############################################################################################################
  #### Clustering Time series
  
  # In the following call to tsclust, specifying the value of k indicates the number of desired
  # clusters, so that the cutree function is called internally. Additionally, the shape extraction
  # function is provided in the centroid argument so that, once the k clusters are obtained, their
  # prototypes are extracted. Therefore, the series are z-normalized by means of the zscore
  # function. The seed is provided because of the randomness in shape extraction when choosing
  # a reference series
  
  #data("uciCT")
  distMatrix <- proxy::dist(t(data.pre),
                            method="dtw_basic")
  
  hc1 <- tsclust(t(data.pre), 
                type = "h", 
                k = n.grp,
                seed = 899,
                distance = "dtw_basic", 
                centroid = "pam",
                control = hierarchical_control(method = "average"))
  
  names(hc1) <- paste0("k_", n.grp)
  cvi.est<-as.data.frame(t(sapply(hc1, cvi, type = "internal"))) 
  cvi.est<-cvi.est[,-which(colnames(cvi.est)=="SF")]
  
  for(i in 1:ncol(cvi.est)){
    if(colnames(cvi.est)[i] %in% c("DB","DBstar","COP")){ ## minimize
      cvi.est[,i] <- 1-cvi.est[,i]/max(cvi.est[,i])
    }else{ ## maximize
      cvi.est[,i] <- cvi.est[,i]/max(cvi.est[,i])
    }
  }
  
  png(paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,"-ncluster_cvi.png",sep=""),
      width=6,
      height=4.5,
      units="in",
      pointsize=10,
      res=300)
  plot(n.grp,cvi.est[,1]/max(cvi.est[,1]),
       type="l",xlim=range(n.grp),ylim=c(0,1.2),las=1,
       xlab="number of cluster",ylab="normalized cluster validity index",lwd=2,
       col=rainbow(6)[1],
       main=target.var)
  for(i in 2:ncol(cvi.est)){
    lines(n.grp,cvi.est[,i]/max(cvi.est[,i]),
          col=rainbow(6)[i],lwd=2)      
  }
  legend(10,1.2,
         lty=1,
         col=c("black",rainbow(6)),
         ncol=4,
         legend = c("MEAN",colnames(cvi.est)),
         lwd=1.5,bty="n")
  lines(n.grp,apply(cvi.est,1,mean),lwd=3)
  abline(v=n.grp[which(apply(cvi.est,1,mean)==max(apply(cvi.est,1,mean)))],lty=4)
  dev.off()
  
  ## Use cluster number from CH (Calinski-Harabasz) index
  n.grp.opt<-n.grp[which(apply(cvi.est,1,mean)==max(apply(cvi.est,1,mean)))]
  
  hc <- tsclust(t(data.pre), 
                 type = "h", 
                 k = n.grp.opt,
                 seed = 899,
                 distance = "dtw_basic", 
                 centroid = "pam",
                 control = hierarchical_control(method = "average"))
  
  grp.ls<-cutree(hc,k=n.grp.opt)
  
  write.csv(as.matrix(distMatrix),
            paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,"-distMatrix.csv",sep=""),
            quote = T)
  
  # return group ID to full list
  full.ls<-data.frame(full.ls,
                      tmp=NA)
  
  for(i3 in 1:length(grp.ls)){
    full.ls$tmp[which(full.ls$SITE_ID==paste(names(grp.ls)[i3]))]<-grp.ls[i3]
  }
  colnames(full.ls)[which(colnames(full.ls)=="tmp")]<-paste(target.var,"_clust_group",sep="")
  
  ## plot trees
  png(paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,"-tree.png",sep=""),
      width=9,
      height=9,
      units="in",
      pointsize=10,
      res=300)
  par(mfrow=c(1,1),mar=c(4.5,4.5,4.5,4.5))
  plot(as.phylo(hc),
       type="fan",
       tip.color=rainbow(n.grp.opt)[grp.ls],
       label.offset=0.5, 
       #show.node.label=T,
       cex=0.9)
  #text(0,0,)
  dev.off()
  
  ## plot grouped diurnal-seasonal plots
  png(paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,".png",sep=""),
      width=11,
      height=9,
      units="in",
      pointsize=9,
      res=300)
  
  par(mfrow=c(5,ceiling(n.grp.opt/5)),
      mar=c(0,0,0,0),
      oma=c(3,5,0.5,0.5))
  
  for(j2 in 1:n.grp.opt){
    
    plot(0,0,
         xlim=c(-20,len.ts+20),
         ylim=range(data.pre),
         type="n",
         xlab="",ylab="",
         xaxt="n",yaxt="n")
    
    data.pre.sub<-as.data.frame(data.pre[,which(grp.ls==j2)])
    
    data.pre.sub.mean<-apply(data.pre.sub,1,na.mean)
    
    for(j1 in 1:ncol(data.pre.sub)){
      lines(data.pre.sub[,j1],
            col="grey")
    }
    legend(-20,range(data.pre)[2],
           paste(names(grp.ls[grp.ls==j2])),
           col="darkgrey",
           cex=0.7,
           bty="n",
           ncol=ceiling(ncol(data.pre.sub)/22))
    
    if((j2-1)==floor((j2-1)/ceiling(n.grp.opt/5))*ceiling(n.grp.opt/5)){
      axis(2,seq(round(range(data.pre)[1],digits=0),
                 round(range(data.pre)[2],digits=0),
                 length.out=6),
           las=2)
    }
    #if(ncol(data.pre.sub)>1) 
    lines(data.pre.sub.mean,
          col=rainbow(n.grp.opt)[j2],
          lwd=1.5)
  }
  mtext(side=1,"Date-Hour",outer = T,line=2,cex=1.5)
  mtext(side=2,paste(target.var),outer=T,line=3.2,cex=1.5)
  dev.off()
  
  ### subgroup plot
  sub.grp.ls<-as.numeric(names(table(grp.ls))[table(grp.ls)>5])
  
  if(length(sub.grp.ls)>0){
    for(i in 1:length(sub.grp.ls)){
      
      data.pre.sub<-as.data.frame(data.pre[,which(grp.ls==sub.grp.ls[i])])
      data.pre.sub.name<-colnames(data.pre.sub)
      data.pre.sub.mean<-apply(data.pre.sub,1,na.mean)
      
      if(ncol(data.pre.sub)<=25){
        data.pre.sub.ls<-c(1:ncol(data.pre.sub))
      }else{
        data.pre.sub.ls<-ceiling(c(1:ncol(data.pre.sub))/(ceiling(ncol(data.pre.sub)/25)))
        data.pre.sub.ls[data.pre.sub.ls<1]<-1
        data.pre.sub.ls[data.pre.sub.ls>25]<-25
      }
      
      png(paste(path.out,"AMF-diurnal-seasonal-cluster-",target.var,"_subgroup_",sub.grp.ls[i],".png",sep=""),
          width=11,height=9,units="in",pointsize=9,res=300)
      par(mfrow=c(5,5),mar=c(0,0,0,0),oma=c(3,5,0.5,0.5))
      
      for(j in 1:min(c(max(data.pre.sub.ls),25))){
        
        plot(0,0,
             xlim=c(0,len.ts),
             ylim=range(data.pre.sub),
             type="n",ylab="",xaxt="n",yaxt="n",xlab="")
        
        lines(data.pre.sub.mean,col="black",lwd=1.5)
        
        data.pre.sub.sub<-as.data.frame(data.pre.sub[,which(data.pre.sub.ls==j)])
        colnames(data.pre.sub.sub)<-data.pre.sub.name[which(data.pre.sub.ls==j)]
        
        for(j1 in 1:ncol(data.pre.sub.sub)){
          lines(data.pre.sub.sub[,j1],
                col=rainbow(ncol(data.pre.sub.sub))[j1])
        }
        legend(0,range(data.pre.sub)[2],
               paste(colnames(data.pre.sub.sub)),
               lty=c(1),
               col=rainbow(ncol(data.pre.sub.sub)),
               cex=1,bty="n")
        if(j==1|j==6|j==11|j==16|j==21){
          axis(2,seq(round(range(data.pre.sub)[1],digits=0),
                     round(range(data.pre.sub)[2],digits=0),
                     length.out=6),las=2)
        }
      }
      
      mtext(side=1,"Date-Hour",outer = T,line=1,cex=1.5)
      mtext(side=2,paste(target.var),outer=T,line=3.2,cex=1.5)
      dev.off()
    }
  }
  
  print(paste("################################################"))
  print(paste("                                                "))
  
}

write.csv(full.ls,paste(path.out,"ALL_BASE_site_short_list2.csv",sep=""),quote=F,row.names=F)
