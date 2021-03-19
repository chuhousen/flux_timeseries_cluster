
#### v3  search for the latest BASE, based on olaf BASE-BADM folders
        # don't generate anything from newer BASE X-2, which has slightly diff header (no need to skipfirst 2 lines)
        #                                              & diff variables     
        # generate both HH and HR version if both BASE exist
        
#### v2  convert VPD from kPa to hPa

rm(list=ls()) # clean working memory

require("zoo")

path<-"D:\\Housen\\Flux\\Data-exploring\\01_AMF_Processing\\02_local_working_hchu\\diurnal-seasionality\\"
#path.in<-"D:\\Housen\\Flux\\Data-exploring\\01_AMF_Processing\\AMF_BASE_data\\"
path.in<-"D:\\Housen\\Flux\\Data-exploring\\01_AMF_Processing\\00_olaf_data_ameriflux\\BASE-BADM\\"
path.out<-"D:\\Housen\\Flux\\Data-exploring\\01_AMF_Processing\\02_local_working_hchu\\diurnal-seasionality\\20180503-version\\"

l.wd<-30
n.wd<-floor(365/l.wd)

na.median<-function(x){ifelse(sum(!is.na(x))>=0.25*length(x),median(x,na.rm=T),NA)}
upp.bd<-function(x){ifelse(sum(!is.na(x))>=0.25*length(x),quantile(x,probs=0.75,na.rm=T),NA)}
low.bd<-function(x){ifelse(sum(!is.na(x))>=0.25*length(x),quantile(x,probs=0.25,na.rm=T),NA)}
upp.bd2<-function(x){ifelse(sum(!is.na(x))>=0.25*length(x),quantile(x,probs=0.975,na.rm=T),NA)}
low.bd2<-function(x){ifelse(sum(!is.na(x))>=0.25*length(x),quantile(x,probs=0.025,na.rm=T),NA)}

limit.ls<-read.csv(paste(path,"Physical_Range_20170710.csv",sep=""),header=T,na=NA)

sel.var<-c("USTAR","TA",
           "WS","NEE_PI","FC","SC","H","SH","LE","SLE",
           "G","TS","RH","PA",
           "CO2","VPD_PI","SWC","NETRAD","PPFD_IN","SW_IN","SW_DIF",    
           "PPFD_OUT","SW_OUT","LW_IN","LW_OUT","H2O",   
           "APAR","PPFD_DIF","FAPAR","ZL")

src.list.in<-list.files(path.in)
src.list.in<-src.list.in[which(substr(src.list.in,start=nchar(src.list.in)-3,stop=nchar(src.list.in))==".zip")]

src.list<-substr(src.list.in,start=1,stop=nchar(src.list.in)-4)

### keep only X-1 version 
src.list<-src.list[substr(src.list,start=nchar(src.list),stop=nchar(src.list))=="1"]

pb<-txtProgressBar(min=0,max=length(src.list),style=3) ## progress bar
for(k in 151:length(src.list)){
  
  case.name<-src.list[k]
  file.grab<-unzip(paste(path.in,case.name,".zip",sep=""),list=TRUE)[,"Name"]
  
  case.ls<-file.grab[which(substr(file.grab,start=12,stop=15)=="BASE")]
  case.ls<-substr(case.ls,start=1,stop=nchar(case.ls)-4)
  
  ### handle multiple BASE at each site (e.g., HH & HR)
  for(m in 1:length(case.ls)){
    case<-case.ls[m]
   
    d.hr<-ifelse(substr(case,start=17,stop=18)=="HH",48,24)
    hr<-ifelse(d.hr==48,30,60)
    
    #if(substr(case,start=nchar(case),stop=nchar(case))=="1"){
    data.in<-read.table(unz(paste(path.in,case.name,".zip",sep=""),paste(case,".csv",sep="")),na=c("-9999"),header=T,sep=",",skip=2)  
    #}else{  ## for newer BASE
    #  data.in<-read.table(unz(paste(path.in,case.name,".zip",sep=""),paste(case,".csv",sep="")),na=c("-9999"),header=T,sep=",",skip=0)  
    #}
    
    var.name<-names(data.in)
    
    #data.in$TIMESTAMP_START<-strptime(data.in$TIMESTAMP_START,format="%Y%m%d%H%M")  
    #data.in$TIMESTAMP_END<-strptime(data.in$TIMESTAMP_END,format="%Y%m%d%H%M")  
    #data.in$TIMESTAMP<-strptime(data.in$TIMESTAMP_START+0.5*hr*60,format="%Y-%m-%d %H:%M:%S")
    
    if(length(which(colnames(data.in)=="TIMESTAMP_START"))==1){
      
      data.in$TIMESTAMP<-strptime(data.in$TIMESTAMP_START,format="%Y%m%d%H%M",tz="UTC")  
      data.in$TIMESTAMP<-strptime(TIMESTAMP+0.5*hr*60,format="%Y-%m-%d %H:%M:%S",tz="UTC")
      
    }else if(length(which(colnames(data.in)=="TIMESTAMP_END"))==1){
      
      data.in$TIMESTAMP<-strptime(data.in$TIMESTAMP_END,format="%Y%m%d%H%M",tz="UTC")  
      data.in$TIMESTAMP<-strptime(TIMESTAMP-0.5*hr*60,format="%Y-%m-%d %H:%M:%S",tz="UTC")
      print("warning: can't find TIMESTAMP_START columns")
      
    }else{
      
      print("error: can't find both TIMESTAMP columns")
      
    }   
    
    
    data.in$DOY2<-floor(data.in$TIMESTAMP$yday/l.wd)+1
    data.in$DOY2[which(data.in$DOY2==max(data.in$DOY2,na.rm=T))]<-n.wd
    data.in$HR2<-data.in$TIMESTAMP$hour+data.in$TIMESTAMP$min/60
    
    DOY2<-rep(c(1:(n.wd)),each=d.hr)
    HR2<-rep(c(1:d.hr)/(60/hr)-(hr/60)/2,times=n.wd)
    
    ### convert VPD from kPa to hPa
    data.in[,which(substr(var.name,start=1,stop=3)=="VPD")]<-data.in[,which(substr(var.name,start=1,stop=3)=="VPD")]*10
    
    ## pre-filter by physical limits
    for(l in 1:(ncol(data.in))){
      if(length(which(limit.ls$Name==names(data.in)[l]))>0){
        var.upp<-limit.ls$Max[which(limit.ls$Name==names(data.in)[l])]
        var.low<-limit.ls$Min[which(limit.ls$Name==names(data.in)[l])]
        if(!is.na(var.upp)){
          data.in[which(data.in[,l]>var.upp),l]<-NA
        }
        if(!is.na(var.low)){
          data.in[which(data.in[,l]<var.low),l]<-NA
        }
      }
    }
    
    temp.out1<-array(NA,dim=c(d.hr,n.wd))
    temp.out2<-array(NA,dim=c(d.hr,n.wd))
    temp.out3<-array(NA,dim=c(d.hr,n.wd))
    temp.out4<-array(NA,dim=c(d.hr,n.wd))
    temp.out5<-array(NA,dim=c(d.hr,n.wd))
    
    data.out1<-cbind(DOY2,HR2)
    data.out2<-cbind(DOY2,HR2)
    data.out3<-cbind(DOY2,HR2)
    data.out4<-cbind(DOY2,HR2)
    data.out5<-cbind(DOY2,HR2)
    
    for(j in 1:length(sel.var)){
      if(sel.var[j]!="CO2"&sel.var[j]!="TS"&sel.var[j]!="SWC"){
        # single culumn in original files
        for(i in 1:n.wd){
          temp.out1[,i]<-c(tapply(data.in[data.in$DOY2==i,sel.var[j]],data.in$HR2[data.in$DOY2==i],low.bd))
          temp.out2[,i]<-c(tapply(data.in[data.in$DOY2==i,sel.var[j]],data.in$HR2[data.in$DOY2==i],na.median))
          temp.out3[,i]<-c(tapply(data.in[data.in$DOY2==i,sel.var[j]],data.in$HR2[data.in$DOY2==i],upp.bd))
          temp.out4[,i]<-c(tapply(data.in[data.in$DOY2==i,sel.var[j]],data.in$HR2[data.in$DOY2==i],low.bd2))
          temp.out5[,i]<-c(tapply(data.in[data.in$DOY2==i,sel.var[j]],data.in$HR2[data.in$DOY2==i],upp.bd2))
        }  
      }else{
        for(i in 1:n.wd){
          temp.out1[,i]<-c(tapply(c(data.in[data.in$DOY2==i,paste(sel.var[j],"_1",sep="")],data.in[data.in$DOY2==i,paste(sel.var[j],"_2",sep="")]),
                                  c(data.in$HR2[data.in$DOY2==i],data.in$HR2[data.in$DOY2==i]),low.bd))
          temp.out2[,i]<-c(tapply(c(data.in[data.in$DOY2==i,paste(sel.var[j],"_1",sep="")],data.in[data.in$DOY2==i,paste(sel.var[j],"_2",sep="")]),
                                  c(data.in$HR2[data.in$DOY2==i],data.in$HR2[data.in$DOY2==i]),na.median))
          temp.out3[,i]<-c(tapply(c(data.in[data.in$DOY2==i,paste(sel.var[j],"_1",sep="")],data.in[data.in$DOY2==i,paste(sel.var[j],"_2",sep="")]),
                                  c(data.in$HR2[data.in$DOY2==i],data.in$HR2[data.in$DOY2==i]),upp.bd))
          temp.out4[,i]<-c(tapply(c(data.in[data.in$DOY2==i,paste(sel.var[j],"_1",sep="")],data.in[data.in$DOY2==i,paste(sel.var[j],"_2",sep="")]),
                                  c(data.in$HR2[data.in$DOY2==i],data.in$HR2[data.in$DOY2==i]),low.bd2))
          temp.out5[,i]<-c(tapply(c(data.in[data.in$DOY2==i,paste(sel.var[j],"_1",sep="")],data.in[data.in$DOY2==i,paste(sel.var[j],"_2",sep="")]),
                                  c(data.in$HR2[data.in$DOY2==i],data.in$HR2[data.in$DOY2==i]),upp.bd2))
        }
      }
      
      ## fill the gaps in each wd
      for(i in 1:n.wd){
        if(sum(!is.na(temp.out3[,i]))!=d.hr&sum(!is.na(temp.out3[,i]))>0.33*d.hr){
          temp.out1[,i]<-na.approx(temp.out1[,i],rule=2)
          temp.out2[,i]<-na.approx(temp.out2[,i],rule=2)
          temp.out3[,i]<-na.approx(temp.out3[,i],rule=2)
          temp.out4[,i]<-na.approx(temp.out4[,i],rule=2)
          temp.out5[,i]<-na.approx(temp.out5[,i],rule=2)
        }
      }
      
      out.name<-ifelse(paste(sel.var[j])!="NEE_PI"&paste(sel.var[j])!="VPD_PI",paste(sel.var[j]),substr(paste(sel.var[j]),start=1,stop=3))
      
      data.out1<-cbind(data.out1,tmp=c(temp.out1))
      colnames(data.out1)[which(colnames(data.out1)=="tmp")]<-out.name
      data.out2<-cbind(data.out2,tmp=c(temp.out2))
      colnames(data.out2)[which(colnames(data.out2)=="tmp")]<-out.name
      data.out3<-cbind(data.out3,tmp=c(temp.out3))
      colnames(data.out3)[which(colnames(data.out3)=="tmp")]<-out.name
      data.out4<-cbind(data.out4,tmp=c(temp.out4))
      colnames(data.out4)[which(colnames(data.out4)=="tmp")]<-out.name
      data.out5<-cbind(data.out5,tmp=c(temp.out5))
      colnames(data.out5)[which(colnames(data.out5)=="tmp")]<-out.name
    }
    
    png(filename=paste(path.out,case,"-%d.png",sep=""),width=11,height=9,units="in",pointsize=10,res=300)  
    par(mfrow=c(4,2),oma=c(2.5,0.5,2.5,0.5),mar=c(3.5,4.5,3.5,0.5))
    for(j in 3:ncol(data.out1)){
      if(sum(!is.na(data.out2[,j]))>0){
        rng<-range(c(data.out1[,j],data.out2[,j],data.out3[,j],data.out4[,j],data.out5[,j]),na.rm=T)
        plot(data.out2[,j],ylim=rng,ylab="",pch=16,cex=0.9,xaxt="n",xlab="")
        lines(data.out1[,j])
        lines(data.out3[,j])
        points(data.out1[,j],pch=16,cex=0.5)
        points(data.out3[,j],pch=16,cex=0.5)
        lines(data.out4[,j],col="darkgrey")
        lines(data.out5[,j],col="darkgrey")
        points(data.out4[,j],pch=16,cex=0.5,col="darkgrey")
        points(data.out5[,j],pch=16,cex=0.5,col="darkgrey")
        axis(side=3,at=c(0.5*d.hr+d.hr*c(0:(n.wd-1))),labels=c(0.5*l.wd+l.wd*c(0:(n.wd-1))))
        mtext(side=2,paste(colnames(data.out1)[j]),line=2.5,font=2)
        mtext(side=3,"DOY",line=0,outer=T)
        abline(v=d.hr*c(0:n.wd)+0.5,col="grey",lty=4)
        axis(side=1,at=seq(1,nrow(data.out1),by=d.hr/4)-0.5,labels=rep(c(0,6,12,18),times=nrow(data.out1)/d.hr),cex.axis=0.5)
        mtext(side=1,"HOUR",line=0,outer=T)
        mtext(side=3,paste("N=",floor(nrow(data.in)/17520),"year"),adj=0,line=0,outer=T)      
      }
    }
    dev.off()
    
    write.csv(data.out1,file=paste(path.out,case,"_LOWER1.csv",sep=""),quote=F,na="-9999",row.names=F)
    write.csv(data.out2,file=paste(path.out,case,"_MEDIAN.csv",sep=""),quote=F,na="-9999",row.names=F)
    write.csv(data.out3,file=paste(path.out,case,"_UPPER1.csv",sep=""),quote=F,na="-9999",row.names=F)
    write.csv(data.out4,file=paste(path.out,case,"_LOWER2.csv",sep=""),quote=F,na="-9999",row.names=F)
    write.csv(data.out5,file=paste(path.out,case,"_UPPER2.csv",sep=""),quote=F,na="-9999",row.names=F)
  }
  # update progress bar
  Sys.sleep(0.001)
  setTxtProgressBar(pb,k)
}

