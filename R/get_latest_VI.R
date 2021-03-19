get_latest_VI<-function(path.out){
  
  library(RCurl) 
  
  url<-"ftp://ftp.fluxdata.org/.ameriflux_downloads/measurement_height/"
  
  filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames = paste(strsplit(filenames, "\r*\n")[[1]], sep = "")

  if(sum(nchar(filenames)==35)==length(filenames)){
    
    time.ls<-as.numeric(substr(filenames,start=24,stop=31))
    var.info.ver<-time.ls[which(time.ls==max(time.ls))]
    
    latest<-filenames[which(time.ls==max(time.ls))]
    
    download.file(paste(url,latest,sep="/"),
                  paste(path.out,latest,sep=""))
    
    var.info<-read.csv(paste(path.out,latest,sep=""),header=T,na.strings=c(""))
    var.info$Height<-as.numeric(as.character(var.info$Height))
    
    ## clean for a bug in earlier version
    var.info$Variable<-gsub("_PI_PI","_PI",var.info$Variable)
    
  }else{
    
    print(paste("[Error] trouble locating the ltest VI file"))
    print(paste(filenames,collapse=" "))
    
    var.info<-NULL
    var.info.ver<-NULL
  }
  
  return(list(var.info,
              var.info.ver))
}