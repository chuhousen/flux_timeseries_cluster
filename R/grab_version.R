
grab.version<-function(x){
  x.length<-nchar(x)
  x.ini<-as.vector(regexpr("BADM_",x))
  x.cut<-NULL
  for(i in 1:length(x)){
    x.cut<-c(x.cut,
             gregexpr("-",x[i])[[1]][length(gregexpr("-",x[i])[[1]])])
  }
  data.ver<-as.numeric(substr(x,start=x.ini+5,stop=x.cut-1))
  proc.ver<-as.numeric(substr(x,start=x.cut+1,stop=nchar(x)))
  return(list(data.ver,proc.ver))
}

latest.two<-function(x,y){
  tmp.xy<-as.numeric(paste(x,".",y,sep=""))
  two.to.grab<-sort(tmp.xy,decreasing=T)
  pre.version<-paste(x[which(tmp.xy==two.to.grab[2])],"-",y[which(tmp.xy==two.to.grab[2])],sep="")
  new.version<-paste(x[which(tmp.xy==two.to.grab[1])],"-",y[which(tmp.xy==two.to.grab[1])],sep="")
  return(c(new.version,pre.version))
}


