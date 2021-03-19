BIF_working<-function(data.in,
                      sel.grp){
  
  names(data.in)<-c("SITE_ID","GROUP_ID","VARIABLE_GROUP","VARIABLE","DATAVALUE")
  
  data.in$GROUP_ID<-as.character(data.in$GROUP_ID)
  data.in$VARIABLE_GROUP<-as.character(data.in$VARIABLE_GROUP)
  #table(data.in$VARIABLE_GROUP)
  data.in$DATAVALUE<-as.character(data.in$DATAVALUE)
  data.in$SITE_ID<-as.character(data.in$SITE_ID)
  
  badm.extract<-function(data.in,sel.grp){
    data.in3<-data.in[data.in$VARIABLE_GROUP==sel.grp,]
    data.in3$VARIABLE<-factor(data.in3$VARIABLE)
    var.ls<-levels(data.in3$VARIABLE)
    
    data.in3$GROUP_ID<-as.numeric(data.in3$GROUP_ID)
    entry.ls<-tapply(data.in3$GROUP_ID,data.in3$GROUP_ID,mean)
    
    data.out3<-data.frame(GROUP_ID=entry.ls,
                          SITE_ID=NA)
    
    tmp<-matrix(nrow=length(entry.ls),ncol=length(var.ls))
    data.out3<-data.frame(data.out3,tmp)
    names(data.out3)[which(names(data.out3)=="X1"):ncol(data.out3)]<-paste(var.ls)
    
    for(i in 1:length(entry.ls)){
      data.out3$SITE_ID[i]<-data.in3$SITE_ID[data.in3$GROUP_ID==entry.ls[i]][1]
      for(j in 1:length(var.ls)){
        if(length(which(data.in3$GROUP_ID==entry.ls[i]&data.in3$VARIABLE==paste(var.ls[j])))>0){
          data.out3[i,paste(var.ls[j])]<-data.in3$DATAVALUE[data.in3$GROUP_ID==entry.ls[i]&data.in3$VARIABLE==paste(var.ls[j])]  
        }
      }
    }
    data.out3<-data.out3[order(data.out3$SITE_ID),]
    return(data.out3)
  }
  
  data.out1.1<-badm.extract(data.in=data.in,sel.grp=sel.grp)
  
  return(data.out1.1)
}