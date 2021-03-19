filter_physical_range<-function(data.in,
                                limit.ls,
                                basename_decode,
                                echo=T,
                                loose.filter=0.1    ## loosen the filtering criteria by +/- loose.filter*100%*full range
                                                    #  if NA, then don't loosen
                                ){
  
  ## re-order the data to match the order as in basename_decode
  var.order<-NULL
  for(i in 1:nrow(basename_decode)){
    var.order<-c(var.order,which(colnames(data.in)==paste(basename_decode$variable_names[i])))
  }
  data.in<-data.in[,var.order] 
  
  ## pre-filter by physical limits
  for(l in 1:(ncol(data.in))){
    
    if(length(which(limit.ls$Name==basename_decode$basename[l]))>0){
      
      # locate corresponding criteria via the parsed basename
      var.upp<-limit.ls$Max[which(limit.ls$Name==basename_decode$basename[l])]
      var.low<-limit.ls$Min[which(limit.ls$Name==basename_decode$basename[l])]
      
      # adjust +/- 10% for loose filtering
      if(!is.na(var.upp)&!is.na(var.low)&!is.na(loose.filter)){
        var.upp<-var.upp+0.1*abs(var.upp-var.low)
        var.low<-var.low-0.1*abs(var.upp-var.low)
      }
      
      if(!is.na(var.upp)){
        
        if(echo&length(which(data.in[,l]>var.upp))>0){
          print(paste("[Info]",basename_decode$variable_names[l],"filterd:",length(which(data.in[,l]>var.upp))))  
        }
        data.in[which(data.in[,l]>var.upp),l]<-NA
        
      }
      if(!is.na(var.low)){
        
        if(echo&length(which(data.in[,l]<var.low))>0){
          print(paste("[Info]",basename_decode$variable_names[l],"filterd:",length(which(data.in[,l]<var.low))))  
        }
        
        data.in[which(data.in[,l]<var.low),l]<-NA
      }
    }
  }
  
  return(data.in)
}