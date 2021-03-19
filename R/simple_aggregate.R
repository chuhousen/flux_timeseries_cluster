simple_aggregate<-function(basename_decode.work,
                           data.work,
                           exclude.ls=NULL,
                           #skip.horizontal.average=c(),
                           echo=T){
  
  #### II. handle the 1step/2step aggregation
  #   1step: already replicate-averaged & not unique, 
  #             -> need to aggregate to layer
  #   2step: quadruplet provided & not unique, 
  #             -> need 1 or 2 steps aggragation 
  
  ## known issue
  #     1. _2_1_1 may not be at the same height as _1_1_1
  #     2. import height info from MeasurementHeightInfo
  #     3. differentiate above-/below-soil var ??
  
  ## order the data to match the order as in basename_decode
  var.order<-NULL
  for(i in 1:nrow(basename_decode.work)){
    var.order<-c(var.order,which(colnames(data.work)==paste(basename_decode.work$variable_names[i])))
  }
  data.work<-data.work[,var.order] 
  
  ## locating/preparing  
  aggre.ls<-which(basename_decode.work$aggregate_method=="1step"|basename_decode.work$aggregate_method=="2step")
  
  if(length(aggre.ls)>0){
    
    ## holder for outputs
    basename_decode.out<-basename_decode.work[-aggre.ls,]
    data.out<-data.work[,-aggre.ls]
    data.out.org.n<-ncol(data.out)
    
    ## temporary files for working, with all variables/data that need aggregation
    basename_decode.work.tmp<-basename_decode.work[aggre.ls,]
    data.work.tmp<-data.work[,aggre.ls]
    
    ## list of variables (basename) that need aggregation
    var.ls<-unique(basename_decode.work.tmp$basename)
    
    if(!is.null(exclude.ls)){
      for(i1 in 1:length(exclude.ls)){ 
        var.ls<-var.ls[-which(var.ls==exclude.ls[i1])]  
      }
    }
    
    ## work through the variable list
    for(j in 1:length(var.ls)){
      
      ## holder for new variable names/data, per target basename variable
      add.name<-NULL
      add.data<-data.frame(TIMESTAMP_START=data.work$TIMESTAMP_START) 
      
      ## same as above, except these are optional only for 2step aggregation  
      add.name.tmp<-NULL
      add.data.tmp<-data.frame(TIMESTAMP_START=data.work$TIMESTAMP_START) 
      
      ## TEmporary working files for specific target basename variable
      var.ls1<-which(basename_decode.work.tmp$basename==paste(var.ls[j]))
      basename_decode.work.tmp1<-basename_decode.work.tmp[var.ls1,]
      data.work.tmp1<-data.work.tmp[,var.ls1]
      
      #### only accept either all 1step or all 2step per basename variable, exception case handled by manual aggregation
      if(length(unique(basename_decode.work.tmp1$aggregate_method))==1){ 
        
        ### handle replicate averaging, only excute if 2step         
        if(unique(basename_decode.work.tmp1$aggregate_method)=="2step"){
          
          ## work through each _H_V location
          pos_id2<-paste("_",basename_decode.work.tmp1$H_index,"_",basename_decode.work.tmp1$V_index,sep="")
          pos_id2.ls<-table(pos_id2)
          
          for(j2 in 1:length(pos_id2.ls)){
            pos_id2.loc<-which(pos_id2==names(pos_id2.ls)[j2])
            
            ## aggregate from _H_V_R to _H_V_A, simply equal for 1-to-1 aggregation, mean for more-to-1 aggregation
            add.name.tmp<-c(add.name.tmp,
                            paste(basename_decode.work.tmp1$basename[1],names(pos_id2.ls)[j2],"_A",sep=""))  
            if(pos_id2.ls[j2]==1){
              add.data.tmp<-cbind.data.frame(add.data.tmp,
                                             tmp=data.work.tmp1[,pos_id2.loc])
            }else{
              add.data.tmp<-cbind.data.frame(add.data.tmp,
                                             tmp=apply(data.work.tmp1[,pos_id2.loc],1,na.mean))
            }
            colnames(add.data.tmp)[which(colnames(add.data.tmp)=="tmp")]<-paste(add.name.tmp[length(add.name.tmp)])
          }
          
          ## replace with replicate-aggregated data/basename_decode
          data.work.tmp1<-as.data.frame(add.data.tmp[,-1])
          basename_decode.work.tmp1<-basename_pharse(var.name=add.name.tmp,
                                                     FP.ls=FP.ls)
        }
        
        ### handle horizontal averaging, excute for both 1step/2step 
        pos_id<-paste(basename_decode.work.tmp1$V_index)
        pos_id.ls<-table(pos_id)
        
        ## work through _V location
        for(j1 in 1:length(pos_id.ls)){
          
          pos_id.loc<-which(pos_id==names(pos_id.ls)[j1])
          
          ## aggregate from _H_V_A to _#, simply equal for 1-to-1 aggregation, mean for more-to-1 aggregation
          add.name<-c(add.name,
                      paste(basename_decode.work.tmp1$basename[1],"_",names(pos_id.ls)[j1],sep=""))  
          
          if(pos_id.ls[j1]==1){
            add.data<-cbind.data.frame(add.data,
                                       tmp=data.work.tmp1[,pos_id.loc])
          }else{
            add.data<-cbind.data.frame(add.data,
                                       tmp=apply(data.work.tmp1[,pos_id.loc],1,na.mean))
          }
          colnames(add.data)[which(colnames(add.data)=="tmp")]<-paste(add.name[length(add.name)])
        }
        
        ## return agregated ones
        basename_decode.out<-rbind.data.frame(basename_decode.out,
                                              basename_pharse(var.name=add.name,
                                                              FP.ls=FP.ls)) 
        data.out.name<-c(colnames(data.out),
                         colnames(add.data)[-1])
        data.out<-cbind.data.frame(data.out,
                                   add.data[,-1])
        colnames(data.out)<-data.out.name
        
      }else{
        print(paste("[Warning]",var.ls[j],"has mixed 1step/2step"))
      }
    }
    
    if(echo){
      print(paste("[Info] Aggregate from:",paste(basename_decode.work$variable_names[aggre.ls],collapse=" ")))
      print(paste("[Info] To:",paste(colnames(data.out)[c((data.out.org.n+1):ncol(data.out))],collapse=" ")))
    }
    
  }else{
    if(echo) print("[Info] No 1step/2step aggrgation")    
    
    ## holder for outputs
    basename_decode.out<-basename_decode.work
    data.out<-data.work
    
  }
  
  return(list(data.out,
              basename_decode.out))
}