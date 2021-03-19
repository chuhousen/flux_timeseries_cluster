rename_aggregate<-function(basename_decode,
                           data.in){
  
  #### I. handle the equal/naive aggregation, i.e., only renaming var
  #   equal: unique variable, no aggregation needed, 
  #             -> use base name 
  #             e.g., ignore any _1_1_1 or _1 if provided
  #   naive: already layer-aggregated & not unique, 
  #             -> use provided layer index
  
  ## re-order the data to match the order as in basename_decode
  var.order<-NULL
  for(i in 1:nrow(basename_decode)){
    var.order<-c(var.order,which(colnames(data.in)==paste(basename_decode$variable_names[i])))
  }
  data.work<-data.in[,var.order] 
  
  ## drop gap-filled ones
  if(sum(basename_decode$is_gapfill)>0){
    data.work<-data.work[,!basename_decode$is_gapfill]
    basename_decode.work<-basename_decode[!basename_decode$is_gapfill,]
  }else{
    basename_decode.work<-basename_decode
  }
  
  ## handle "equal/naive" aggregation/renaming
  equal.ls<-which(basename_decode.work$aggregate_method=="equal"|basename_decode.work$aggregate_method=="naive")
  colnames(data.work)[equal.ls]<-as.character(basename_decode.work$proposed_name)[equal.ls]
  
  basename_decode.work$variable_names<-as.character(basename_decode.work$variable_names)
  basename_decode.work$variable_names[equal.ls]<-basename_decode.work$proposed_name[equal.ls]
  
  return(list(data.work,
              basename_decode.work))
}