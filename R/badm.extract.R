###################################################################################################
badm.extract <- function(data.in, sel.grp) {
  
  if (length(which(data.in$VARIABLE_GROUP == sel.grp)) > 0) {
    
    data.in3 <- data.in[data.in$VARIABLE_GROUP == sel.grp, ]
    data.in3$VARIABLE <- factor(data.in3$VARIABLE)
    var.ls <- levels(data.in3$VARIABLE)
    
    data.in3$GROUP_ID <- as.numeric(data.in3$GROUP_ID)
    entry.ls <- tapply(data.in3$GROUP_ID, data.in3$GROUP_ID, mean)
    
    entry.ls.loc.i <- tapply(data.in3$SITE_ID,
                             data.in3$GROUP_ID,
                             function(x)
                               paste(x[1]))
    
    #entry.ls %in% data.out3$GROUP_ID
    
    data.out3 <- data.frame(
      GROUP_ID = tapply(data.in3$GROUP_ID,
                        data.in3$GROUP_ID,
                        function(x)
                          paste(x[1])),
      SITE_ID = tapply(data.in3$SITE_ID,
                       data.in3$GROUP_ID,
                       function(x)
                         paste(x[1])),
      stringsAsFactors = F
    )
    
    for (j in 1:length(var.ls)) {
      data.in3.tmp <-
        data.in3[data.in3$VARIABLE == paste(var.ls[j]), c("GROUP_ID", "DATAVALUE")]
      data.out3 <- merge.data.frame(data.out3,
                                    data.in3.tmp,
                                    by = "GROUP_ID", all = T)
      
      colnames(data.out3)[ncol(data.out3)] <- paste(var.ls[j])
    }
    
    data.out3 <- data.out3[order(data.out3$SITE_ID), ]
    
  } else{
    data.out3 <- NULL
    
  }
  
  return(data.out3)
}
