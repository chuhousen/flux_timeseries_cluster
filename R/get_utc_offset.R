get_utc_offset <- function(site,
                           badm.path,
                           badm.file){
  
  badm.tmp.ls <-
    dir(badm.path)[which(grepl(badm.file, dir(badm.path)))]
  
  data.in <- read_excel(
    paste0(badm.path, badm.tmp.ls),
    sheet = 1,
    col_names = T,
    na = "-9999"
  )
  
  data.in$GROUP_ID <- as.character(data.in$GROUP_ID)
  data.in$VARIABLE_GROUP <- as.character(data.in$VARIABLE_GROUP)
  data.in$DATAVALUE <- as.character(data.in$DATAVALUE)
  data.in$SITE_ID <- as.character(data.in$SITE_ID)
  
  data.out <- badm.extract(data.in = data.in, sel.grp = "GRP_UTC_OFFSET")
 
  data.out <- data.out[!duplicated(data.out$SITE_ID), ]
    
  return(data.out) 
}