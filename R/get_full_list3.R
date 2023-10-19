get_full_list3 <- function(base.in = base.in,
                           target = "FULLSET") {
  list.tmp <- dir(base.in)
  list.tmp <-
    list.tmp[substr(list.tmp,
                    start = nchar(list.tmp) - 3,
                    stop = nchar(list.tmp)) == ".zip"]
  
  list.tmp <- list.tmp[grep(target, list.tmp)]
  
  #x.length <- nchar(list.tmp)
  list.tmp <- gsub("FLX_", "", list.tmp)
  list.tmp <- gsub("AMF_", "", list.tmp)
  
  x.ini <- 1
  x.cut <- x.ini + 5
  
  site.id <- substr(list.tmp, start = x.ini, stop = x.cut)
  site.id <- site.id[!duplicated(site.id)]
  
  return(site.id)
}