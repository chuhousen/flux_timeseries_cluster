get_full_list2 <- function(base.in = base.in) {
  list.tmp <- dir(base.in)
  list.tmp <-
    list.tmp[substr(list.tmp,
                    start = nchar(list.tmp) - 3,
                    stop = nchar(list.tmp)) == ".zip"]
  list.tmp <- list.tmp[substr(list.tmp, start = 24, stop = 30) == "FULLSET"]
  
  x.length <- nchar(list.tmp)
  x.ini <- as.vector(regexpr("FLX_", list.tmp))
  x.cut <- x.ini + 5
  site.id <- substr(list.tmp, start = x.ini + 4, stop = x.cut + 4)
  site.id <- site.id[!duplicated(site.id)]
  
  return(site.id)
}