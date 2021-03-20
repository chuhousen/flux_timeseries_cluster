get_ustar_threshold <- function(data.in,
                                fc.loc,
                                sw.loc,
                                ta.loc,
                                ustar.loc,
                                lat,
                                long
                                ){
    data.tmp <-
      data.in[, c(which(colnames(data.in) == "TIMESTAMP_END"),
                  fc.loc,
                  sw.loc,
                  ta.loc,
                  ustar.loc)]
    colnames(data.tmp) <-
      c("TIMESTAMP_END", "NEE", "Rg", "Tair", "Ustar")
    
    data.tmp$TIMESTAMP_END <-
      strptime(data.tmp$TIMESTAMP_END, format = "%Y%m%d%H%M", tz = "UTC")
    
    #Setting the basic info and the columns we are interested in analyzing.
    eddyC <-
      sEddyProc$new(
        'work',
        data.tmp,
        c("NEE", "Tair", "Rg", "Ustar"),
        ColPOSIXTime = "TIMESTAMP_END",
        LatDeg = lat,
        LongDeg = long
      )
    
    uStarTh <- eddyC$sEstUstarThreshold()$uStarTh
    
    return(uStarTh)
    
  }