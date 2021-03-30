get_ustar_threshold <- function(data.work,
                                fc.loc,
                                sw_pot.loc,
                                ta.loc,
                                ustar.loc,
                                lat,
                                long,
                                d.hr.org){
    data.tmp <-
      data.work[, c(which(colnames(data.work) == "TIMESTAMP_END"),
                  fc.loc,
                  sw_pot.loc,
                  ta.loc,
                  ustar.loc)]
    colnames(data.tmp) <-
      c("TIMESTAMP_END", "NEE", "Rg", "Tair", "Ustar")
    
    data.tmp$TIMESTAMP_END <-
      strptime(data.tmp$TIMESTAMP_END, format = "%Y%m%d%H%M", tz = "UTC")
    
    ## check data availability for each year
    data.tmp2 <- na.omit(data.tmp)
    year.sum <- as.numeric(names(which(table(data.tmp2$TIMESTAMP_END$year + 1900) > d.hr.org * 100))) ## at least 100 days per year
    
    if(length(year.sum) > 0){
      
      ## work on any gap year
      year.check <-c(min(year.sum): max(year.sum))
      year.gap <- year.check[!year.check %in% year.sum]
      
      if(length(year.gap) > 0){
        year.seg <- c(year.sum[1]-1, year.gap, year.sum[length(year.sum)]+1)
        uStarTh <- NULL  
        uStarTh.single <- NULL
        
        ## work on each segment of consecutive years
        for(yy in seq_len(length(year.seg) - 1)) {
          # find valid year between gap years
          if (year.seg[yy + 1] - year.seg[yy] > 1) {
            data.tmp3 <-
              data.tmp[data.tmp$TIMESTAMP_END$year + 1900 > year.seg[yy] &
                         data.tmp$TIMESTAMP_END$year + 1900 < year.seg[yy + 1] ,]
            
            #Setting the basic info and the columns we are interested in analyzing.
            eddyC <-
              sEddyProc$new(
                'work',
                data.tmp3,
                c("NEE", "Tair", "Rg", "Ustar"),
                ColPOSIXTime = "TIMESTAMP_END",
                LatDeg = lat,
                LongDeg = long,
                DTS = d.hr.org
              )
            
            uStarTh.tmp <- eddyC$sEstUstarThreshold()$uStarTh
            uStarTh <- rbind.data.frame(uStarTh,
                                        uStarTh.tmp[uStarTh.tmp$aggregationMode != "single", ])
            uStarTh.single <- c(uStarTh.single, 
                                uStarTh.tmp$uStar[uStarTh.tmp$aggregationMode == "single"])
          }
        }
        
        ## recalculate single ustar threshold        
        uStarTh <- rbind.data.frame(uStarTh,
                                    data.frame(aggregationMode = "single",
                                               seasonYear = NA,
                                               season = NA,
                                               uStar = mean(uStarTh.single, na.rm = T)))
        
      }else{
        ## no gap year
        data.tmp3 <- data.tmp[data.tmp$TIMESTAMP_END$year + 1900 >= min(year.sum) & 
                               data.tmp$TIMESTAMP_END$year + 1900 <= max(year.sum) , ]  
        
        #Setting the basic info and the columns we are interested in analyzing.
        eddyC <-
          sEddyProc$new('work',
                        data.tmp3,
                        c("NEE", "Tair", "Rg", "Ustar"),
                        ColPOSIXTime = "TIMESTAMP_END",
                        LatDeg = lat,
                        LongDeg = long,
                        DTS = d.hr.org
          )
        
        uStarTh <- eddyC$sEstUstarThreshold()$uStarTh
      }

    }else{
      
      uStarTh <- NULL
      
    }

    return(uStarTh)
    
  }