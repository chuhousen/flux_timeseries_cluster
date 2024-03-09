rm(list = ls()) # clean working memory

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
path.in.root <- paste0(path, "diurnal-seasonality\\")
path.out.root <- paste0(path, "cluster-output\\")
RDir <- paste0(path, "R\\")
map.in <- paste0(path, "ecoregion\\")

source(paste0(RDir, "math_utility.R"))
source(paste0(RDir, "get_ecomap.R"))
source(paste0(RDir, "cluster_color_panel2.R"))
source(paste0(RDir, "ecoregion_table.R"))
source(paste0(RDir, "igbp_table.R"))

###### Control parameter
ver <- "20231019-12h7d"
dist.ls <- c("dtw")#, "euclidean")
use.prescribed.col <- T

len.ts <- 52 * 12  ## length of time series
len.ts.tck <- 24
# 24 * 24 for 20231019-24h15d
# 52 * 12 for 20231019-12h7d
# 73 * 8  for 20231019-8h5d

# for plotting
target.eco.ls <- c("Eastern Temperate Forests",
                   "Great Plains",
                   "Mediterranean California",
                   "Northern Forests",
                   "Northwestern Forested Mountains",
                   "Noth american Deserts")

# target variables to run uni-variate clustering
target.var.outname.ls <- c("ALL_FLUX")
target.var.ls <- list(c("NEE", "LE", "H", "USTAR"))
target.lab.ls <- list(expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
                      expression(LE~'('*W~m^{-2}*')'),
                      expression(H~'('*W~m^{-2}*')'),
                      expression(USTAR~'('*m~s^{-1}*')')
                      )
target.rng.ls <- list(seq(-60, 40, length.out = 6),
                      seq(-130, 520, length.out = 6),
                      seq(-130, 520, length.out = 6),
                      seq(0, 1.6, length.out = 5))

path.in <- paste0(path.in.root, ver, "\\")

if (!dir.exists(paste(path.out.root, ver, sep = "")))
  dir.create(paste(path.out.root, ver, sep = ""))
path.out <- paste(path.out.root, ver, "\\", sep = "")

###################################################################################################################
## extract full site list
full.ls <-
  read.csv(
    paste(path.out, "ALL_BASE_site_short_list3.csv", sep = ""),
    na = "-9999",
    header = T,
    stringsAsFactors = F
  )
full.ls <- full.ls[!is.na(full.ls$SITE_ID), ]
rownames(full.ls) <- paste(full.ls$SITE_ID)

## find a list of pre-processed diurnal-seasonal files
src.list.in <- list.files(path.in)
src.list.in <- src.list.in[which(substr(
  src.list.in,
  start = nchar(src.list.in) - 3,
  stop = nchar(src.list.in)
) == ".csv")]
src.list.in <- src.list.in[which(substr(
  src.list.in,
  start = nchar(src.list.in) - 9,
  stop = nchar(src.list.in) - 4
) == "MEDIAN")]

src.list.in <- src.list.in[which(substr(src.list.in,
                                        start = 1,
                                        stop = 6) %in% full.ls$SITE_ID)]

### Drop problematic cases
## drop.ls<-c("BR-CST","CA-SCB","US-Elm","US-Esm","US-Hn2","US-NC2",)

## add ecoregion/IGBP color panels
if(!"igbp.color" %in% colnames(full.ls)){
  full.ls$igbp.color <- NA
  full.ls$igbp.color2 <- NA
  full.ls$ecoregion.color <- NA
  full.ls$ecoregion.color2 <- NA
  for (i in 1:nrow(full.ls)) {
    full.ls$igbp.color[i] <-
      igbp.color.table$col[which(grepl(paste(full.ls$IGBP[i]),
                                       igbp.color.table$igbp))]
    full.ls$igbp.color2[i] <-
      igbp.color.table$col2[which(grepl(paste(full.ls$IGBP[i]),
                                        igbp.color.table$igbp))]
    
    if (!is.na(full.ls$eco_L1[i]) & full.ls$eco_L1[i] != "NA") {
      full.ls$ecoregion.color[i] <-
        ecoregion.color.table$col[which(grepl(
          paste(full.ls$eco_L1[i]),
          ecoregion.color.table$ecoregion
        ))]
      full.ls$ecoregion.color2[i] <-
        ecoregion.color.table$col2[which(grepl(
          paste(full.ls$eco_L1[i]),
          ecoregion.color.table$ecoregion
        ))]
    } else {
      full.ls$ecoregion.color[i] <- "gray"
      full.ls$ecoregion.color2[i] <- "gray"
    }
  }
}

## Work by variables
for (l1 in 1:length(target.var.ls)) {
  target.var <- target.var.ls[[1]]
  
  src.ls <- sub("_MEDIAN.csv", "", src.list.in)
  src.ls <- src.ls[!duplicated(src.ls)]
  
  site.ls <- substr(src.ls, start = 1, stop = 6)
  res.ls <- substr(src.ls, start = 8, stop = 9)
  ver.ls <- substr(src.ls, start = 11, stop = nchar(src.ls))
  
  if(sum(duplicated(site.ls)) > 0)
    stop(paste(
      "[error] multiple composite time series files for following",
      site.ls[which(duplicated(site.ls))]
    ))
  
  #########################################################################################################
  #### work on prepare site list
  
  ## read in pre-processed diurnal-seasonal files
  drop.case <- NULL
  for (i0 in 1:length(src.ls)) {
    data.tmp <-
      read.csv(
        paste(path.in, src.ls[i0], "_MEDIAN.csv", sep = ""),
        header = T,
        na.strings = "-9999",
        stringsAsFactors = F
      )
    data.tmp <- data.tmp[, target.var]
    
    ## check whether missing values
    # drop SWC from US-BZF as the gap-filling is problematic
    if (sum(is.na(data.tmp)) != 0 |
        (site.ls[i0] == "US-BZF" & "SWC" %in% target.var.ls[[1]])) {
      drop.case <- c(drop.case, site.ls[i0])
    }
  }
  
  print(paste(
    "[Info]",
    length(src.ls) - length(drop.case),
    "has sufficient data"
  ))
  
  if (length(drop.case) > 0) {
    for (i1 in 1:length(drop.case)) {
      src.ls <- src.ls[site.ls != drop.case[i1]]
      ver.ls <- ver.ls[site.ls != drop.case[i1]]
      res.ls <- res.ls[site.ls != drop.case[i1]]
      site.ls <- site.ls[site.ls != drop.case[i1]]
    }
    print(paste("[Info]", "drop", length(drop.case), "sites because of gaps"))
    print(paste("[Info]", "drop", paste(drop.case, collapse = " ")))
  }
  lat.ls <- as.numeric(full.ls[site.ls, "LOCATION_LAT"])
  lon.ls <- as.numeric(full.ls[site.ls, "LOCATION_LONG"])
  
  full.ls <- full.ls[which(full.ls$SITE_ID %in% site.ls),]
  
  ################################################################################################################
  #### prepare full diurnal-seasonal time series
  
  num.ts <- length(src.ls) # number of time series
  num.dim <- length(target.var) # number of dimensions

  data.pre <- list() #array(dim = c(len.ts, num.ts, num.dim))
  
  for (i2 in 1:length(src.ls)) {
    data.tmp <- read.csv(
      paste0(path.in, src.ls[i2], "_MEDIAN.csv"),
      header = T,
      na.strings = "-9999",
      stringsAsFactors = F
    )
    data.tmp <- data.matrix(data.tmp[, target.var])
    data.pre[[i2]] <- data.tmp
    names(data.pre)[i2] <- paste(site.ls[i2])
  }
  
  ###### Loop through ecoregions
  for(dd in 1: length(target.eco.ls)){
    
    print(paste("################################################"))
    print(paste("######  ", target.eco.ls[dd], "  ##################"))
    
    sub.ls <- full.ls[full.ls$eco_L1_name == target.eco.ls[dd],]
    igbp.ls <- unique(sub.ls$IGBP)
    grp.ls <- as.numeric(table(sub.ls$ALL_FLUX_clust_group_dtw))
    grp.name.ls <- names(table(sub.ls$ALL_FLUX_clust_group_dtw))
    data.pre.sub <- data.pre[which(full.ls$eco_L1_name == target.eco.ls[dd])]

    ############################################################################
    plot.n <- ifelse(length(grp.ls) > 14, 2, 1)
    
    for(gg in 1:plot.n){
      
      grp.ls.s1 <- grp.ls[c((1 + 14 * (gg - 1)) : min(c((14 + 14 * (gg - 1)),length(grp.ls))))]
      grp.name.ls.s1 <- grp.name.ls[c((1 + 14 * (gg - 1)) : min(c((14 + 14 * (gg - 1)),length(grp.ls))))]
      
      sub.ls.s1 <- sub.ls[sub.ls$ALL_FLUX_clust_group_dtw %in% grp.name.ls.s1,]
      data.pre.sub.s1 <- data.pre.sub[which(sub.ls$ALL_FLUX_clust_group_dtw %in% grp.name.ls.s1)]
      
      #### plot by ALL_FLUX groups
      png(
        paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.eco.ls[dd], "-time-series-ALLFLUX",
               "-", gg,".png"),
        width = 7,
        height = 0.5 + length(grp.ls.s1) * 5.5 / 9,
        units = "in",
        pointsize = 9,
        res = 300
      )
      
      par(
        mfrow = c(length(grp.ls.s1), 4),
        mar = c(0.2, 3, 0.2, 0.1),
        oma = c(3.5, 0, 1.5, 0.5)
      )
      
      for (j2 in 1:length(grp.ls.s1)) {
        for(j3 in 1:4){
          
          plot(
            0,
            0,
            xlim = c(0, len.ts),
            ylim = c(target.rng.ls[[j3]][1],
                     target.rng.ls[[j3]][length(target.rng.ls[[j3]])]),
            type = "n",
            xlab = "",
            ylab = "",
            xaxt = "n",
            yaxt = "n"
          )
          
         #if(length(which(sub.ls.s1$ALL_FLUX_clust_group_dtw == grp.name.ls.s1[j2])) > 1){
          data.pre.sub.s2 <- data.pre.sub.s1[which(sub.ls.s1$ALL_FLUX_clust_group_dtw == grp.name.ls.s1[j2])]  
          sub.ls.s2 <- sub.ls.s1[which(sub.ls.s1$ALL_FLUX_clust_group_dtw == grp.name.ls.s1[j2]),]

          mean.tmp <- NULL
          for (j1 in 1:length(data.pre.sub.s2)) {
            lines(data.pre.sub.s2[[j1]][,j3],
                  col = "darkgrey",
                  lwd = 0.5)
            if(j1 == 1){
              mean.tmp <- data.pre.sub.s2[[j1]][,j3]  
            }else{
              mean.tmp <- mean.tmp + data.pre.sub.s2[[j1]][,j3]    
            }
          }
          mean.tmp <- mean.tmp / length(data.pre.sub.s2)
          if(length(data.pre.sub.s2) >= 5){
            lines(mean.tmp,
                  col = sub.ls.s1$ALL_FLUX_color_dtw[which(sub.ls.s1$ALL_FLUX_clust_group_dtw == grp.name.ls.s1[j2])][1],
                  lwd = 0.75)  
          }
          
          axis(2, target.rng.ls[[j3]],las = 2)
          
          if(j3 == 1) {
            text(
              0,
              target.rng.ls[[j3]][length(target.rng.ls[[j3]])] - 0.1 * (target.rng.ls[[j3]][length(target.rng.ls[[j3]])] - target.rng.ls[[j3]][1]),
              paste0("Group:", grp.name.ls.s1[j2], " (", length(data.pre.sub.s2) , ")"),
              adj = c(0, 1),
              cex = 0.9, 
              font =2
            )
            
            if(length(data.pre.sub.s2) <= 12){
              legend(
                "bottomright",
                paste0(names(data.pre.sub.s2), " (", sub.ls.s2$IGBP, ")"),
                col = "darkgrey",
                cex = 0.7,
                bty = "n",
                ncol = ifelse(length(data.pre.sub.s2) > 6, 2, 1)
              )
            }
          }
          
          if(j2 == length(grp.ls.s1)){
            axis(1, at = seq(0, len.ts, by = len.ts.tck), labels = FALSE)
          }
        }
      }
      mtext(
        side = 1,
        "Hour / Window",
        outer = T,
        line = 1.5,
        cex = 1.2
      )
      mtext(
        side = 3,
        target.lab.ls[[1]],
        outer = T,
        line = 0,
        cex = 0.75,
        adj = c(0.1)
      )
      mtext(
        side = 3,
        target.lab.ls[[2]],
        outer = T,
        line = 0,
        cex = 0.75,
        adj = c(0.4)
      )
      mtext(
        side = 3,
        target.lab.ls[[3]],
        outer = T,
        line = 0,
        cex = 0.75,
        adj = c(0.65)
      )
      mtext(
        side = 3,
        target.lab.ls[[4]],
        outer = T,
        line = 0,
        cex = 0.75,
        adj = c(0.95)
      )
      
      
      dev.off()
      
    }
    
    
    ## plot grouped diurnal-seasonal plots
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.eco.ls[dd], "-time-series.png"),
      width = 7,
      height = 0.5 + length(igbp.ls) * 6 / 9,
      units = "in",
      pointsize = 9,
      res = 300
    )
    
    par(
      mfrow = c(length(igbp.ls), 4),
      mar = c(0.2, 3, 0.2, 0.1),
      oma = c(3.5, 0, 1.5, 0.5)
    )
    
    for (j2 in 1:length(igbp.ls)) {
      for(j3 in 1:4){
        
        plot(
          0,
          0,
          xlim = c(0, len.ts),
          ylim = c(target.rng.ls[[j3]][1],
                   target.rng.ls[[j3]][length(target.rng.ls[[j3]])]),
          type = "n",
          xlab = "",
          ylab = "",
          xaxt = "n",
          yaxt = "n"
        )
        
        data.pre.sub2 <- data.pre.sub[which(sub.ls$IGBP == igbp.ls[j2])]
        
        mean.tmp <- NULL
        for (j1 in 1:length(data.pre.sub2)) {
          lines(data.pre.sub2[[j1]][,j3],
                col = "darkgrey",
                lwd = 0.5)
          if(j1 == 1){
            mean.tmp <- data.pre.sub2[[j1]][,j3]  
          }else{
            mean.tmp <- mean.tmp + data.pre.sub2[[j1]][,j3]    
          }
        }
        mean.tmp <- mean.tmp / length(data.pre.sub2)
        if(length(data.pre.sub2) > 5){
          lines(mean.tmp,
                col = "black",
                lwd = 0.75)  
        }
        
        axis(2, target.rng.ls[[j3]],las = 2)
        
        if(j3 == 1) {
          text(
            0,
            target.rng.ls[[j3]][length(target.rng.ls[[j3]])] - 0.1 * (target.rng.ls[[j3]][length(target.rng.ls[[j3]])] - target.rng.ls[[j3]][1]),
            paste0(igbp.ls[j2], " (", length(data.pre.sub2) , ")"),
            adj = c(0, 1)
          )
          
          if(length(data.pre.sub2) <= 5){
            legend(
              "bottomright",
              paste(names(data.pre.sub2)),
              col = "darkgrey",
              cex = 0.9,
              bty = "n"
            )
          }
        }

        if(j2 == length(igbp.ls)){
          axis(1, at = seq(0, len.ts, by = len.ts.tck), labels = FALSE)
        }
      }
    }
    mtext(
      side = 1,
      "Hour / Window",
      outer = T,
      line = 1.5,
      cex = 1.2
    )
    mtext(
      side = 3,
      target.lab.ls[[1]],
      outer = T,
      line = 0,
      cex = 0.75,
      adj = c(0.1)
    )
    mtext(
      side = 3,
      target.lab.ls[[2]],
      outer = T,
      line = 0,
      cex = 0.75,
      adj = c(0.4)
    )
    mtext(
      side = 3,
      target.lab.ls[[3]],
      outer = T,
      line = 0,
      cex = 0.75,
      adj = c(0.65)
    )
    mtext(
      side = 3,
      target.lab.ls[[4]],
      outer = T,
      line = 0,
      cex = 0.75,
      adj = c(0.95)
    )
    
    
    dev.off()
  
  print(paste("################################################"))
  print(paste("                                                "))
  
  }
}

