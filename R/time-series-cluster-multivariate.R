rm(list = ls()) # clean working memory

require("pdc")
require("proxy")
require("dtw")
require("rgdal")
require("ape") ##  for displaying more appealing trees
require("dtwclust")


path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
path.in.root <- paste0(path, "diurnal-seasonality\\")
path.out.root <- paste0(path, "cluster-output\\")
RDir <- paste0(path, "R\\")
map.in <- paste0(path, "ecoregion\\")

source(paste0(RDir, "math_utility.R"))
source(paste0(RDir, "get_ecomap.R"))

###### Control parameter
ver <- "20210331"

len.ts <- 24 * 24  ## length of time series
n.grp <- c(10:30)  # searching window for the n of clusters/branches to keep

# target variables to run multi-variate clustering
target.var.ls <- list(c("SW_IN", "TA", "VPD", "SWC"),
                      c("NEE", "LE", "H", "USTAR"),
                      c("NEE", "LE", "H", "USTAR", "SW_IN", "TA", "VPD", "SWC"))
target.var.outname.ls <- c("ALL_MET", "ALL_FLUX", "ALL_FLUXMET")

path.in <- paste0(path.in.root, ver, "\\")
path.out <- paste0(path.out.root, ver, "\\")

###################################################################################################################
## extract full site list
full.ls <-
  read.csv(
    paste(path.out, "ALL_BASE_site_short_list2.csv", sep = ""),
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


## Work by variable groups
for (l1 in 1:length(target.var.ls)) {
  target.var <- target.var.ls[[l1]]
  
  src.ls <- sub("_MEDIAN.csv", "", src.list.in)
  src.ls <- src.ls[!duplicated(src.ls)]
  
  site.ls <- substr(src.ls, start = 1, stop = 6)
  res.ls <- substr(src.ls, start = 8, stop = 9)
  ver.ls <- substr(src.ls, start = 11, stop = nchar(src.ls))
  
  print(paste("################################################"))
  print(paste("######  ", target.var.outname.ls[l1], "  ##################"))
  
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
    if (sum(is.na(data.tmp)) != 0) {
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
  lat.ls <- as.numeric(full.ls[site.ls, "GRP_LOCATION.LOCATION_LAT"])
  lon.ls <- as.numeric(full.ls[site.ls, "GRP_LOCATION.LOCATION_LONG"])
  
  ################################################################################################################
  #### prepare full diurnal-seasonal time series
  
  num.ts <- length(src.ls) # number of time series
  num.dim <- length(target.var) # number of dimensions

  data.pre <- list() #array(dim = c(len.ts, num.ts, num.dim))
  data.mean <- NULL
  data.sd <- NULL
    
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
    
    if(i2 == 1){
      data.mean <- data.tmp
      data.sd <- data.tmp
    }else{
      data.mean <- data.mean + data.tmp
      data.sd <- rbind(data.sd,
                       data.tmp) 
    }
  }
  # network composite mean
  data.mean <- data.mean / length(src.ls)
  
  # for normalize time series across sites
  data.mean.mean <- matrix(rep(apply(data.mean, 2, mean), each = len.ts), 
                                      byrow = FALSE, nrow = len.ts) 
  data.mean.sd <- matrix(rep(apply(data.sd, 2, sd), each = len.ts), 
                                    byrow = FALSE, nrow = len.ts)

  ## normalize time series by network mean & sd  
  for(i3 in 1:length(data.pre)){
    data.pre[[i3]] <- (data.pre[[i3]] - data.mean.mean) / data.mean.sd
    
  }
  
  write.csv(
    data.pre,
    paste0(
      path.out,
      "AMF-diurnal-seasonal-cluster-",
     target.var.outname.ls[l1],
      "-composite_mean.csv"
    ),
    row.names = F
  )
  
  #############################################################################################################
  #### Clustering time series

  # save distance matrix for output
  distMatrix <- proxy::dist(data.pre,
                            method = "dtw_basic",
                            step.pattern = symmetric1)
  
  ## main clustering function
  hc1 <- dtwclust::tsclust(
    data.pre,
    type = "hierarchical",
    k = n.grp,
    seed = 1234,
    distance = "dtw_basic",
    control = hierarchical_control(method = "average")
  )
  
  ## look through different n of trees and their CVI
  names(hc1) <- paste0("k_", n.grp)
  cvi.est <- as.data.frame(t(sapply(hc1, cvi, type = "internal")))
  cvi.est <- cvi.est[, -which(colnames(cvi.est) == "SF")]
  
  for (i in 1:ncol(cvi.est)) {
    if (colnames(cvi.est)[i] %in% c("DB", "DBstar", "COP")) {
      ## minimize
      cvi.est[, i] <- 1 - cvi.est[, i] / max(cvi.est[, i])
    } else{
      ## maximize
      cvi.est[, i] <- cvi.est[, i] / max(cvi.est[, i])
    }
  }
  
  png(
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", 
          target.var.outname.ls[l1],
           "-ncluster_cvi.png"),
    width = 6,
    height = 4.5,
    units = "in",
    pointsize = 10,
    res = 300
  )
  plot(
    n.grp,
    cvi.est[, 1] / max(cvi.est[, 1]),
    type = "l",
    xlim = range(n.grp),
    ylim = c(0, 1.2),
    las = 1,
    xlab = "number of cluster",
    ylab = "normalized cluster validity index",
    lwd = 2,
    col = rainbow(6)[1],
    main = target.var
  )
  for (i in 2:ncol(cvi.est)) {
    lines(n.grp,
          cvi.est[, i] / max(cvi.est[, i]),
          col = rainbow(6)[i],
          lwd = 2)
  }
  legend(
    10,
    1.2,
    lty = 1,
    col = c("black", rainbow(6)),
    ncol = 4,
    legend = c("MEAN", colnames(cvi.est)),
    lwd = 1.5,
    bty = "n"
  )
  lines(n.grp, apply(cvi.est, 1, mean), lwd = 3)
  abline(v = n.grp[which(apply(cvi.est, 1, mean) == max(apply(cvi.est, 1, mean)))], lty = 4)
  dev.off()
  
  ## Choose optimal tree number based on CH (Calinski-Harabasz) index
  n.grp.opt <-
    n.grp[which(apply(cvi.est, 1, mean) == max(apply(cvi.est, 1, mean)))]
  
  hc <- dtwclust::tsclust(
    data.pre,
    type = "hierarchical",
    k = n.grp.opt,
    seed = 1234,
    distance = "dtw_basic",
    control = hierarchical_control(method = "average")
  )
  
  grp.ls <- cutree(hc, k = n.grp.opt)
  
  write.csv(
    as.matrix(distMatrix),
    paste0(path.out, "AMF-diurnal-seasonal-cluster-",
          target.var.outname.ls[l1],
           "-distMatrix.csv"),
    quote = T
  )
  
  # return cluster ID to full list
  full.ls <- data.frame(full.ls,
                        tmp = NA)
  
  for (i3 in 1:length(grp.ls)) {
    full.ls$tmp[which(full.ls$SITE_ID == paste(names(grp.ls)[i3]))] <-
      grp.ls[i3]
  }
  colnames(full.ls)[which(colnames(full.ls) == "tmp")] <-
    paste(target.var.outname.ls[l1], "_clust_group", sep = "")
  
  ## plot trees
  png(
    paste0(path.out, "AMF-diurnal-seasonal-cluster-",
          target.var.outname.ls[l1],
           "-tree.png"),
    width = 9,
    height = 9,
    units = "in",
    pointsize = 10,
    res = 300
  )
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 4.5, 4.5))
  plot(
    ape::as.phylo(hc),
    type = "fan",
    tip.color = rainbow(n.grp.opt)[grp.ls],
    label.offset = 0.5,
    #show.node.label=T,
    cex = 0.9
  )
  #text(0,0,)
  dev.off()
  
  print(paste("################################################"))
  print(paste("                                                "))
  
}

write.csv(
  full.ls,
  paste0(path.out, "ALL_BASE_site_short_list2.csv"),
  quote = T,
  row.names = F
)
