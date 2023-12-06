rm(list = ls()) # clean working memory

require("pdc")
require("proxy")
require("dtw")
#require("rgdal")
require("ape") ##  for displaying more appealing trees
require("dtwclust")
require("nnet") # for multinomial logistic model
require("DescTools") # for calculating Pseudo R-square
library(maps)

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
path.in.root <- paste0(path, "diurnal-seasonality\\")
path.out.root <- paste0(path, "cluster-output\\")
RDir <- paste0(path, "R\\")
map.in <- paste0(path, "ecoregion\\")

source(paste0(RDir, "math_utility.R"))
source(paste0(RDir, "get_ecomap.R"))
source(paste0(RDir, "cluster_color_panel2.R"))
source(paste0(RDir, "igbp_table.R"))
source(paste0(RDir, "ecoregion_table.R"))

###### Control parameter
ver <- "20231019-12h7d"
dist.ls <- c("dtw")#, "euclidean")
use.prescribed.col <- T

len.ts <- 52 * 12  ## length of time series
n.grp <- seq(30, 75, by = 1)  # searching window for the n of clusters/branches to keep 

# target variables to run multi-variate clustering
target.var.ls <- list(c("NETRAD", "TA", "VPD", "SWC"),
                      c("NEE", "LE", "H", "USTAR"),
                      c("NEE", "LE", "H", "USTAR", "NETRAD", "TA", "VPD", "SWC"))
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
    if (sum(is.na(data.tmp)) != 0 |
        (site.ls[i0] == "US-BZF" & "SWC" %in% target.var.ls[l1])) {
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
  
  ###### Loop through diff distance metric
  for(dd in 1: length(dist.ls)){
    
    target.dist <- dist.ls[dd]
    target.dist.get <- ifelse(target.dist == "dtw", "dtw_basic", 
                              ifelse(target.dist == "euclidean", "Euclidean",
                                     ifelse(target.dist == "gak", "gak", NA)))
    
    #############################################################################################################
    #### Clustering time series
    
    # save distance matrix for output
    # distMatrix <- proxy::dist(data.pre,
    #                           method = target.dist.get,
    #                           step.pattern = symmetric1)
    
    if(target.dist.get == "dtw_basic"){
      ## main clustering function
      hc1 <- dtwclust::tsclust(
        data.pre,
        type = "hierarchical",
        k = n.grp,
        seed = 1234,
        distance = target.dist.get,
        centroid = "pam"
        #control = hierarchical_control(method = "average")
      )
      
    }else{
      ## main clustering function
      hc1 <- dtwclust::tsclust(
        lapply(data.pre, as.vector),
        type = "hierarchical",
        k = n.grp,
        seed = 1234,
        distance = target.dist.get,
        centroid = "pam"
        #control = hierarchical_control(method = "average")
      )
    }
    
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
    
    write.csv(
      cvi.est,
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", 
             target.var.outname.ls[l1], "-CVI-", target.dist,".csv"),
      quote = F
    )
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", 
             target.var.outname.ls[l1],
             "-ncluster_cvi-", target.dist, ".png"),
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
    
    if(target.dist.get == "dtw_basic"){
      hc <- dtwclust::tsclust(
        data.pre,
        type = "hierarchical",
        k = n.grp.opt,
        seed = 1234,
        distance = target.dist.get,
        centroid = "pam"
        #control = hierarchical_control(method = "average")
      )
      
    }else{
      hc <- dtwclust::tsclust(
        lapply(data.pre, as.vector),
        type = "hierarchical",
        k = n.grp.opt,
        seed = 1234,
        distance = target.dist.get,
        centroid = "pam"
        #control = hierarchical_control(method = "average")
      )
    }  
    
    new.grp.ls <- grp.ls <- cutree(hc, k = n.grp.opt)
    
    ## reorder group ID to follow relative location in a tree
    sort.ls <- hc$labels[hc$order]
    new.grp.ls <- rep(NA, length(new.grp.ls))
    names(new.grp.ls) <- names(grp.ls)
    
    for(sss in 1:length(sort.ls)){
      if(is.na(new.grp.ls[which(names(new.grp.ls) == sort.ls[sss])])){
        
        new.grp.ls[which(names(new.grp.ls) %in% 
                           names(grp.ls[which(grp.ls == grp.ls[names(grp.ls) == sort.ls[sss]])]))] <-
          ifelse(is.finite(suppressWarnings(max(new.grp.ls, na.rm = T))),
                 max(new.grp.ls, na.rm = T) + 1,
                 1)
      }
    }
    
    ## assign group color
    ## assign group color
    if(use.prescribed.col){
      col_pan <- cluster_color_panel2(target.var = target.var.outname.ls[l1])
      
    }else{
      col_pan <- data.frame(r = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 1],
                            g = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 2],
                            b = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 3],
                            tmp = sort(unique(new.grp.ls)))
      colnames(col_pan)[4] <- target.var.outname.ls[l1]
    }

    # return cluster ID to full list
    full.ls <- data.frame(full.ls,
                          tmp = NA,
                          tmp2 = NA,
                          centroid = FALSE,
                          color = NA)
    col_pan_get <- rep(which(is.na(col_pan[, target.var.outname.ls[l1]])), nrow(full.ls))
    
    for (i3 in 1:length(new.grp.ls)) {
      full.ls$tmp[which(full.ls$SITE_ID == paste(names(new.grp.ls)[i3]))] <- new.grp.ls[i3]
      full.ls$centroid[which(full.ls$SITE_ID == paste(names(hc@centroids)[i3]))] <- TRUE
      
      if(new.grp.ls[i3] %in% col_pan[, target.var.outname.ls[l1]]){
        col_pan_get[i3] <- which(col_pan[, target.var.outname.ls[l1]] == new.grp.ls[i3])
        
      }else if(!is.na(new.grp.ls[i3])){
        col_pan_get[i3] <- which(col_pan[, target.var.outname.ls[l1]] == 99)
        
      }
      
      full.ls$color[which(full.ls$SITE_ID ==  paste(names(new.grp.ls)[i3]))] <- rgb(col_pan$r[col_pan_get[i3]],
                                                                                    col_pan$g[col_pan_get[i3]],
                                                                                    col_pan$b[col_pan_get[i3]], 
                                                                                    maxColorValue = 255)
      
    }
    
    for(i33 in 1:length(sort.ls)){
      full.ls$tmp2[which(full.ls$SITE_ID == sort.ls[i33])] <- i33
    }
    
    colnames(full.ls)[which(colnames(full.ls) == "tmp")] <-
      paste(target.var.outname.ls[l1], "_clust_group_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "tmp2")] <-
      paste(target.var.outname.ls[l1], "_group_order_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "centroid")] <-
      paste(target.var.outname.ls[l1], "_centroid_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "color")] <-
      paste(target.var.outname.ls[l1], "_color_", target.dist, sep = "")
    
    
    ## plot trees
    hc.plot.tree <- ape::as.phylo(hc)
    scale.facor <-
      (range(hc.plot.tree$edge.length)[2] - range(hc.plot.tree$edge.length)[1]) *
      0.01
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-", target.dist, ".png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "fan",
      tip.color = rgb(col_pan$r[col_pan_get],
                      col_pan$g[col_pan_get],
                      col_pan$b[col_pan_get], maxColorValue = 255),
      label.offset = scale.facor,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-", target.dist, "-ecoregion.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "fan",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "ecoregion.color"],
      label.offset = scale.facor,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    legend(
      "center",
      legend = ecoregion.color.table$plot.text[ecoregion.color.table$ecoregion %in% unique(full.ls$eco_L1)],
      fill = ecoregion.color.table$col[ecoregion.color.table$ecoregion %in% unique(full.ls$eco_L1)],
      border = NA,
      box.lty = 0,
      #bty = "n",
      cex = 1,
      ncol = 2,
      bg = rgb(255, 255, 255, 175, maxColorValue = 255)
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-", target.dist, "-igbp.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "fan",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "igbp.color"],
      label.offset = scale.facor,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    legend(
      "center",
      legend = igbp.color.table$plot.text[igbp.color.table$igbp %in% unique(full.ls$IGBP)],
      fill = igbp.color.table$col[igbp.color.table$igbp %in% unique(full.ls$IGBP)],
      border = NA,
      box.lty = 0,
      #bty = "n",
      cex = 1,
      ncol = 2,
      bg = rgb(255, 255, 255, 175, maxColorValue = 255)
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-radial-", target.dist, ".png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = rgb(col_pan$r[col_pan_get],
                      col_pan$g[col_pan_get],
                      col_pan$b[col_pan_get], maxColorValue = 255),
      label.offset = scale.facor * 0.0002,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-radial-", target.dist, "-ecoregion.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "ecoregion.color"],
      label.offset = scale.facor * 0.0002,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    legend(
      "center",
      legend = ecoregion.color.table$plot.text[ecoregion.color.table$ecoregion %in% unique(full.ls$eco_L1)],
      fill = ecoregion.color.table$col[ecoregion.color.table$ecoregion %in% unique(full.ls$eco_L1)],
      border = NA,
      box.lty = 0,
      #bty = "n",
      cex = 1,
      ncol = 2,
      bg = rgb(255, 255, 255, 175, maxColorValue = 255)
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var.outname.ls[l1],
             "-tree-radial-", target.dist, "-igbp.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "igbp.color"],
      label.offset = scale.facor * 0.0002,
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    legend(
      "center",
      legend = igbp.color.table$plot.text[igbp.color.table$igbp %in% unique(full.ls$IGBP)],
      fill = igbp.color.table$col[igbp.color.table$igbp %in% unique(full.ls$IGBP)],
      border = NA,
      box.lty = 0,
      #bty = "n",
      cex = 1,
      ncol = 2,
      bg = rgb(255, 255, 255, 175, maxColorValue = 255)
    )
    #text(0,0,)
    dev.off()
    
   
    ############################################################################
    ## Map 
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var.outname.ls[l1], "-map-", target.dist,".png"),
      height = 5,
      width = 4.6,
      res = 400,
      units = "in",
      pointsize = 10
    )
    par(
      fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0)
    )
    plot(
      0,
      0,
      type = "n",
      xlim = c(-170, -30),
      ylim = c(-60, 75),
      yaxt = "n",
      xaxt = "n",
      xlab = "",
      ylab = "",
      xaxs = "i",
      yaxs = "i",
      bty = "n"
    )
    maps::map(
      "world",
      fill = TRUE,
      border = "white",
      col = "grey80",
      bg = "white",
      add = T,
      xlim = c(-170, -30),
      ylim = c(-60, 75),
      mar = c(0, 0, 0, 0)
    )
    polygon(c(-60, -30, -30, -60),
            c(60, 60, 75, 75),
            border = "white",
            col= "white")
    polygon(c(-70, -30, -30, -70),
            c(75, 75, 90, 90),
            border = "white",
            col= "white")
    points(
      full.ls$LOCATION_LONG[!is.na(full.ls[,paste0(target.var.outname.ls[l1], "_clust_group_", target.dist)])],
      full.ls$LOCATION_LAT[!is.na(full.ls[,paste0(target.var.outname.ls[l1], "_clust_group_", target.dist)])],
      bg = full.ls[,paste0(target.var.outname.ls[l1], "_color_", target.dist)][!is.na(full.ls[,paste0(target.var.outname.ls[l1], "_clust_group_", target.dist)])],
      col = "white",
      cex = 0.75,
      pch = 21,
      lwd = 0.4
    )
    legend(
      -165,
      10,
      title = "Group",
      title.adj = 0,
      legend = col_pan[, 4],
      fill = rgb(col_pan$r,
                 col_pan$g,
                 col_pan$b, maxColorValue = 255),
      border = NA,
      bty = "n",
      cex = 0.75,
      ncol = 4
    )
    dev.off()
    
    
    
    ###### explore cross table btw cluster & IGBP, Eco region
    short.ls <-
      data.frame(
        group = as.factor(full.ls[, grep(paste0(target.var.outname.ls[l1], "_clust_group_", target.dist), colnames(full.ls))]),
        igbp = as.factor(full.ls$IGBP),
        ecoregion = as.factor(full.ls$eco_L1_name),
        eco_igbp = as.factor(paste0(full.ls$eco_L1_name, "_", full.ls$IGBP))
      )
    
    write.csv(table(short.ls$eco_igbp,
                    short.ls$group),
              file = paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var.outname.ls[l1], "-", target.dist, "eco_igbp.csv"),
              row.names = T) 
    
    ####### Logistic regression to explore the predictability of IGBP & Ecoregion
    sink(file = paste(
      path.out,
      "AMF-diurnal-seasonal-cluster-", target.var.outname.ls[l1], "-", target.dist,
      "-logistic-summary.txt",
      sep = ""
    ))
    
    print(paste("IGBP + Ecoregion"))
    multi_mo <- multinom(group ~  igbp + ecoregion, 
                         data = na.omit(short.ls),
                         MaxNWts = 10000000,
                         model = TRUE)
    #summary(multi_mo)
    multi_mo_R2 <- PseudoR2(multi_mo, which = c("Nagelkerke"))
    print(paste("Pseudo R2:", round(multi_mo_R2, digits = 3)))
    print(paste("                                                "))
    
    print(paste("IGBP only"))
    multi_mo1 <- multinom(group ~  igbp, 
                          data = na.omit(short.ls),
                          MaxNWts = 10000000,
                          model = TRUE)
    #summary(multi_mo1)
    multi_mo1_R2 <- PseudoR2(multi_mo1, which = c("Nagelkerke"))
    print(paste("Pseudo R2", round(multi_mo1_R2, digits = 3)))
    print(paste("                                                "))
    
    print(paste("Ecoregion only"))
    multi_mo2 <- multinom(group ~  ecoregion, 
                          data = na.omit(short.ls),
                          MaxNWts = 10000000,
                          model = TRUE)
    #summary(multi_mo2)
    multi_mo2_R2 <- PseudoR2(multi_mo2, which = c("Nagelkerke"))
    print(paste("Pseudo R2", round(multi_mo2_R2, digits = 3)))
    print(paste("                                                "))
    
    sink()
    
    
    
  }
  
  print(paste("################################################"))
  print(paste("                                                "))
  
}

write.csv(
  full.ls,
  paste0(path.out, "ALL_BASE_site_short_list3.csv"),
  quote = T,
  row.names = F
)
