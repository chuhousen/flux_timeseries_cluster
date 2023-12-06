rm(list = ls()) # clean working memory

require("pdc")
require("proxy")
require("dtw")
require("sp")
#require("rgdal")
require("terra")
require("ape") ##  for displaying more appealing trees
require("dtwclust")
#require("gmodels") # for cross table
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

# searching window for the n of clusters/branches to keep 
n.grp <- seq(30, 75, by = 1)  
#n.grp2 <- seq(31, 75, by = 1)  
        #c(10:30) for cluster ms

## original colnames from diurnal-seasonal files, revising the colname
sel.var <- c("FC", "SC", "NEE", "LE", "H",
             #"FCH4", "SW_IN", "WTD", "P",
             "USTAR", "NETRAD", "TA", "VPD", 
             "SWC")

drop.ls<-c("FC", "SC")

# target variables to run uni-variate clustering
target.var.ls <- c("NEE", "LE", "H", "USTAR",
                   "NETRAD", "TA", "VPD", "SWC"
                   #"FCH4", "P", "WTD"
                   )
target.lab.ls <- list(expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
                      expression(LE~'('*W~m^{-2}*')'),
                      expression(H~'('*W~m^{-2}*')'),
                      expression(USTAR~'('*m~s^{-1}*')'),
                      expression(NETRAD~'('*W~m^{-2}*')'),
                      expression(TA~'('*degree~C*')'),
                      expression(VPD~'('*kPa*')'),
                      expression(SWC~'('*'%'*')')
                      #expression(FCH4~'('*nmole~m^{-2}~s^{-1}*')'),
                      #expression(P~'('*mm*')'),
                      #expression(WTD~'('*m*')')
                      )
target.rng.ls <- list(seq(-60, 40, length.out = 6),
                      seq(-130, 520, length.out = 6),
                      seq(-130, 520, length.out = 6),
                      seq(0, 1.6, length.out = 5),
                      seq(-300, 900, length.out = 5),
                      seq(-40, 40, length.out = 5),
                      seq(0, 60, length.out = 5),
                      seq(0, 100, length.out = 5)
                      #seq(-100, 500, length.out = 6),
                      #seq(0, 1500, length.out = 6),
                      #seq(-15, 15, length.out = 5)
                      )

path.in <- paste0(path.in.root, ver, "\\")

if (!dir.exists(paste(path.out.root, ver, sep = "")))
  dir.create(paste(path.out.root, ver, sep = ""))
path.out <- paste(path.out.root, ver, "\\", sep = "")

###################################################################################################################
## extract full site list
full.ls <-
  read.csv(
    paste(path.in, "ALL_BASE_site_short_list.csv", sep = ""),
    na = "-9999",
    header = T,
    stringsAsFactors = F
  )
full.ls <- full.ls[!is.na(full.ls$SITE_ID), ]
rownames(full.ls) <- paste(full.ls$SITE_ID)

for (l in 1:length(sel.var)) {
  if (sel.var[l] %in% drop.ls) {
    full.ls <- full.ls[, -which(colnames(full.ls) == sel.var[l])]
    full.ls <-
      full.ls[, -which(colnames(full.ls) == paste0(sel.var[l], "_F"))]
  } else{
    colnames(full.ls)[which(colnames(full.ls) == sel.var[l])] <-
      paste0(sel.var[l], "_availability")
    colnames(full.ls)[which(colnames(full.ls) == paste0(sel.var[l], "_F"))] <-
      paste0(sel.var[l], "_availability_after_filling")
  }
}


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


### read in Ecoregion map
ecomap1 <- terra::vect(
  paste(map.in, "sa_eco_l3", sep = ""),
  layer = "sa_eco_l3"
)
ecomap3 <- terra::vect(
  paste(map.in, "NA_CEC_Eco_Level3", sep = ""),
  layer = "NA_CEC_Eco_Level3"
)

ecomap1 <- terra::project(ecomap1,
                          "+proj=longlat +datum=WGS84")
ecomap3 <- terra::project(ecomap3,
                          "+proj=longlat +datum=WGS84")

### work on extraction of ecoregion for each site
pts <- terra::vect(cbind(longitude = full.ls$LOCATION_LONG,
                         latitude = full.ls$LOCATION_LAT),
                   crs = "+proj=longlat +datum=WGS84")
full.ls <- cbind.data.frame(full.ls,
                            terra::extract(ecomap1, pts)[, c(6:8)],
                            terra::extract(ecomap3, pts)[, c(2, 4, 6)])

full.ls$eco_L1 <- full.ls$NA_L1CODE
full.ls$eco_L1[is.na(full.ls$eco_L1)] <- full.ls$LEVEL1[is.na(full.ls$eco_L1)]
full.ls$eco_L2 <- full.ls$NA_L2CODE
full.ls$eco_L2[is.na(full.ls$eco_L2)] <- full.ls$LEVEL2[is.na(full.ls$eco_L2)]
full.ls$eco_L3 <- full.ls$NA_L3CODE
full.ls$eco_L3[is.na(full.ls$eco_L3)] <- full.ls$LEVEL3[is.na(full.ls$eco_L3)]

### manually fill in 
full.ls$eco_L1[which(full.ls$SITE_ID == "AR-TF1")] <- 
  full.ls$eco_L1[which(full.ls$SITE_ID == "AR-TF2")]

full.ls$eco_L1_name <- NA
for(ee in 1:nrow(ecoregion.code.l1)){
  full.ls$eco_L1_name[which(as.character(full.ls$eco_L1) == ecoregion.code.l1$code[ee])] <-
    ecoregion.code.l1$name[ee]
}
full.ls$eco_L1_name[which(full.ls$SITE_ID %in% c("US-SuS", "US-SuW", "US-SuM", "US-xPU"))] <- "Hawaii"
full.ls$eco_L1[which(full.ls$SITE_ID %in% c("US-SuS", "US-SuW", "US-SuM", "US-xPU"))] <- 25

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
  target.var <- target.var.ls[l1]
  
  src.ls <- sub("_MEDIAN.csv", "", src.list.in)
  src.ls <- src.ls[!duplicated(src.ls)]
  
  site.ls <- substr(src.ls, start = 1, stop = 6)
  res.ls <- substr(src.ls, start = 8, stop = 9)
  ver.ls <- substr(src.ls, start = 11, stop = nchar(src.ls))
  
  print(paste("################################################"))
  print(paste("######  ", target.var, "  ##################"))
  
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
        (target.var == "SWC" & site.ls[i0] == "US-BZF")) {
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

  data.pre <- NULL #array(dim = c(len.ts, num.ts, num.dim))
  
  for (i2 in 1:length(src.ls)) {
    data.tmp <- read.csv(
      paste0(path.in, src.ls[i2], "_MEDIAN.csv"),
      header = T,
      na.strings = "-9999",
      stringsAsFactors = F
    )
    data.tmp <- data.tmp[, target.var]
    data.pre <- cbind(data.pre, data.tmp)
    colnames(data.pre)[which(colnames(data.pre) == "data.tmp")] <-
      paste(site.ls[i2])
  }
  
  write.csv(
    data.pre,
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-compiled_ts.csv"),
    row.names = F
  )
  
  ###### Loop through diff distance metric
  for(dd in 1: length(dist.ls)){
    
    target.dist <- dist.ls[dd]
    target.dist.get <- ifelse(target.dist == "dtw", "dtw_basic", 
                              ifelse(target.dist == "euclidean", "Euclidean",
                                     ifelse(target.dist == "gak", "gak", NA)))
    
    #############################################################################################################
    #### Clustering Time series
    
    # In the following call to tsclust, specifying the value of k indicates the number of desired
    # clusters, so that the cutree() is called internally. Additionally, the shape extraction
    # function is provided in the centroid argument so that, once the k clusters are obtained, their
    # prototypes are extracted. Therefore, the series are z-normalized by means of the z-score
    # function. The seed is provided because of the randomness in shape extraction when choosing
    # a reference series
    
    # save distance matrics for output
    distMatrix <- proxy::dist(t(data.pre),
                              method = target.dist.get)
    
    write.csv(
      as.matrix(distMatrix),
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-distMatrix-", target.dist,".csv"),
      quote = T
    )
    
    ## main clustering function
    hc1 <- dtwclust::tsclust(
      t(data.pre),
      type = "hierarchical",
      k = n.grp,
      seed = 1234,
      distance = target.dist.get,
      centroid = "pam"
      #control = hierarchical_control(method = "average")
    )
    
    
    ##############################################################
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
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-CVI-", target.dist,".csv"),
      quote = F
    )
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-ncluster_cvi_", target.dist, ".png"),
      width = 6,
      height = 4.5,
      units = "in",
      pointsize = 10,
      res = 300
    )
    plot(
      c(n.grp),
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
      "topleft",
      lty = 1,
      col = c("black", rainbow(6)),
      ncol = 4,
      legend = c("MEAN", colnames(cvi.est)),
      lwd = 1.5,
      bty = "n"
    )
    lines(c(n.grp), apply(cvi.est, 1, mean), lwd = 3)
    abline(v = n.grp[which(apply(cvi.est, 1, mean) == max(apply(cvi.est, 1, mean)))], lty = 4)
    dev.off()
    
    ##############################################################
    ## Choose optimal tree number based on CH (Calinski-Harabasz) index
    
    n.grp.opt <- n.grp[which(apply(cvi.est, 1, mean) == max(apply(cvi.est, 1, mean)))]
    
    hc <- dtwclust::tsclust(
      t(data.pre),
      type = "hierarchical",
      k = n.grp.opt,
      seed = 1234,
      distance = target.dist.get,
      centroid = "pam"
      #control = hierarchical_control(method = "average")
    )

    ########################################################################
    ## reorder group ID to follow relative location in a tree    
    ## prepare group list, colors, centroids
    new.grp.ls <- grp.ls <- cutree(hc, k = n.grp.opt)
    
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
    #plot(new.grp.ls)
    
    ## assign group color
    if(use.prescribed.col){
      col_pan <- cluster_color_panel2(target.var = target.var)
      
    }else{
      col_pan <- data.frame(r = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 1],
                            g = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 2],
                            b = t(col2rgb(rainbow(length(unique(new.grp.ls)))))[, 3],
                            tmp = sort(unique(new.grp.ls)))
      colnames(col_pan)[4] <- target.var
    }
    
    # return cluster ID to full list
    full.ls <- data.frame(full.ls,
                          tmp = NA,
                          tmp2 = NA,
                          centroid = FALSE,
                          color = NA)
    col_pan_get <- rep(which(is.na(col_pan[, target.var])), nrow(full.ls))
    
    for (i3 in 1:length(new.grp.ls)) {
      full.ls$tmp[which(full.ls$SITE_ID == paste(names(new.grp.ls)[i3]))] <- new.grp.ls[i3]
      full.ls$centroid[which(full.ls$SITE_ID == paste(names(hc@centroids)[i3]))] <- TRUE
      
      if(new.grp.ls[i3] %in% col_pan[, target.var]){
        col_pan_get[i3] <- which(col_pan[, target.var] == new.grp.ls[i3])
        
      }else if(!is.na(new.grp.ls[i3])){
        col_pan_get[i3] <- which(col_pan[, target.var] == 99)
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
      paste(target.var, "_clust_group_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "tmp2")] <-
      paste(target.var, "_group_order_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "centroid")] <-
      paste(target.var, "_centroid_", target.dist, sep = "")
    colnames(full.ls)[which(colnames(full.ls) == "color")] <-
      paste(target.var, "_color_", target.dist, sep = "")
 
    ############################################################################
    ## Figure preparation
    ##  plot trees
    hc.plot.tree <- ape::as.phylo(hc)
    scale.facor <-
      (range(hc.plot.tree$edge.length)[2] - range(hc.plot.tree$edge.length)[1]) *
      0.01
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-tree-", target.dist,".png"),
      width = 4.5,
      height = 12,
      units = "in",
      pointsize = 8,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(0.1, 0.1, 0.1, 0.1), oma= c(0, 0, 0, 0))
    plot(
      hc.plot.tree,
      #type = "unrooted",
      tip.color = rgb(col_pan$r[col_pan_get],
                      col_pan$g[col_pan_get],
                      col_pan$b[col_pan_get], maxColorValue = 255),
      label.offset = scale.facor,#ifelse(l1 == 4, 0.5,15),
      #show.node.label=T,
      cex = 0.4
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(
        path.out,
        "AMF-diurnal-seasonal-cluster-",
        target.var,
        "-tree-fan-",
        target.dist,
        ".png"
      ),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
    plot(
      hc.plot.tree,
      type = "fan",
      tip.color = rgb(col_pan$r[col_pan_get],
                      col_pan$g[col_pan_get],
                      col_pan$b[col_pan_get], maxColorValue = 255),
      label.offset = scale.facor,#ifelse(l1 == 4, 0.5,15),
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var,
             "-tree-fan-", target.dist, "-ecoregion.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
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
             target.var,
             "-tree-fan-", target.dist, "-igbp.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
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
             target.var,
             "-tree-radial-", target.dist, ".png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = rgb(col_pan$r[col_pan_get],
                      col_pan$g[col_pan_get],
                      col_pan$b[col_pan_get], maxColorValue = 255),
      label.offset = scale.facor * ifelse(target.var %in% c("NETRAD"),
                                          0.00005, 
                                          ifelse(target.var %in% c("LE", "H", "SWC"), 
                                                 0.0001, 0.0002)),
      #show.node.label=T,
      cex = 0.65,
      no.margin = T
    )
    #text(0,0,)
    dev.off()
    
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-",
             target.var,
             "-tree-radial-", target.dist, "-ecoregion.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "ecoregion.color"],
      label.offset = scale.facor * ifelse(target.var %in% c("NETRAD"),
                                          0.00005, 
                                          ifelse(target.var %in% c("LE", "H", "SWC"), 
                                                 0.0001, 0.0002)),
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
             target.var,
             "-tree-radial-", target.dist, "-igbp.png"),
      width = 9,
      height = 9,
      units = "in",
      pointsize = 10,
      res = 300
    )
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 3.5), oma= c(0, 0, 0, 0))
    plot(
      hc.plot.tree,
      type = "radial",
      tip.color = full.ls[hc.plot.tree$tip.label,][, "igbp.color"],
      label.offset = scale.facor * ifelse(target.var %in% c("NETRAD"),
                                          0.00005, 
                                          ifelse(target.var %in% c("LE", "H", "SWC"), 
                                                 0.0001, 0.0002)),
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
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-map-", target.dist,".png"),
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
      full.ls$LOCATION_LONG[!is.na(full.ls[,paste0(target.var, "_clust_group_", target.dist)])],
      full.ls$LOCATION_LAT[!is.na(full.ls[,paste0(target.var, "_clust_group_", target.dist)])],
      col = "white", 
      bg = full.ls[,paste0(target.var, "_color_", target.dist)][!is.na(full.ls[,paste0(target.var, "_clust_group_", target.dist)])],
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
    
    ############################################################################
    ## plot grouped diurnal-seasonal plots
    png(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-", target.dist,".png"),
      width = 7,
      height = 0.5 + ceiling(n.grp.opt / 5) * 6 / 9,
      units = "in",
      pointsize = 9,
      res = 300
    )
    
    par(
      mfrow = c(ceiling(n.grp.opt / 5), 5),
      mar = c(0.1, 0.1, 0.1, 0.1),
      oma = c(3.5, 5.5, 1, 1)
    )
    
    # sort by number of sites
    #new.grp.ls.sort <- as.numeric(rev(names(sort(table(new.grp.ls)))))
    
    for (j2 in 1:n.grp.opt) {
      plot(
        0,
        0,
        xlim = c(0, len.ts),
        ylim = c(target.rng.ls[[l1]][1],
                 target.rng.ls[[l1]][length(target.rng.ls[[l1]])]),
        type = "n",
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n"
      )
      
      #data.pre.sub <- as.data.frame(data.pre[, which(new.grp.ls == new.grp.ls.sort[j2])])
      data.pre.sub <- as.data.frame(data.pre[, which(new.grp.ls == j2)])
      
      ### find prototype, either average or centroid
      #data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)
      data.pre.sub.mean <- hc@centroids[[which(names(hc@centroids) %in% names(new.grp.ls[new.grp.ls == j2]))]]
      
      for (j1 in 1:ncol(data.pre.sub)) {
        lines(data.pre.sub[, j1],
              col = "grey",
              lwd = 0.5)
      }
      
     
      if (j2 - 1 == ceiling((j2 - 1 ) / 5) * 5 ) {
        axis(2, target.rng.ls[[l1]],las = 2)
      }
      
      #if (ncol(data.pre.sub) > 1)
        if (j2 %in% col_pan[, target.var]) {
          lines(
            data.pre.sub.mean,
            col = rgb(col_pan$r[which(col_pan[, target.var] == j2)],
                      col_pan$g[which(col_pan[, target.var] == j2)],
                      col_pan$b[which(col_pan[, target.var] == j2)],
                      maxColorValue = 255),
            lwd = 0.8
          )
        }
      
      if(ncol(data.pre.sub) < 5){
        legend(
          "topleft",
          paste(names(new.grp.ls[new.grp.ls == j2])),
          col = "darkgrey",
          cex = 1,
          bty = "n"
        )
      }
      
      text(len.ts * 0.8,
           target.rng.ls[[l1]][length(target.rng.ls[[l1]])] - 0.1 * (
             target.rng.ls[[l1]][length(target.rng.ls[[l1]])] - target.rng.ls[[l1]][1]
           ),
           paste0("Group:", j2, "(", ncol(data.pre.sub) ,")"))
      
      if(j2 >= n.grp.opt - 4){
        axis(1, at = seq(0, len.ts, by = len.ts.tck), labels = FALSE)
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
      side = 2,
      target.lab.ls[[l1]],
      outer = T,
      line = 3,
      cex = 1.2
    )
    
    dev.off()
    
    
    ############################################################################
    ### subgroup plot
    sub.new.grp.ls <- as.numeric(names(table(new.grp.ls))[table(new.grp.ls) > 5])
    
    if (length(sub.new.grp.ls) > 0) {
      for (i in 1:length(sub.new.grp.ls)) {
        data.pre.sub1 <- as.data.frame(data.pre[, which(new.grp.ls == sub.new.grp.ls[i])])
        data.pre.sub.name1 <- colnames(data.pre.sub1)
        
        #data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)
        get.prototype <-
          data.pre.sub.name1[full.ls[data.pre.sub.name1, 
                                    paste(target.var, "_centroid_", target.dist, sep = "")]]
        data.pre.sub.mean <- data.pre.sub1[, get.prototype]
        
        get.order <-
          full.ls[data.pre.sub.name1, 
                  paste(target.var, "_group_order_", target.dist, sep = "")]
        get.order <- get.order - min(get.order) + 1
        
        
        data.pre.sub <- data.frame(tmp = data.pre.sub1[, which(get.order == 1)])
        colnames(data.pre.sub)[which(colnames(data.pre.sub) == "tmp")] <- data.pre.sub.name1[which(get.order == 1)]
        
        for(ii in 2:length(get.order)){
          data.pre.sub <- cbind.data.frame(data.pre.sub,
                                           data.frame(tmp = data.pre.sub1[, which(get.order == ii)]))
          colnames(data.pre.sub)[which(colnames(data.pre.sub) == "tmp")] <- data.pre.sub.name1[which(get.order == ii)]
        }
        data.pre.sub.name <- colnames(data.pre.sub)
        
        if (ncol(data.pre.sub) <= 25) {
          data.pre.sub.ls <- c(1:ncol(data.pre.sub))
        } else{
          data.pre.sub.ls <-
            ceiling(c(1:ncol(data.pre.sub)) / (ceiling(ncol(
              data.pre.sub
            ) / 25)))
          data.pre.sub.ls[data.pre.sub.ls < 1] <- 1
          data.pre.sub.ls[data.pre.sub.ls > 25] <- 25
        }
        
        png(
          paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "_subgroup_", sub.new.grp.ls[i], "-", target.dist, ".png"),
          width = 11,
          height = 9,
          units = "in",
          pointsize = 9,
          res = 300
        )
        par(
          mfrow = c(5, 5),
          mar = c(0, 0, 0, 0),
          oma = c(3, 5, 0.5, 0.5)
        )
        
        for (j in 1:min(c(max(data.pre.sub.ls), 25))) {
          plot(
            0,
            0,
            xlim = c(0, len.ts),
            ylim = range(data.pre.sub),
            type = "n",
            ylab = "",
            xaxt = "n",
            yaxt = "n",
            xlab = ""
          )
          
          #lines(data.pre.sub.mean, col = "black", lwd = 1.5)
          
          data.pre.sub.sub <-
            as.data.frame(data.pre.sub[, which(data.pre.sub.ls == j)])
          colnames(data.pre.sub.sub) <-
            data.pre.sub.name[which(data.pre.sub.ls == j)]
          
          for (j1 in 1:ncol(data.pre.sub.sub)) {
            lines(data.pre.sub.sub[, j1],
                  col = rainbow(ncol(data.pre.sub.sub))[j1])
          }
          legend(
            0,
            range(data.pre.sub)[2],
            paste(colnames(data.pre.sub.sub)),
            lty = c(1),
            col = rainbow(ncol(data.pre.sub.sub)),
            cex = 1,
            bty = "n"
          )
          if (j == 1 | j == 6 | j == 11 | j == 16 | j == 21) {
            axis(2, seq(
              round(range(data.pre.sub)[1], digits = 0),
              round(range(data.pre.sub)[2], digits = 0),
              length.out = 6
            ), las = 2)
          }
        }
        
        mtext(
          side = 1,
          "Date-Hour",
          outer = T,
          line = 1,
          cex = 1.5
        )
        mtext(
          side = 2,
          paste(target.var),
          outer = T,
          line = 3.2,
          cex = 1.5
        )
        dev.off()
      }
    }
    
    ############################################################################
    ###### explore cross table btw cluster & IGBP, Eco region
    short.ls <-
      data.frame(
        group = as.factor(full.ls[, grep(paste0(target.var, "_clust_group_", target.dist),
                                         colnames(full.ls))]),
        igbp = as.factor(full.ls$IGBP),
        ecoregion = as.factor(full.ls$eco_L1_name),
        eco_igbp = as.factor(paste0(full.ls$eco_L1_name, "_", full.ls$IGBP))
      )
    
    write.csv(
      table(short.ls$eco_igbp,
            short.ls$group),
      file = paste0(
        path.out,
        "AMF-diurnal-seasonal-cluster-",
        target.var,
        "-",
        target.dist,
        "eco_igbp.csv"
      ),
      row.names = T
    )
    
    ############################################################################
    ####### Logistic regression to explore the predictability of IGBP & Ecoregion
    sink(file = paste(
      path.out,
      "AMF-diurnal-seasonal-cluster-", target.var, "-", target.dist,
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
  paste0(path.out, "ALL_BASE_site_short_list2.csv"),
  quote = T,
  row.names = F
)
