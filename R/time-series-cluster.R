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
source(paste0(RDir, "cluster_color_panel.R"))
source(paste0(RDir, "ecoregion_table.R"))

###### Control parameter
ver <- "20230726"

len.ts <- 52 * 12  ## length of time series
n.grp <- c(15:50)  # searching window for the n of clusters/branches to keep 
        #c(10:30) for cluster ms

## original colnames from diurnal-seasonal files, revising the colname
sel.var <- c("FC", "SC", "NEE", "LE", "H",
             "USTAR", "NETRAD", "TA", "VPD", 
             "SWC")

drop.ls<-c("FC", "SC")

# target variables to run uni-variate clustering
target.var.ls <- c("NEE", "LE", "H", "USTAR",
                   "NETRAD", 
                   "TA", "VPD", "SWC"#,
                   #"FCH4", "P"
                   )
target.lab.ls <- list(expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
                      expression(LE~'('*W~m^{-2}*')'),
                      expression(H~'('*W~m^{-2}*')'),
                      expression(USTAR~'('*m~s^{-1}*')'),
                      expression(NETRAD~'('*W~m^{-2}*')'),
                      expression(TA~'('*degree~C*')'),
                      expression(VPD~'('*kPa*')'),
                      expression(SWC~'('*'%'*')')#,
                      #expression(FCH4~'('*nmole~m^{-2}~s^{-1}*')'),
                      #expression(P~'('*mm*')')
                      )
target.rng.ls <- list(seq(-60, 20, length.out = 5),
                      seq(-130, 520, length.out = 6),
                      seq(-130, 520, length.out = 6),
                      seq(0, 1.6, length.out = 5),
                      seq(-300, 900, length.out = 5),
                      seq(-40, 40, length.out = 5),
                      seq(0, 60, length.out = 5),
                      seq(0, 100, length.out = 5)#,
                      #seq(-100, 500, length.out = 6),
                      #seq(0, 1500, length.out = 6)
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
ecomap1 <- rgdal::readOGR(
  paste(map.in, "wwf_terr_ecos_oRn", sep = ""),
  layer = c("wwf_terr_ecos_oRn"),
  verbose = F
)
ecomap3 <- rgdal::readOGR(
  paste(map.in, "NA_CEC_Eco_Level3", sep = ""),
  layer = c("NA_CEC_Eco_Level3"),
  verbose = F
)

ecomap1 <- sp::spTransform(ecomap1,
                       CRS("+proj=longlat +datum=WGS84"))
ecomap3 <- sp::spTransform(ecomap3,
                       CRS("+proj=longlat +datum=WGS84"))

### work on extraction of ecoregion for each site
pts <- sp::SpatialPoints(full.ls[, c("LOCATION_LONG",
                                 "LOCATION_LAT")],
                     proj4string = CRS("+proj=longlat +datum=WGS84"))
full.ls <- cbind.data.frame(full.ls,
                            sp::over(pts, ecomap1)[, c(4, 5, 6, 7)],
                            sp::over(pts, ecomap3)[, c(1:6)])

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
                            method = "dtw_basic")
  
  ## main clustering function
  hc1 <- dtwclust::tsclust(
    t(data.pre),
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
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-ncluster_cvi.png"),
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
    t(data.pre),
    type = "hierarchical",
    k = n.grp.opt,
    seed = 1234,
    distance = "dtw_basic",
    control = hierarchical_control(method = "average")
  )
  
  grp.ls <- cutree(hc, k = n.grp.opt)
  
  write.csv(
    as.matrix(distMatrix),
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-distMatrix.csv"),
    quote = T
  )
  
  #col_pan <- cluster_color_panel(target.var = target.var)
  col_pan <- data.frame(r = t(col2rgb(rainbow(length(unique(grp.ls)))))[, 1],
                        g = t(col2rgb(rainbow(length(unique(grp.ls)))))[, 2],
                        b = t(col2rgb(rainbow(length(unique(grp.ls)))))[, 3],
                        tmp = unique(grp.ls))
  colnames(col_pan)[4] <- target.var
  
  # return cluster ID to full list
  full.ls <- data.frame(full.ls,
                        tmp = NA)
  col_pan_get <- rep(which(is.na(col_pan[, target.var])), nrow(full.ls))

  for (i3 in 1:length(grp.ls)) {
    full.ls$tmp[which(full.ls$SITE_ID == paste(names(grp.ls)[i3]))] <- grp.ls[i3]

    if(grp.ls[i3] %in% col_pan[, target.var]){
      col_pan_get[i3] <- which(col_pan[, target.var] == grp.ls[i3])

    }else if(!is.na(grp.ls[i3])){
      col_pan_get[i3] <- which(col_pan[, target.var] == 99)

    }
  }
  colnames(full.ls)[which(colnames(full.ls) == "tmp")] <-
    paste(target.var, "_clust_group", sep = "")

  ## plot trees
  png(
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-tree.png"),
    width = 4.5,
    height = 12,
    units = "in",
    pointsize = 8,
    res = 300
  )
  par(mfrow = c(1, 1), mar = c(0.1, 0.1, 0.1, 0.1), oma= c(0, 0, 0, 0))
  plot(
    ape::as.phylo(hc),
    #type = "unrooted",
    tip.color = rgb(col_pan$r[col_pan_get],
                    col_pan$g[col_pan_get],
                    col_pan$b[col_pan_get], maxColorValue = 255),
    label.offset = ifelse(l1 == 4, 0.5,15),
    #show.node.label=T,
    cex = 0.4
  )
  #text(0,0,)
  dev.off()
  
  png(
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-tree-fan.png"),
    width = 9,
    height = 9,
    units = "in",
    pointsize = 10,
    res = 300
  )
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 4.5, 4.5), oma= c(0, 0, 0, 0))
  plot(
    ape::as.phylo(hc),
    type = "fan",
    tip.color = rgb(col_pan$r[col_pan_get],
                    col_pan$g[col_pan_get],
                    col_pan$b[col_pan_get], maxColorValue = 255),
    label.offset = ifelse(l1 == 4, 0.5,15),
    #show.node.label=T,
    cex = 0.65
  )
  #text(0,0,)
  dev.off()

  ## plot grouped diurnal-seasonal plots
  png(
    paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, ".png"),
    width = 7,
    height = 0.5 + ceiling(n.grp.opt / 3) * 8.5 / 9,
    units = "in",
    pointsize = 9,
    res = 300
  )

  par(
    mfrow = c(ceiling(n.grp.opt / 3), 3),
    mar = c(0.2, 0.2, 0.2, 0.2),
    oma = c(3.5, 5.5, 1, 1)
  )

  # sort by number of sites
  grp.ls.sort <- as.numeric(rev(names(sort(table(grp.ls)))))

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

    data.pre.sub <- as.data.frame(data.pre[, which(grp.ls == grp.ls.sort[j2])])

    data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)

    for (j1 in 1:ncol(data.pre.sub)) {
      lines(data.pre.sub[, j1],
            col = "grey",
            lwd = 0.5)
    }
    if(ncol(data.pre.sub) < 5){
      legend(
        0,
        target.rng.ls[[l1]][length(target.rng.ls[[l1]])],
        paste(names(grp.ls[grp.ls == grp.ls.sort[j2]])),
        col = "darkgrey",
        cex = 1,
        bty = "n"
      )
    }
    if (j2 - 1 == ceiling((j2 - 1 ) / 3) * 3 ) {
      axis(2, target.rng.ls[[l1]],las = 2)
    }

    if (ncol(data.pre.sub) > 1)
      if (grp.ls.sort[j2] %in% col_pan[, target.var]) {
        lines(
          data.pre.sub.mean,
          col = rgb(col_pan$r[which(col_pan[, target.var] == grp.ls.sort[j2])],
                    col_pan$g[which(col_pan[, target.var] == grp.ls.sort[j2])],
                    col_pan$b[which(col_pan[, target.var] == grp.ls.sort[j2])],
                    maxColorValue = 255),
          lwd = 2
        )
      }

    if(j2 >= n.grp.opt - 2){
      axis(1, at = seq(0, 24 * 24, by = 24), labels = FALSE)
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

  ### subgroup plot
  sub.grp.ls <- as.numeric(names(table(grp.ls))[table(grp.ls) > 5])

  if (length(sub.grp.ls) > 0) {
    for (i in 1:length(sub.grp.ls)) {
      data.pre.sub <- as.data.frame(data.pre[, which(grp.ls == sub.grp.ls[i])])
      data.pre.sub.name <- colnames(data.pre.sub)
      data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)

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
        paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "_subgroup_", sub.grp.ls[i], ".png"),
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

        lines(data.pre.sub.mean, col = "black", lwd = 1.5)

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
  
  print(paste("################################################"))
  print(paste("                                                "))
  
}

write.csv(
  full.ls,
  paste0(path.out, "ALL_BASE_site_short_list2.csv"),
  quote = T,
  row.names = F
)
