rm(list = ls()) # clean working memory

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
path.in.root <- paste0(path, "diurnal-seasonality\\")
path.out.root <- paste0(path, "cluster-output\\")
sum.out.root <- paste0(path, "summary\\")

RDir <- paste0(path, "R\\")
map.in <- paste0(path, "ecoregion\\")
cru.dir <- paste0(path, "cru\\clean\\")

source(paste0(RDir, "math_utility.R"))
#source(paste0(RDir, "get_ecomap.R"))
source(paste0(RDir, "cluster_color_panel2.R"))
#source(paste0(RDir, "ecoregion_table.R"))
#source(paste0(RDir, "igbp_table.R"))

###### Control parameter
ver <- "20231019-12h7d"
ver2 <- "20241004-GPP"
dist.ls <- c("dtw")#, "euclidean")
use.prescribed.col <- T

# target variables to run uni-variate clustering
target.var.ls <- c("NEE", "LE", "H", "USTAR", "ALL_FLUX")
target.var.ls2 <- c("GPP_NT_VUT_REF",
                   "RECO_NT_VUT_REF")
target.name.ls2 <- c("GPP",
                     "RECO")
target.lab.ls <- #list(expression(GPP~'('*mu*mole~m^{-2}~s^{-1}*')'),
                #      expression(RECO~'('*mu*mole~m^{-2}~s^{-1}*')'))
                 list(expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
                       expression(LE~'('*W~m^{-2}*')'),
                       expression(H~'('*W~m^{-2}*')'),
                       expression(USTAR~'('*m~s^{-1}*')'),
                      NULL)
target.lab.ls2 <- list(expression(GPP~'('*mu*mole~m^{-2}~s^{-1}*')'),
                       expression(RECO~'('*mu*mole~m^{-2}~s^{-1}*')'))

path.in <- paste0(path.in.root, ver, "\\")

if (!dir.exists(paste(path.out.root, ver, sep = "")))
  dir.create(paste(path.out.root, ver, sep = ""))
path.out <- paste(path.out.root, ver, "\\", sep = "")
path.out2 <- paste(path.out.root, ver2, "\\", sep = "")

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

full.ls2 <-
  read.csv(
    paste(path.out2, "ALL_FLUXNET_site_list2.csv", sep = ""),
    na = "-9999",
    header = T,
    stringsAsFactors = F
  )
full.ls2 <- full.ls2[!is.na(full.ls2$SITE_ID), ]
rownames(full.ls2) <- paste(full.ls2$SITE_ID)

## prepare for CRU mapping
full.ls$LAT <-
  ceiling((full.ls$LOCATION_LAT + 90) / 0.5)  ## location index for CRU matrix
full.ls$LONG <-
  ceiling((full.ls$LOCATION_LONG + 180) / 0.5)  ## location index for CRU matrix
full.ls2$LAT <-
  ceiling((full.ls2$LOCATION_LAT + 90) / 0.5)  ## location index for CRU matrix
full.ls2$LONG <-
  ceiling((full.ls2$LOCATION_LONG + 180) / 0.5)  ## location index for CRU matrix

tmp.yr.mn2 <-
  read.csv(paste0(cru.dir, "cru.4.05_world_year_tmp_1981_2020.csv"),
           header = T)
pre.yr.mn2 <-
  read.csv(paste0(cru.dir, "cru.4.05_world_year_pre_1981_2020.csv"),
           header = T)
pet.yr.mn2 <-
  read.csv(paste0(cru.dir, "cru.4.05_world_year_pet_1981_2020.csv"),
           header = T)

full.ls$tmp2 <- NA # 1981-2013
full.ls$pre2 <- NA
full.ls$pet2 <- NA
full.ls2$tmp2 <- NA # 1981-2013
full.ls2$pre2 <- NA
full.ls2$pet2 <- NA

for (k in 1:nrow(full.ls)) {
  full.ls$tmp2[k] <-
    tmp.yr.mn2[full.ls$LAT[k], full.ls$LONG[k]]
  full.ls$pre2[k] <-
    pre.yr.mn2[full.ls$LAT[k], full.ls$LONG[k]]
  full.ls$pet2[k] <-
    pet.yr.mn2[full.ls$LAT[k], full.ls$LONG[k]]
}

for (k in 1:nrow(full.ls2)) {
  full.ls2$tmp2[k] <-
    tmp.yr.mn2[full.ls2$LAT[k], full.ls2$LONG[k]]
  full.ls2$pre2[k] <-
    pre.yr.mn2[full.ls2$LAT[k], full.ls2$LONG[k]]
  full.ls2$pet2[k] <-
    pet.yr.mn2[full.ls2$LAT[k], full.ls2$LONG[k]]
}

################################################################################
##
## Work for all flux
##
################################################################################
for (l1 in 1:length(target.var.ls)) {
  
  target.var <- target.var.ls[l1]
  
  ## get annual sum
  full.ls$mean <- NA
  if(target.var != "ALL_FLUX"){
    data.pre <- read.csv(
      paste0(path.out, "AMF-diurnal-seasonal-cluster-", target.var, "-compiled_ts.csv"),
      header = T)
    colnames(data.pre) <- gsub("\\.", "\\-", colnames(data.pre))
    
    for(l11 in 1:ncol(data.pre)){
      full.ls$mean[which(full.ls$SITE_ID == colnames(data.pre)[l11])] <- mean(data.pre[, l11], na.rm = T)
    }
  }
  
  sub.ls <- full.ls[, c("SITE_ID",
                        paste0(target.var, "_clust_group_dtw"),
                        paste0(target.var, "_color_dtw"),
                        "tmp2",
                        "pre2",
                        "pet2",
                        "mean")]
  #sub.ls <- na.omit(sub.ls)
  colnames(full.ls)[which(colnames(full.ls) == "mean")] <- 
    paste0(target.var, "_mean")
  
  sub.ls[, 2] <- as.numeric(paste(sub.ls[, 2]))
  
  ####
  cluster_color_panel_get <- cluster_color_panel2(target.var)
  cluster_color_panel_get[which(cluster_color_panel_get[, 4] == 99), 4] <- "others"
  cluster_color_panel_get <- cluster_color_panel_get[-which(is.na(cluster_color_panel_get[, 4])), ]
    
  ## Annual sum distribution
  if (target.var != "ALL_FLUX") {
    png(
      paste0(path.out, "Annual_sum_", target.var , "_univariate.png", sep = ""),
      height = 4,
      width = 3.5,
      res = 300,
      units = "in",
      pointsize = 10
    )
    par(
      oma = c(4, 4, 0.2, 0.2),
      mar = c(0, 0, 0, 0),
      mfrow = c(1, 1)
    )
    xbox <-
      boxplot(
        sub.ls$mean[sub.ls[, 2] %in% cluster_color_panel_get[, 4]] ~
          sub.ls[, 2][sub.ls[, 2] %in% cluster_color_panel_get[, 4]],
        col = rgb(
          cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
          maxColorValue = 255
        ),
        ylim = quantile(sub.ls$mean, c(0.01, 0.99), na.rm = T),
        xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
        frame.plot = T,
        xlab = "",
        ylab = "",
        outline = F,
        border = rgb(
          cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
          maxColorValue = 255
        ),
        lty = 1,
        las = 1,
        #xaxt = "n",
        #yaxt = "n"
      )
    
    points(
      rep(nrow(cluster_color_panel_get),
          length(sub.ls$mean[!sub.ls[, 2] %in% cluster_color_panel_get[, 4]])),
      sub.ls$mean[!sub.ls[, 2] %in% cluster_color_panel_get[, 4]],
      col = rgb(
        cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
        cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
        cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
        maxColorValue = 255
      ),
      cex = 0.8,
      pch = 20
    )
    
    points(
      xbox$group,
      xbox$out,
      pch = 4,
      cex = 0.7,
      col = "black"
    )
    axis(1, 
         at = nrow(cluster_color_panel_get),
         labels = "others",
         cex.axis = 0.9)
    mtext(
      target.lab.ls[[l1]],
      side = 2,
      line = 2,
      font = 1,
      outer = T,
      cex = 1
    )
    mtext(
      "Cluster",
      side = 1,
      line = 2,
      font = 1,
      outer = T,
      cex = 1
    )
    dev.off()
    
  }
  
  ## Climatic map + IGBP
  png(
    paste0(path.out, "Climatic_space_", target.var , "_cluster.png", sep = ""),
    height = 4,
    width = 4,
    res = 300,
    units = "in",
    pointsize = 10
  )
  par(
    oma = c(3, 4.5, 0.5, 0.5),
    mar = c(0, 0, 0, 0),
    mfrow = c(1, 1)
  )
  plot(
    sub.ls$tmp2,
    sub.ls$pre2,
    cex = 1,
    pch = 16,
    col = sub.ls[, paste0(target.var, "_color_dtw")],
    ylim = c(0, 4000),
    xlim = c(-15, 31),
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )
  text(
    -16,
    4100,
    paste0("(", letters[l1], ") ", target.var),
    adj = c(0, 1),
    font = 1,
    cex = 0.9
  )
  axis(1, at = seq(-20, 30, by = 10), cex.axis = 0.9,
       mgp = c(1, 0.7, 0))
  mtext(
    expression(Precipitation ~ (mm ~ year ^ {
      -1
    })),
    side = 2,
    line = 3,
    font = 1,
    outer = T,
    cex = 1
  )
  legend(
    -17,
    3900,
    legend = paste("Cluster", cluster_color_panel_get[, 4]),
    fill = rgb(cluster_color_panel_get$r,
               cluster_color_panel_get$g,
               cluster_color_panel_get$b,
               maxColorValue = 255),
    border = NA,
    bty = "n",
    cex = 0.85,
    ncol = 3
  )
  axis(
    2,
    at = seq(0, 4000, by = 1000),
    las = 1,
    cex.axis = 0.9,
    mgp = c(3, 0.8, 0)
  )
  axis(1, at = seq(-10, 30, by = 10), cex.axis = 0.9,
       mgp = c(1, 0.7, 0))
  mtext(
    "Mean Air Temperature (°C)",
    side = 1,
    line = 2,
    font = 1,
    outer = T,
    cex = 1
  )
  
  dev.off()  
  
  print(paste("################################################"))
  print(paste("                                                "))
}

################################################################################
## Work for GPP & RECO
##
################################################################################

for (l1 in 1:length(target.var.ls2)) {
  
  target.var <- target.var.ls2[l1]
  
  ## get annual sum
  full.ls2$mean <- NA
  if(target.var != "ALL_FLUX"){
    data.pre <- read.csv(
      paste0(path.out2, "AMF-diurnal-seasonal-cluster-", target.var, "-compiled_ts.csv"),
      header = T)
    colnames(data.pre) <- gsub("\\.", "\\-", colnames(data.pre))
    
    for(l11 in 1:ncol(data.pre)){
      full.ls2$mean[which(full.ls2$SITE_ID == colnames(data.pre)[l11])] <- mean(data.pre[, l11], na.rm = T)
    }
  }
  
  sub.ls2 <- full.ls2[, c("SITE_ID",
                        paste0(target.var, "_clust_group_dtw"),
                        paste0(target.var, "_color_dtw"),
                        "tmp2",
                        "pre2",
                        "pet2",
                        "mean")]
  #sub.ls <- na.omit(sub.ls)
  colnames(full.ls2)[which(colnames(full.ls2) == "mean")] <- 
    paste0(target.var, "_mean")
  
  sub.ls2[, 2] <- as.numeric(paste(sub.ls2[, 2]))
  
  ####
  cluster_color_panel_get <- cluster_color_panel2(target.var)
  cluster_color_panel_get[which(cluster_color_panel_get[, 4] == 99), 4] <- "others"
  cluster_color_panel_get <- cluster_color_panel_get[-which(is.na(cluster_color_panel_get[, 4])), ]
  
  ## Annual sum distribution
  if (target.var != "ALL_FLUX") {
    png(
      paste0(path.out, "Annual_sum_", target.var , "_univariate.png", sep = ""),
      height = 4,
      width = 3.5,
      res = 300,
      units = "in",
      pointsize = 10
    )
    par(
      oma = c(4, 4, 0.2, 0.2),
      mar = c(0, 0, 0, 0),
      mfrow = c(1, 1)
    )
    xbox <-
      boxplot(
        sub.ls2$mean[sub.ls2[, 2] %in% cluster_color_panel_get[, 4]] ~
          sub.ls2[, 2][sub.ls2[, 2] %in% cluster_color_panel_get[, 4]],
        col = rgb(
          cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
          maxColorValue = 255
        ),
        ylim = quantile(sub.ls2$mean, c(0.01, 0.99), na.rm = T),
        xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
        frame.plot = T,
        xlab = "",
        ylab = "",
        outline = F,
        border = rgb(
          cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
          cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
          maxColorValue = 255
        ),
        lty = 1,
        las = 1
        #xaxt = "n",
        #yaxt = "n"
      )
    
    points(
      rep(nrow(cluster_color_panel_get),
          length(sub.ls2$mean[!sub.ls2[, 2] %in% cluster_color_panel_get[, 4]])),
      sub.ls2$mean[!sub.ls2[, 2] %in% cluster_color_panel_get[, 4]],
      col = rgb(
        cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
        cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
        cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
        maxColorValue = 255
      ),
      cex = 0.8,
      pch = 20
    )
    
    points(
      xbox$group,
      xbox$out,
      pch = 4,
      cex = 0.7,
      col = "black"
    )
    axis(1, 
         at = nrow(cluster_color_panel_get),
         labels = "others",
         cex.axis = 0.9)
    mtext(
      target.lab.ls2[[l1]],
      side = 2,
      line = 2,
      font = 1,
      outer = T,
      cex = 1
    )
    mtext(
      "Cluster",
      side = 1,
      line = 2,
      font = 1,
      outer = T,
      cex = 1
    )
    dev.off()
    
  }
  
  ## Climatic map + IGBP
  png(
    paste0(path.out, "Climatic_space_", target.var , "_cluster.png", sep = ""),
    height = 4,
    width = 4,
    res = 300,
    units = "in",
    pointsize = 10
  )
  par(
    oma = c(3, 4.5, 0.5, 0.5),
    mar = c(0, 0, 0, 0),
    mfrow = c(1, 1)
  )
  plot(
    sub.ls2$tmp2,
    sub.ls2$pre2,
    cex = 1,
    pch = 16,
    col = sub.ls2[, paste0(target.var, "_color_dtw")],
    ylim = c(0, 4000),
    xlim = c(-15, 31),
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )
  text(
    -16,
    4100,
    paste0("(", letters[l1 + 4], ") ", target.name.ls2[l1]),
    adj = c(0, 1),
    font = 1,
    cex = 0.9
  )
  axis(1, at = seq(-20, 30, by = 10), cex.axis = 0.9,
       mgp = c(1, 0.7, 0))
  mtext(
    expression(Precipitation ~ (mm ~ year ^ {
      -1
    })),
    side = 2,
    line = 3,
    font = 1,
    outer = T,
    cex = 1
  )
  legend(
    -17,
    3900,
    legend = paste("Cluster", cluster_color_panel_get[, 4]),
    fill = rgb(cluster_color_panel_get$r,
               cluster_color_panel_get$g,
               cluster_color_panel_get$b,
               maxColorValue = 255),
    border = NA,
    bty = "n",
    cex = 0.85,
    ncol = 3
  )
  axis(
    2,
    at = seq(0, 4000, by = 1000),
    las = 1,
    cex.axis = 0.9,
    mgp = c(3, 0.8, 0)
  )
  axis(1, at = seq(-10, 30, by = 10), cex.axis = 0.9,
       mgp = c(1, 0.7, 0))
  mtext(
    "Mean Air Temperature (°C)",
    side = 1,
    line = 2,
    font = 1,
    outer = T,
    cex = 1
  )
  
  dev.off()  
  
  print(paste("################################################"))
  print(paste("                                                "))
}

cluster_color_panel_get <- cluster_color_panel2("ALL_FLUX")
cluster_color_panel_get[which(cluster_color_panel_get[, 4] == 99), 4] <- "others"
cluster_color_panel_get <- cluster_color_panel_get[-which(is.na(cluster_color_panel_get[, 4])), ]

################################################################################
## Figure 4 in the main text 
## Annual sum vs cluster mutivariate
##
################################################################################
jpeg(
  paste0(path.out, "Annual_sum_All_FLUX_cluster.jpg", sep = ""),
  height = 8,
  width = 8,
  res = 300,
  units = "in",
  pointsize = 10
)
par(
  oma = c(3, 4, 0.3, 0.3),
  mar = c(0, 0, 0.5, 0.5),
  fig = c(0, 0.3, 0.6, 0.9)
)
plot(
  full.ls$NEE_mean,
  full.ls$LE_mean,
  xlim = c(-4, 2),
  ylim= c(0, 120),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(-4, 120, "(b)", adj = c(0, 1))
axis(2, 
     at = seq(0, 120, by = 20), 
     cex.axis = 0.9,
     mgp = c(3, 0.8, 0),
     las = 1,)
mtext(
  expression(LE ~ (W ~ m ^ { -2 })),
  side = 2,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

par(
  fig = c(0, 0.3, 0.9, 1.0),
  new = T
)
xbox <-
  boxplot(
    full.ls$NEE_mean[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX] ~ 
      full.ls$ALL_FLUX_clust_group_dtw[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX],
    col = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    ylim = c(-4, 2),
    xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
    frame.plot = F,
    xlab = "",
    ylab = "",
    horizontal = T,
    outline = F,
    border = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    lty = 1,
    las = 1,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )
points(
  full.ls$NEE_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]],
  rep(nrow(cluster_color_panel_get),
      length(full.ls$NEE_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]])),
  col = rgb(
    cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
    maxColorValue = 255
  ),
  cex = 0.7,
  pch = 20
)
points(
  xbox$out,
  xbox$group,
  pch = 4,
  cex = 0.5,
  col = "black"
)
text(-4, nrow(cluster_color_panel_get) + 0.5, "(a)", adj = c(0, 1))

par(
  fig = c(0.3, 0.4, 0.6, 0.9),
  new = T
)
xbox <-
  boxplot(
    full.ls$LE_mean[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX] ~ 
      full.ls$ALL_FLUX_clust_group_dtw[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX],
    col = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    ylim = c(0, 120),
    xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
    frame.plot = F,
    xlab = "",
    ylab = "",
    horizontal = F,
    outline = F,
    border = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    lty = 1,
    las = 1,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )

points(
  rep(nrow(cluster_color_panel_get),
      length(full.ls$LE_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]])),
  full.ls$LE_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]],
  col = rgb(
    cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
    maxColorValue = 255
  ),
  cex = 0.7,
  pch = 20
)
points(
  xbox$group,
  xbox$out,
  pch = 4,
  cex = 0.5,
  col = "black"
)
text(0.5, 120, "(c)", adj = c(0, 1))

par(
  fig = c(0, 0.3, 0.3, 0.6),
  new = T
)
plot(
  full.ls$NEE_mean,
  full.ls$H_mean,
  xlim = c(-4, 2),
  ylim = c(0, 100),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(-4, 100, "(d)", adj = c(0, 1))
axis(2, 
     at = seq(0, 100, by = 20), 
     cex.axis = 0.9,
     mgp = c(3, 0.8, 0),
     las = 1,)
mtext(
  expression(H ~ (W ~ m ^ { -2 })),
  side = 2,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

par(
  fig = c(0.3, 0.6, 0.3, 0.6),
  new = T
)
plot(
  full.ls$LE_mean,
  full.ls$H_mean,
  ylim = c(0, 100),
  xlim = c(0, 120),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(0, 100, "(e)", adj = c(0, 1))

par(
  fig = c(0.6, 0.7, 0.3, 0.6),
  new = T
)
xbox <-
  boxplot(
    full.ls$H_mean[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX] ~ 
      full.ls$ALL_FLUX_clust_group_dtw[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX],
    col = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    ylim = c(0, 100),
    xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
    frame.plot = F,
    xlab = "",
    ylab = "",
    horizontal = F,
    outline = F,
    border = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    lty = 1,
    las = 1,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )
points(
  rep(nrow(cluster_color_panel_get),
      length(full.ls$H_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]])),
  full.ls$H_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]],
  col = rgb(
    cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
    maxColorValue = 255
  ),
  cex = 0.7,
  pch = 20
)
points(
  xbox$group,
  xbox$out,
  pch = 4,
  cex = 0.5,
  col = "black"
)
text(0.5, 100, "(f)", adj = c(0, 1))

### 3rd row
par(
  fig = c(0, 0.3, 0, 0.3),
  new = T
)
plot(
  full.ls$NEE_mean,
  full.ls$USTAR_mean,
  xlim = c(-4, 2),
  ylim = c(0.1, 0.75),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(-4, 0.75, "(g)", adj = c(0, 1))
axis(2, 
     at = seq(0.1, 0.7, by = 0.1), 
     cex.axis = 0.9,
     mgp = c(3, 0.8, 0),
     las = 1,)
mtext(
  expression(USTAR ~ (m ~ s ^ { -1 })),
  side = 2,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)
axis(1, 
     at = seq(-4, 3, by = 1), 
     cex.axis = 0.9,
     mgp = c(1, 0.7, 0))
mtext(
  expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
  side = 1,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

par(
  fig = c(0.3, 0.6, 0, 0.3),
  new = T
)
plot(
  full.ls$LE_mean,
  full.ls$USTAR_mean,
  xlim = c(0, 120),
  ylim = c(0.1, 0.75),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(0, 0.75, "(h)", adj = c(0, 1))
axis(1, 
     at = seq(0, 100, by = 20), 
     cex.axis = 0.9,
     mgp = c(1, 0.7, 0))
mtext(
  expression(LE ~ (W ~ m ^ { -2 })),
  side = 1,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

par(
  fig = c(0.6, 0.9, 0, 0.3),
  new = T
)
plot(
  full.ls$H_mean,
  full.ls$USTAR_mean,
  ylim = c(0.1, 0.75),
  xlim = c(0, 100),
  cex = 0.8,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xaxs = "i",
  yaxs = "i"
)
text(0, 0.75, "(i)", adj = c(0, 1))
axis(1, 
     at = seq(0, 100, by = 20), 
     cex.axis = 0.9,
     mgp = c(1, 0.7, 0))
mtext(
  expression(H ~ (W ~ m ^ { -2 })),
  side = 1,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

par(
  fig = c(0.9, 1.0, 0.0, 0.3),
  new = T
)
xbox <-
  boxplot(
    full.ls$USTAR_mean[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX] ~ 
      full.ls$ALL_FLUX_clust_group_dtw[full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get$ALL_FLUX],
    col = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    ylim = c(0.1, 0.75),
    xlim = c(0.5, nrow(cluster_color_panel_get) + 0.5),
    frame.plot = F,
    xlab = "",
    ylab = "",
    horizontal = F,
    outline = F,
    border = rgb(
      cluster_color_panel_get$r[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$g[-nrow(cluster_color_panel_get)],
      cluster_color_panel_get$b[-nrow(cluster_color_panel_get)],
      maxColorValue = 255
    ),
    lty = 1,
    las = 1,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )
points(
  rep(nrow(cluster_color_panel_get),
      length(full.ls$USTAR_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]])),
  full.ls$USTAR_mean[!full.ls$ALL_FLUX_clust_group_dtw %in% cluster_color_panel_get[, 4]],
  col = rgb(
    cluster_color_panel_get$r[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$g[nrow(cluster_color_panel_get)],
    cluster_color_panel_get$b[nrow(cluster_color_panel_get)],
    maxColorValue = 255
  ),
  cex = 0.7,
  pch = 20
)
points(
  xbox$group,
  xbox$out,
  pch = 4,
  cex = 0.5,
  col = "black"
)
text(0.5, 0.75, "(j)", adj = c(0, 1))

par(
  fig = c(0.4, 0.6, 0.6, 0.9),
  new = T
)
plot(NULL,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     bty = "n")
legend(
  "topleft",
  legend = paste("Cluster", cluster_color_panel_get[, 4]),
  fill = rgb(cluster_color_panel_get$r,
             cluster_color_panel_get$g,
             cluster_color_panel_get$b,
             maxColorValue = 255),
  border = NA,
  bty = "n",
  cex = 1,
  ncol = 1
)

par(
  fig = c(0.62, 0.92, 0.62, 0.92),
  new = T
)
plot(
  full.ls$tmp2,
  full.ls$pre2,
  pch = 16,
  col = full.ls$ALL_FLUX_color_dtw,
  ylim = c(0, 4000),
  xlim = c(-15, 31),
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  xec = 0.8 
)
text(
  -16,
  4100,
  "(k)",
  adj = c(0, 1)
)
axis(3, 
     at = seq(-20, 30, by = 10), 
     cex.axis = 0.9,
     mgp = c(1, 0.7, 0))
mtext(
  expression(Precipitation ~ (mm ~ year ^ {
    -1
  })),
  side = 4,
  line = 3,
  font = 1,
  outer = F,
  cex = 1
)
axis(
  4,
  at = seq(0, 4000, by = 1000),
  las = 1,
  cex.axis = 0.9,
  mgp = c(3, 0.8, 0)
)
mtext(
  "Mean Air Temperature (°C)",
  side = 3,
  line = 2,
  font = 1,
  outer = F,
  cex = 1
)

dev.off()  

################################################################################
### Prepare supplement Table S1
################################################################################
out.ls <- merge.data.frame(full.ls,
                           full.ls2[, c("SITE_ID",
                                        "RECO_NT_VUT_REF_availability",
                                        "RECO_NT_VUT_REF_availability_after_filling",
                                        "GPP_NT_VUT_REF_availability",
                                        "GPP_NT_VUT_REF_availability_after_filling",
                                        "GPP_NT_VUT_REF_clust_group_dtw",
                                        "RECO_NT_VUT_REF_clust_group_dtw",
                                        "GPP_NT_VUT_REF_mean",
                                        "RECO_NT_VUT_REF_mean")],
                           by = "SITE_ID",
                           all.x = T)

out.ls <- out.ls[, c("SITE_ID", 
                     "data_year_length",
                     "IGBP", "eco_L1_name",
                     "LOCATION_LAT", "LOCATION_LONG", "LOCATION_ELEV",
                     "tmp2", "pre2", 
                     "NEE_availability", "NEE_availability_after_filling",   
                     "LE_availability", "LE_availability_after_filling",    
                     "H_availability", "H_availability_after_filling",     
                     "USTAR_availability", "USTAR_availability_after_filling", 
                     "NETRAD_availability", "NETRAD_availability_after_filling",
                     "TA_availability", "TA_availability_after_filling",    
                     "VPD_availability", "VPD_availability_after_filling",   
                     "SWC_availability", "SWC_availability_after_filling",
                     "RECO_NT_VUT_REF_availability",
                     "RECO_NT_VUT_REF_availability_after_filling",
                     "GPP_NT_VUT_REF_availability",
                     "GPP_NT_VUT_REF_availability_after_filling",
                     "NEE_clust_group_dtw",
                     "LE_clust_group_dtw",
                     "H_clust_group_dtw",
                     "USTAR_clust_group_dtw",
                     "NETRAD_clust_group_dtw",
                     "TA_clust_group_dtw",
                     "VPD_clust_group_dtw",
                     "SWC_clust_group_dtw",
                     "ALL_FLUX_clust_group_dtw",
                     "GPP_NT_VUT_REF_clust_group_dtw",
                     "RECO_NT_VUT_REF_clust_group_dtw",
                     "NEE_mean",
                     "LE_mean",
                     "H_mean",
                     "USTAR_mean",
                     "GPP_NT_VUT_REF_mean",
                     "RECO_NT_VUT_REF_mean")]
colnames(out.ls) <- c("SITE_ID", 
                      "Data_year_length",
                      "IGBP", "ECOREGION",
                      "LOCATION_LAT", "LOCATION_LONG", "LOCATION_ELEV",
                      "MAT_CRU", "MAP_CRU", 
                      "NEE_availability", "NEE_availability_after_filling",   
                      "LE_availability", "LE_availability_after_filling",    
                      "H_availability", "H_availability_after_filling",     
                      "USTAR_availability", "USTAR_availability_after_filling", 
                      "NETRAD_availability", "NETRAD_availability_after_filling",
                      "TA_availability", "TA_availability_after_filling",    
                      "VPD_availability", "VPD_availability_after_filling",   
                      "SWC_availability", "SWC_availability_after_filling",
                      "RECO_availability",
                      "RECO_availability_after_filling",
                      "GPP_availability",
                      "GPP_availability_after_filling",
                      "NEE_cluster_group",
                      "LE_cluster_group",
                      "H_cluster_group",
                      "USTAR_cluster_group",
                      "NETRAD_cluster_group",
                      "TA_cluster_group",
                      "VPD_cluster_group",
                      "SWC_cluster_group",
                      "ALL_FLUX_cluster_group",
                      "GPP_cluster_group",
                      "RECO_cluster_group",
                      "NEE_mean",
                      "LE_mean",
                      "H_mean",
                      "USTAR_mean",
                      "GPP_mean",
                      "RECO_mean")  

write.csv(out.ls,
          paste0(sum.out.root, ver2, "\\",
                 "TableS1_site_list_harmonic_202411rev.csv"),
          row.names = F,
          quote = F)