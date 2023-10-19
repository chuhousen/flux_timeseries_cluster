rm(list = ls()) # clean working memory

# require("pdc")
# require("proxy")
# require("dtw")
require("rgdal")
require('amerifluxr')
# require("ape") ##  for displaying more appealing trees
# require("dtwclust")


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
ver <- "20230117"

len.ts <- 24 * 12
#len.ts <- 24 * 24  ## length of time series

# target variables to run uni-variate clustering
target.var.ls <- "LE"
target.lab.ls <- list(expression(ET~'('*W~m^{-2}*')'))
target.rng.ls <- list(seq(0, 400, length.out = 5))

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

sgi <- amerifluxr::amf_site_info()
sgi <- sgi[sgi$DATA_POLICY == "CCBY4.0",]

## original colnames from diurnal-seasonal files, revising the colname
sel.var <- c("LE", "TA")
drop.ls<-c("FC","SC","WS")

## keep only CCBY4.0
full.ls <- full.ls[full.ls$SITE_ID %in% c(sgi$SITE_ID, "BR-Sa1"), ]

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

### Drop problematic cases
## drop.ls<-c("BR-CST","CA-SCB","US-Elm","US-Esm","US-Hn2","US-NC2",)


### read in Ecoregion map
ecomap1 <- rgdal::readOGR(
  paste(map.in, "sa_eco_l3", sep = ""),
  layer = c("sa_eco_l3"),
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
                            sp::over(pts, ecomap1)[, c(5:7)],
                            sp::over(pts, ecomap3)[, c(1:6)])
full.ls$eco_L1 <- full.ls$NA_L1CODE
full.ls$eco_L1[is.na(full.ls$eco_L1)] <- full.ls$LEVEL1[is.na(full.ls$eco_L1)]
full.ls$eco_L2 <- full.ls$NA_L2CODE
full.ls$eco_L2[is.na(full.ls$eco_L2)] <- full.ls$LEVEL2[is.na(full.ls$eco_L2)]
full.ls$eco_L3 <- full.ls$NA_L3CODE
full.ls$eco_L3[is.na(full.ls$eco_L3)] <- full.ls$LEVEL3[is.na(full.ls$eco_L3)]

full.ls$eco_L1[full.ls$SITE_ID == "AR-TF1"] <- full.ls$eco_L1[full.ls$SITE_ID == "AR-TF2"]
full.ls <- full.ls[!is.na(full.ls$eco_L1), ]

full.ls$ecoregion.color <- NA
full.ls$ecoregion.color2 <- NA
for (i in 1:nrow(full.ls)) {
  if (!is.na(full.ls$eco_L1[i])) {
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
  }
}

for(ee in 1:nrow(ecoregion.color.table)){
  ecoregion.color.table$n.site[ee] <-
    nrow(full.ls[full.ls$eco_L1 == as.character(ecoregion.color.table$ecoregion[ee]),])
}

src.list.in <- src.list.in[which(substr(src.list.in,
                                        start = 1,
                                        stop = 6) %in% full.ls$SITE_ID)]

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
    ta.tmp <- data.tmp[, "TA"]
    ta.tmp[is.na(ta.tmp)] <- mean(ta.tmp, na.rm = T)
    data.tmp <- data.tmp[, target.var]
    
    ## unit conversion to mm hr-1
    #data.tmp <- data.tmp / 1000000 / (2.501 - (0.002361) *  ta.tmp) * 3600
    
    data.pre <- cbind(data.pre, data.tmp)
    colnames(data.pre)[which(colnames(data.pre) == "data.tmp")] <-
      paste(site.ls[i2])
  }
 
  # sort by number of sites
  grp.ls <- table(full.ls$eco_L1)
  grp.ls1 <- grp.ls[grp.ls >= 8]
  grp.ls2 <- grp.ls[grp.ls < 8]
  
  #grp.ls.sort <- as.numeric((names(sort(grp.ls))))
  grp.ls1.sort <- c(8, 11, 9, 5, 12, 4, 6, 7, 3, 10, 2)
  grp.ls2.sort <- c(21, 20, 14, 16, 13, 15, 19)
  
  ## plot grouped diurnal-seasonal plots
  png(
    paste0(path.out, Sys.Date(), "_ET-diurnal-seasonal-ecoregion2.png"),
    width = 9,
    height = 6,
    units = "in",
    pointsize = 9,
    res = 300
  )
  par(
    mar = c(0, 2.5, 0.5, 0),
    oma = c(4, 2, 4, 0.5),
    mfrow = c(2, 1)
  )
  plot(
    0,
    0,
    xaxs = "i",
    xlim = c(0, len.ts),
    ylim = c(target.rng.ls[[l1]][1],
             420),
    type = "n",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )
  for (j2 in 1:length(grp.ls1.sort)) {
    grp.ls.site <- full.ls$SITE_ID[full.ls$eco_L1 == grp.ls1.sort[j2]]
    data.pre.sub <- as.data.frame(data.pre[, grp.ls.site])
    data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)
    
    points(
      data.pre.sub.mean,
      pch = 16,
      cex = 0.7,
      col = full.ls$ecoregion.color2[which(full.ls$SITE_ID == colnames(data.pre.sub)[1])],
      lwd = 1.5,
      lty = 1
    )
    lines(
      data.pre.sub.mean,
      col = full.ls$ecoregion.color2[which(full.ls$SITE_ID == colnames(data.pre.sub)[1])],
      lwd = 1.5,
      lty = 1
    )
    print(max(data.pre.sub.mean))
  }
  axis(2, target.rng.ls[[l1]],las = 2)
  axis(
    side = 3,
    at = seq(12, 24*12-12, by = 24),
    labels = seq(1, 12, by = 1),
    cex.axis = 0.9)
  # axis(
  #   side = 3,
  #   at = seq(12, 576-12, by = 24),
  #   labels = seq(8, 353, by = 15),
  #   cex.axis = 0.75)
  abline(v = seq(24, len.ts-24, by = 24), lty = 3, col = "grey")
  text(3, 420, "(a)", adj = c(0, 1))
  legend(
    6,
    420,
    legend = paste0(ecoregion.color.table$plot.text[ecoregion.color.table$ecoregion %in% grp.ls1.sort],
                    " (", ecoregion.color.table$n.site[ecoregion.color.table$ecoregion %in% grp.ls1.sort], ")"),
    fill = ecoregion.color.table$col2[ecoregion.color.table$ecoregion %in% grp.ls1.sort],
    border = NA,
    bty = "n",
    cex = 1,
    ncol = 2
  )
  mtext(side = 3,
        "Month",
        line = 2.5,
        outer = F,
        cex = 1.2)
  # mtext(side = 3,
  #       "DOY (central date)",
  #       line = 2.5,
  #       outer = F,
  #       cex = 1.2)
  mtext(
    side = 2,
    target.lab.ls[[l1]],
    outer = T,
    line = 0.4,
    cex = 1.2
  )
  
  plot(
    0,
    0,
    xaxs = "i",
    xlim = c(0, len.ts),
    ylim = c(target.rng.ls[[l1]][1],
             420),
    type = "n",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )
  for (j3 in 1:length(grp.ls2.sort)) {
    grp.ls.site <- full.ls$SITE_ID[full.ls$eco_L1 == grp.ls2.sort[j3]]
    data.pre.sub <- as.data.frame(data.pre[, grp.ls.site])
    if(ncol(data.pre.sub) == 1){
      data.pre.sub.mean <- data.pre.sub
      colnames(data.pre.sub)[1] <- grp.ls.site
    }else{
      data.pre.sub.mean <- apply(data.pre.sub, 1, na.mean)  
    }

    points(
      data.pre.sub.mean,
      pch = 16,
      cex = 0.7,
      col = full.ls$ecoregion.color2[which(full.ls$SITE_ID == colnames(data.pre.sub)[1])],
      lwd = 1.5,
      lty = 1
    )
    lines(
      data.pre.sub.mean,
      col = full.ls$ecoregion.color2[which(full.ls$SITE_ID == colnames(data.pre.sub)[1])],
      lwd = 1.5,
      lty = 1
    )
    print(max(data.pre.sub.mean))
  }
  axis(2, target.rng.ls[[l1]],las = 2)
  abline(v = seq(24, len.ts-24, by = 24), lty = 3, col = "grey")
  text(3, 420, "(b)", adj = c(0, 1))
  legend(
    len.ts * 0.3,
    420 * 1.05,
    legend = paste0(ecoregion.color.table$plot.text[ecoregion.color.table$ecoregion %in% grp.ls2.sort],
                    " (", ecoregion.color.table$n.site[ecoregion.color.table$ecoregion %in% grp.ls2.sort], ")"),
    fill = ecoregion.color.table$col2[ecoregion.color.table$ecoregion %in% grp.ls2.sort],
    border = NA,
    bty = "n",
    cex = 1,
    ncol = 2
  )
  mtext(
    side = 1,
    "Hour of Day",
    outer = F,
    line = 2.5,
    cex = 1.2
  )
  axis(
    side = 1,
    at = seq(0, 576, by = 6),
    labels = c(rep(c(0, 6, 12, 18), 24), 0),
    cex.axis = 0.7
  )
  axis(1, at = seq(0, 24 * 24, by = 24), labels = FALSE)
  
  
  dev.off()
  
  print(paste("################################################"))
  print(paste("                                                "))
  
}
