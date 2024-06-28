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
source(paste0(RDir, "dtwPlotTwoWay2.R"))

###### Control parameter
ver <- "20231019-12h7d"
dist.ls <- c("dtw")#, "euclidean")
use.prescribed.col <- T

len.ts <- 52 * 12  ## length of time series
len.ts.tck <- 24
# 24 * 24 for 20231019-24h15d
# 52 * 12 for 20231019-12h7d
# 73 * 8  for 20231019-8h5d

## original colnames from diurnal-seasonal files, revising the colname
sel.var <- c("FC", "SC", "NEE", "LE", "H",
             #"FCH4", "SW_IN", "WTD", "P",
             "USTAR", "NETRAD", "TA", "VPD", 
             "SWC")

drop.ls<-c("FC", "SC")

# target variables to run uni-variate clustering
target.var.ls <- c("NEE", "LE", "H", "USTAR"
                   #"NETRAD", "TA", "VPD", "SWC"
                   #"FCH4", "P", "WTD"
                   )
target.lab.ls <- list(expression(NEE~'('*mu*mole~m^{-2}~s^{-1}*')'),
                      expression(LE~'('*W~m^{-2}*')'),
                      expression(H~'('*W~m^{-2}*')'),
                      expression(USTAR~'('*m~s^{-1}*')')#,
                      # expression(NETRAD~'('*W~m^{-2}*')'),
                      # expression(TA~'('*degree~C*')'),
                      # expression(VPD~'('*hPa*')'),
                      # expression(SWC~'('*'%'*')')
                      # #expression(FCH4~'('*nmole~m^{-2}~s^{-1}*')'),
                      #expression(P~'('*mm*')'),
                      #expression(WTD~'('*m*')')
                      )
target.rng.ls <- list(seq(-35, 10, length.out = 2),
                      seq(-50, 400, length.out = 2),
                      seq(-50, 400, length.out = 2),
                      seq(0, 1, length.out = 2)
                      # seq(-300, 900, length.out = 5),
                      # seq(-40, 40, length.out = 5),
                      # seq(0, 60, length.out = 5),
                      # seq(0, 100, length.out = 5)
                      # #seq(-100, 500, length.out = 6),
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

  ###### Loop through diff distance metric

    target.dist <- dist.ls[1]
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
    # distMatrix <- proxy::dist(t(data.pre),
    #                           method = target.dist.get)
    # 
    target1 <- c("US-Oho", "US-Oho", "US-Oho", "US-Oho", "US-Oho")
    target2 <- c("US-Prr", "US-MMS", "US-KL1", "US-KL2", "US-Myb")
    
    for(i in 1: length(target2)){
      ### work on examples
      
      png(paste0(path.out.root, ver, "\\", "DTW_example1_", target.var, "_", target1[i], "_", target2[i], ".png"),
          width = 7.5, height = 4.5, units = "in", res = 300)
      par(fig = c(0.1, 0.99, 0.1, 1), mar= c(1.5, 1.5, 0.5, 0))
      dist.exp1 <- dtw::dtw(data.pre[,which(colnames(data.pre) == target1[i])],
                            data.pre[,which(colnames(data.pre) == target2[i])])
      dtwPlotTwoWay2(dtw::dtw(data.pre[,which(colnames(data.pre) == target1[i])],
                    data.pre[,which(colnames(data.pre) == target2[i])],
                    keep = T),
           #type = "two", 
           col = c('deepskyblue', 'red'),
           lwd = 1.5,
           lty = 1,
           xaxt = "n",
           yaxt = "n",
           main = "",
           #axes = F,
           ylim = target.rng.ls[[l1]])
      legend(0.1,  ifelse(target.var == "NEE", 0.4, 1),
             title = paste("DTW distance:", round(dist.exp1$distance)),
             legend = c(paste0(target1[i], " (Group ", full.ls[target1[i], paste0(target.var, "_clust_group_dtw")], ")"), 
                        paste0(target2[i], " (Group ", full.ls[target2[i], paste0(target.var, "_clust_group_dtw")], ")"), 
                        "DTW path"),
             lty = c(1, 1, 2),
             col = c('deepskyblue', 'red', 'gray'),
             bty = "n")
      
      mtext("Time Series Index", side = 1, font = 2)
      mtext(target.lab.ls[[l1]], side = 2, line = 0, font =2)
      dev.off()
      
      # dist.exp2 <- dtw::dtw(data.pre[,which(colnames(data.pre) == target1[i])],
      #                       data.pre[,which(colnames(data.pre) == target2[i])])
      # par(fig = c(0, 1, 0, 0.5), mar= c(2.5, 0.5, 2, 0), new = T)
      # plot(dtw::dtw(data.pre[,which(colnames(data.pre) == target1[i])],
      #               data.pre[,which(colnames(data.pre) == target2[i])],
      #               #dist.method = target.dist.get,
      #               keep = T),
      #      type = "two", 
      #      col = c('deepskyblue', 'red'),
      #      lwd = 1.5,
      #      lty = 1,
      #      main = paste("Global window, DTW distance:", round(dist.exp2$distance)))
      # dev.off()
    }
    
   

}
