#### v6.2 handle FLUXNET version
#### v6 replace with amerifluxr for handling data / API
#### v4 # start implementation for decoding qualifier, work with old/new BASE
# capable of down-scaling HH to HR
#### v3 # search for the latest BASE, based on olaf BASE-BADM folders
# don't generate anything from newer BASE X-2, which has slightly diff header (no need to skip first 2 lines)
#                                              & diff variables
# generate both HH and HR version if both BASE exist

#### v2 #-1 version has VPD in hPa, need to convert it to hPa

rm(list = ls()) # clean working memory

library(zoo)
library(stringr)
library(httr)
library(REddyProc)
library(readxl)
library(jsonlite)
library(amerifluxr)

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
base.in <- "D:\\AmeriFlux-Data\\ALL_FLUXNET\\"
badm.in <- "D:\\AmeriFlux-Data\\00_olaf_data_ameriflux\\BADM\\"
path.out.root <- paste0(path, "diurnal-seasonality\\")
sum.out <- paste0(path, "summary\\20231019\\")
RDir <- paste0(path, "R\\")

source(paste0(RDir, "math_utility.R"))
source(paste0(RDir, "get_full_list3.R"))
source(paste0(RDir, "grab_version3.R"))
source(paste0(RDir, "basename_parse.R"))
source(paste0(RDir, "rename_aggregate.R"))
source(paste0(RDir, "simple_aggregate.R"))
source(paste0(RDir, "filter_physical_range.R"))
source(paste0(RDir, "plot_diurnal_seasonal.R"))
#source(paste0(RDir, "get_latest_VI.R"))
source(paste0(RDir, "work_on_zm_list3.R"))
source(paste0(RDir, "get_ustar_threshold.R"))
source(paste0(RDir, "ext_radiation_v2.R"))
source(paste0(RDir, "badm.extract.R"))
source(paste0(RDir, "get_utc_offset.R"))


### Workflow control
sink.log.to.file <- F     # Write warning/error messages to file

###   Create a version sub directory
ver <- "20241004-GPP"     # for storing outputs

if (!dir.exists(paste(path.out.root, ver, sep = "")))
  dir.create(paste(path.out.root, ver, sep = ""))
path.out <- paste(path.out.root, ver, "\\", sep = "")

### Target variables to generate diurnal-seasonal summary
sel.var <- c("TA_F", "SW_IN_F", "LW_IN_F", "VPD_F",
             "LE_F_MDS", "H_F_MDS", "NEE_VUT_REF",
             "RECO_NT_VUT_REF", "GPP_NT_VUT_REF",
             "RECO_DT_VUT_REF", "GPP_DT_VUT_REF")

### specify the sites don't run
drop.ls <- NULL

### define the window for calculating diurnal-seasonal summary
l.wd <- 7                    # length of days
n.wd <- floor(365 / l.wd)    # number of windows per year
min.year.to.run <- 2         # minimum data record to include the site

## output time resolution
#   org: maintain original time resolution (HH/HR) in generating diurnal stat
#   aggregate:  convert half-hour resolution to hourly
hr.res <-  "aggregate"
out.res <- 2 # output resolution, in hours

## list of sites need change VPD unit, other than #-1 version
VPD.need.correct <- NULL

## min percentage of non-missing in filling diurnal-seasonal time series, within each window
allow.nongap.within.wd <- 0.80  
## min percentage of non-missing in filling diurnal-seasonal time series, across window
allow.nongap.across.wd <- 0.80  


##############################################################################################
## Site general info
sgi.amf <- amerifluxr::amf_site_info()

## UTC offset
utc.offset <- get_utc_offset(site = sgi.amf$SITE_ID,
                             badm.path = badm.in,
                             badm.file = "AMF_AA-Net_BIF_LEGACY_20230629.xlsx")

sgi.amf <- merge.data.frame(sgi.amf,
                            utc.offset[, -1],
                            by = "SITE_ID", 
                            all = TRUE)

rownames(sgi.amf) <- paste(sgi.amf$SITE_ID)

##############################################################################################
### FP-Standard variable list
### Predefined physical range
FP.ls <- amerifluxr::amf_variables()

########################################################################
## get the latest measurement height from ftp
var.info <- amerifluxr::amf_var_info()

#############################################################################################
### Specify the target site list
## could be full list (initial pull out), pull from get_full_list function
target.site <-
  get_full_list3(base.in,
                 target = "SUBSET")
##  or specified subset for update only, by c(CC-XXX,.....)

##
last.site <- read.csv(paste0(sum.out, "ALL_BASE_site_short_list4.csv"),
                      header = T)

# drop site not in previous run
target.site <- target.site[which(target.site %in% last.site$SITE_ID)]

# drop sites with less than min.year.to.run
# target.site <- target.site[target.site %in% sgi.amf$SITE_ID[!is.na(sgi.amf$DATA_YEAR) &
#                                                               sgi.amf$DATA_YEAR >= min.year.to.run]]

target.res <- rep("HH", length(target.site))
target.res[which(
  target.site %in% c(
    "BR-Sa1",
    "US-MMS",
    "US-Ha1",
    "US-Cop",
    "US-Ne1",
    "US-Ne2",
    "US-Ne3",
    "US-PFa",
    "US-Cwt",
    "US-UMB"   ### specifici for FLUXNET2015, UMB in HR
  )
)] <- "HR"

### Search within BASE-BADM files, looking for latest version
src.list <- list.files(base.in)
src.list <-
  src.list[which(substr(
    src.list,
    start = nchar(src.list) - 3,
    stop = nchar(src.list)
  ) == ".zip")]
src.list <- substr(src.list, start = 1, stop = nchar(src.list) - 4)

#####################################################################################################
## Prepare a full site list
full.ls <- sgi.amf[target.site,
                   c(
                     "SITE_ID",
                     "IGBP",
                     "LOCATION_LAT",
                     "LOCATION_LONG",
                     "LOCATION_ELEV",
                     "MAT",
                     "MAP",
                     "CLIMATE_KOEPPEN",
                     "UTC_OFFSET"
                   )]
full.ls$LOCATION_LAT <- as.numeric(full.ls$LOCATION_LAT)
full.ls$LOCATION_LONG <- as.numeric(full.ls$LOCATION_LONG)
full.ls$LOCATION_ELEV <- as.numeric(full.ls$LOCATION_ELEV)
full.ls$UTC_OFFSET <- as.numeric(full.ls$UTC_OFFSET)
full.ls$data_year_length <- NA
#full.ls$SC_source <- NA
#full.ls$NEE_source <- NA
full.ls.col <- ncol(full.ls)
full.ls <- data.frame(full.ls,
                      matrix(0, nrow = nrow(full.ls), 
                             ncol = length(sel.var) * 2))
colnames(full.ls)[c((full.ls.col + 1):ncol(full.ls))] <-
  paste(rep(sel.var, each = 2), c("", "_F"), sep = "")


if (sink.log.to.file) {
  sink(file = paste(
    path.out,
    substr(Sys.time(), start = 1, stop = 10),
    "-diurnal-seasonal-pre-log.txt",
    sep = ""
  ))
}

######################################################################################################
for (k in 1:length(target.site)) {
  
  ## control var, used in time-handling the original data file
  hr <- ifelse(target.res[k] == "HH", 30, 60)
  ## control var, original data file
  d.hr.org <- ifelse(target.res[k] == "HH", 48, 24)
  
  ## control var, used in preparing the diurnal-seasonal outputs
  d.hr <- ifelse(hr.res == "org",
                 ifelse(target.res[k] == "HH", 48, 24),
                 48 / out.res / 2)  
  
  #################################################################################################################
  ## Deal with file names, read in files
  
  ## locate the latest BASE version
  base.list.tmp <-
    src.list[which(substr(src.list, start = 5, stop = 10) == target.site[k])]
  
  target.base.ver <-
    latest.two(grab.version3(base.list.tmp, target = "SUBSET")[[1]],
               grab.version3(base.list.tmp, target = "SUBSET")[[2]])[1]
  base.list <-
    base.list.tmp[grep(target.base.ver, base.list.tmp)]
  
  ## specify the BASE files under zip to grab
  case.ls <- NULL
  for (j in 1:length(target.base.ver)) {
    file.grab <-
      unzip(paste(base.in, base.list[j], ".zip", sep = ""), list = TRUE)[, "Name"]
    case.ls <-
      c(case.ls, file.grab[which(grepl("SUBSET", file.grab) &
                                   grepl(target.res[k], file.grab))
      ])
  }
  case.ls <- substr(case.ls, start = 1, stop = nchar(case.ls) - 4)
  
  ## newest BASE
  case.name <- base.list[1]
  data.in <-
    read.table(
      unz(
        paste(base.in, case.name, ".zip", sep = ""),
        paste(case.ls[1], ".csv", sep = "")
      ),
      na = c("-9999", "-9999.0", "-9999.00", "-9999.000", "-9999.0000"),
      header = T,
      sep = ",",
      skip = 0,
      stringsAsFactors = F
    )
  
  ## used for outputs
  case.name.short <-
    paste(target.site[k], "_", target.res[k], "_", target.base.ver, sep = "")
  
  print(paste("################################################"))
  print(paste("######  ", target.site[k], "  ##################"))
  
  #################################################################################################################
  ## This part cleans empty columns, VPD unit
  # drop all empty columns
  all.empty <-
    which(apply(data.in[, c(3:ncol(data.in))], 2, sum.notna) == 0) + 2
  if (length(all.empty) > 0) {
    data.in <- data.in[, -all.empty]
  }
  
  # check if variable misinterpreted as characters
  for (l in 3:ncol(data.in)) {
    if (!is.numeric(data.in[, l])) {
      data.in[, l] <- as.numeric(paste(data.in[, l]))
      print(paste("[warning] force ", colnames(data.in)[l], "to numeric"))
    }
  }

 var.name <- names(data.in)

  #################################################################################################################
  ## Deal with time stamps
  # Working on TIMESTAMP index
  data.work <- data.frame(
    TIMESTAMP = NA,
    DATE_ID = NA,
    HOUR_ID = NA,
    data.in
  )
  
  if (length(which(colnames(data.work) == "TIMESTAMP_START")) == 1) {
    data.work$TIMESTAMP <-
      strptime(data.work$TIMESTAMP_START,
               format = "%Y%m%d%H%M",
               tz = "UTC")
    data.work$TIMESTAMP <-
      strptime(data.work$TIMESTAMP + 0.5 * hr * 60,
               format = "%Y-%m-%d %H:%M:%S",
               tz = "UTC")
  } else if (length(which(colnames(data.work) == "TIMESTAMP_END")) == 1) {
    data.work$TIMESTAMP <-
      strptime(data.work$TIMESTAMP_END,
               format = "%Y%m%d%H%M",
               tz = "UTC")
    data.work$TIMESTAMP <-
      strptime(data.work$TIMESTAMP - 0.5 * hr * 60,
               format = "%Y-%m-%d %H:%M:%S",
               tz = "UTC")
    print("[Warning] can't find TIMESTAMP_START columns")
  } else{
    print("[Error] can't find both TIMESTAMP columns")
  }
  
  full.ls$data_year_length[k] <- floor(nrow(data.work) / d.hr.org / 365)
  
  if(full.ls$data_year_length[k] < min.year.to.run){
    
    print(paste("[Warning] Skip the site as less than", min.year.to.run , "year record"))
    
  } else {
    
    
    #################################################################################################################
    ## Create holder for working diurnal-seasonal outputs
    # Create time id for diurnal-seasonal aggregation, apply function
    data.work$DATE_ID <- floor(data.work$TIMESTAMP$yday / l.wd) + 1
    data.work$DATE_ID[which(data.work$DATE_ID == max(data.work$DATE_ID, na.rm =
                                                       T))] <- n.wd
    # wrap all hanging days into last window
    
    if (hr.res == "org") {
      data.work$HOUR_ID <-
        data.work$TIMESTAMP$hour + data.work$TIMESTAMP$min / 60
    } else if (hr.res == "aggregate") {
      data.work$HOUR_ID <-
        (floor(data.work$TIMESTAMP$hour / (24 / d.hr)) * (24 / d.hr)) + 0.5 * (24 / d.hr)
    }
    
    #################################################################################################################################
    
    ## working holder for diel-seasonal output, with length of n.wd * d.hr
    DATE_ID <- rep(c(1:(n.wd)), each = d.hr)
    
    if (hr.res == "org") {
      HOUR_ID <- rep(c(1:d.hr) / (60 / hr) - (hr / 60) / 2, times = n.wd)
    } else if (hr.res == "aggregate") {
      HOUR_ID <- rep(sort(unique(data.work$HOUR_ID)), times = n.wd)
    }
    
    #data.out1<-cbind(DATE_ID,HOUR_ID)
    data.out2 <- cbind(DATE_ID, HOUR_ID)
    data.out3 <- cbind(DATE_ID, HOUR_ID)
    data.out4 <- cbind(DATE_ID, HOUR_ID)
    #data.out5<-cbind(DATE_ID,HOUR_ID)
    
    for (j1 in 1:length(sel.var)) {
      sel.var.loc <- which(colnames(data.work) == sel.var[j1])
      
      ## working holder for diurnal-seasonal output, with length of n.wd * d.hr
      #temp.out1<-array(NA,dim=c(d.hr,n.wd))
      temp.out2 <- array(NA, dim = c(d.hr, n.wd))
      temp.out3 <- array(NA, dim = c(d.hr, n.wd))
      temp.out4 <- array(NA, dim = c(d.hr, n.wd))
      #temp.out5<-array(NA,dim=c(d.hr,n.wd))
      
      if (length(sel.var.loc) == 1) {
        ## if single/unique column
        
        ## get data availability
        full.ls[k, (full.ls.col + 2 * j1 - 1)] <-
          sum.notna(data.work[, sel.var.loc]) / nrow(data.work)
        
        for (i in 1:n.wd) {
          #temp.out1[,i]<-c(tapply(data.work[data.work$DATE_ID==i,sel.var.loc],
          #                        data.work$HOUR_ID[data.work$DATE_ID==i],
          #                        low.bd2))
          temp.out2[, i] <-
            c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc],
                     data.work$HOUR_ID[data.work$DATE_ID ==
                                         i],
                     low.bd))
          temp.out3[, i] <-
            c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc],
                     data.work$HOUR_ID[data.work$DATE_ID ==
                                         i],
                     na.median))
          temp.out4[, i] <-
            c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc],
                     data.work$HOUR_ID[data.work$DATE_ID ==
                                         i],
                     upp.bd))
          #temp.out5[,i]<-c(tapply(data.work[data.work$DATE_ID==i,sel.var.loc],
          #                        data.work$HOUR_ID[data.work$DATE_ID==i],
          #                        upp.bd2))
        }
        
      } else if (length(sel.var.loc) > 0 &
                 sum(basename_decode.work$aggregate_method[sel.var.loc] ==
                     "naive") == length(sel.var.loc)) {
        ## if all layer-aggregted
        
        ## locate the top level
        sel.var.loc.top <-
          sel.var.loc[which(basename_decode.work$layer[sel.var.loc] == "1")]
        
        ## get data availability
        full.ls[k, (full.ls.col + 2 * j1 - 1)] <-
          sum.notna(data.work[, sel.var.loc.top]) / nrow(data.work)
        
        if (length(sel.var.loc.top) == 1) {
          ## if single uniue top-layer
          for (i in 1:n.wd) {
            #temp.out1[,i]<-c(tapply(data.work[data.work$DATE_ID==i,sel.var.loc.top],
            #                        data.work$HOUR_ID[data.work$DATE_ID==i],
            #                        low.bd2))
            temp.out2[, i] <-
              c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc.top],
                       data.work$HOUR_ID[data.work$DATE_ID ==
                                           i],
                       low.bd))
            temp.out3[, i] <-
              c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc.top],
                       data.work$HOUR_ID[data.work$DATE_ID ==
                                           i],
                       na.median))
            temp.out4[, i] <-
              c(tapply(data.work[data.work$DATE_ID == i, sel.var.loc.top],
                       data.work$HOUR_ID[data.work$DATE_ID ==
                                           i],
                       upp.bd))
            #temp.out5[,i]<-c(tapply(data.work[data.work$DATE_ID==i,sel.var.loc.top],
            #                        data.work$HOUR_ID[data.work$DATE_ID==i],
            #                        upp.bd2))
          }
        } else{
          print(paste(
            "[Warning] can't locate top-layer",
            paste(basename_decode.work$variable_names[sel.var.loc], collapse =
                    " ")
          ))
        }
        
      } else if (length(sel.var.loc) > 0) {
        # catch up cases that variables exist, but not unique or aggregated to laters
        
        print(paste(
          "[Warning]",
          paste(basename_decode.work$variable_names[sel.var.loc], collapse =
                  " "),
          "can't generate diel-seasonal stats"
        ))
        
      }
      
      # iterate gap-filling by interpolation within & across windows,
      # relax allow.nongap.within.wd in the 2nd round of gap-filling
      for (i3 in 1:2) {
        ## fill the gaps within each wd
        for (i in 1:n.wd) {
          if (sum(!is.na(temp.out3[, i])) != d.hr &
              sum(!is.na(temp.out3[, i])) > allow.nongap.within.wd * d.hr /
              i3) {
            #temp.out1[,i]<-na.approx(temp.out1[,i],rule=2)
            temp.out2[, i] <- na.approx(temp.out2[, i], rule = 2)
            temp.out3[, i] <- na.approx(temp.out3[, i], rule = 2)
            temp.out4[, i] <- na.approx(temp.out4[, i], rule = 2)
            #temp.out5[,i]<-na.approx(temp.out5[,i],rule=2)
          }
        }
        ## fill the gaps across wd
        for (i2 in 1:d.hr) {
          if (sum(!is.na(temp.out3[i2, ])) != n.wd &
              sum(!is.na(temp.out3[i2, ])) > allow.nongap.across.wd * n.wd) {
            #temp.out1[i2,]<-na.approx(temp.out1[i2,],rule=2)
            temp.out2[i2, ] <- na.approx(temp.out2[i2, ], rule = 2)
            temp.out3[i2, ] <- na.approx(temp.out3[i2, ], rule = 2)
            temp.out4[i2, ] <- na.approx(temp.out4[i2, ], rule = 2)
            #temp.out5[i2,]<-na.approx(temp.out5[i2,],rule=2)
          }
        }
      }
      
      ## get data availability, post filling (approximate, since gap-filling is done at diurnal-sesonal scale)
      full.ls[k, (full.ls.col + 2 * j1)] <-
        sum.notna(c(temp.out3)) / length(temp.out3)
      
      ## need to come back on this ? how to name output variable, qualifier or not??
      out.name <- paste(sel.var[j1])
      
      #data.out1<-cbind(data.out1,tmp=c(temp.out1))
      #colnames(data.out1)[which(colnames(data.out1)=="tmp")]<-out.name
      data.out2 <- cbind(data.out2, tmp = c(temp.out2))
      colnames(data.out2)[which(colnames(data.out2) == "tmp")] <-
        out.name
      data.out3 <- cbind(data.out3, tmp = c(temp.out3))
      colnames(data.out3)[which(colnames(data.out3) == "tmp")] <-
        out.name
      data.out4 <- cbind(data.out4, tmp = c(temp.out4))
      colnames(data.out4)[which(colnames(data.out4) == "tmp")] <-
        out.name
      #data.out5<-cbind(data.out5,tmp=c(temp.out5))
      #colnames(data.out5)[which(colnames(data.out5)=="tmp")]<-out.name
    }
   
    #### plot seasonal-diurnal plots
    plot_diurnal_seasonal(
      data.median = data.out3,
      data.upp1 = data.out2,
      #data.upp2=data.out1,
      data.low1 = data.out4,
      #data.low2=data.out5,
      path.out = path.out,
      case.name = case.name,
      d.hr = d.hr,
      n.wd = n.wd,
      l.wd = l.wd,
      skip.var.loc = c(1:2),
      n.data.yr = floor(nrow(data.work) / d.hr.org /
                          365)
    )
    
    # write.csv(data.out1, file = paste0(path.out, case.name.short, "_LOWER2.csv"), quote = F, na = "-9999", row.names = F)
    # write.csv(data.out2, file = paste0(path.out, case.name.short, "_LOWER1.csv"), quote = F, na = "-9999", row.names = F)
    write.csv(
      data.out3,
      file = paste0(path.out, case.name.short, "_MEDIAN.csv"),
      quote = F,
      na = "-9999",
      row.names = F
    )
    # write.csv(data.out4, file = paste0(path.out, case.name.short, "_UPPER1.csv"), quote = F, na = "-9999", row.names = F)
    # write.csv(data.out5, file = paste0(path.out, case.name.short, "_UPPER2.csv"), quote = F, na = "-9999", row.names = F)

  }

  print(paste("################################################"))
  print(paste("                                                "))
}


##########################################################
write.csv(
  full.ls,
  paste(path.out, "ALL_FLUXNET_site_list.csv", sep = ""),
  quote = F,
  na = "-9999",
  row.names = F
)

if (sink.log.to.file) {
  sink(NULL)
}
