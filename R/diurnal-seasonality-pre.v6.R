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
base.in <- "D:\\AmeriFlux-Data\\00_olaf_data_ameriflux\\BASE-BADM\\"
badm.in <- "D:\\AmeriFlux-Data\\00_olaf_data_ameriflux\\BADM\\"
path.out.root <- paste0(path, "diurnal-seasonality\\")
RDir <- paste0(path, "R\\")

source(paste0(RDir, "math_utility.R"))
source(paste0(RDir, "get_full_list.R"))
source(paste0(RDir, "grab_version.R"))
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
ver <- "20231019-8h5d"     # for storing outputs

if (!dir.exists(paste(path.out.root, ver, sep = "")))
  dir.create(paste(path.out.root, ver, sep = ""))
path.out <- paste(path.out.root, ver, "\\", sep = "")

### Target variables to generate diurnal-seasonal summary
# sel.var <- c("FC", "SC", "NEE", "LE", "H", "FCH4",
#              "USTAR", "NETRAD", "TA", "VPD", "WS", 
#              "SWC", "SW_IN", "LW_IN", "LW_OUT", "G", 
#              "SW_OUT", "WTD", "P" 
# )
sel.var <- c("FC", "SC", "NEE", "LE", "H", "USTAR",
             #"FCH4", "SW_IN", "WTD", "P",
             "NETRAD", "TA", "VPD", "SWC")

### specify the sites don't run
drop.ls <- c("US-GBT", "US-PFa", ## unique site
             "CR-Lse" # no flux
              ) 

## following sites have better coverage in FC+SC than NEE. Use FC + SC instead
skip.nee.ls <- c("PR-xGU", "PR-xLA", "US-Ivo", "US-xBA", "US-xBR",
                 "US-xLE", "US-xMB", "US-xNW", "US-xSL", "US-xWR")
## short tower w/ SC, but poor data coverage, dropping SC
skip.sc.ls <- c("US-Ivo", "US-xNW", "US-xBA")

### define the window for calculating diurnal-seasonal summary
l.wd <- 5                    # length of days
n.wd <- floor(365 / l.wd)    # number of windows per year
min.year.to.run <- 3         # minimum data record to include the site, except the follows
border.year.site <- c("US-CdM", "US-LA2", "AR-TF1", "AR-TF2", "US-INe", "US-MWA", "US-MWW") 
forest.igbp.ls <- c("EBF", "ENF", "DBF", "DNF", "MF", "WSA")

## output time resolution
#   org: maintain original time resolution (HH/HR) in generating diurnal stat
#   aggregate:  convert half-hour resolution to hourly
hr.res <-  "aggregate"
out.res <- 3 # output resolution, in hours

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
                             badm.file = "AMF_AA-Net_BIF_LEGACY_20230922.xlsx")

sgi.amf <- merge.data.frame(sgi.amf,
                            utc.offset[, -1],
                            by = "SITE_ID", 
                            all = TRUE)

sgi.amf$DATA_YEAR <- sgi.amf$DATA_END - sgi.amf$DATA_START + 1

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
 target.site <- get_full_list(base.in)   
##  or specified subset for update only, by c(CC-XXX,.....)
# target.site<-c("BR-Npw","CA-BOU","CA-DB2","CA-DBB",
#                "CA-DSM","CA-HPC","CA-SCB","CA-SCC","CA-SMC",
#                "CA-TVC","PE-QFR","US-A03","US-A10","US-Act",
#                "US-ALQ","US-Bi1","US-Bi2","US-BZB","US-BZF",
#                "US-BZo","US-BZS","US-Cms","US-CRT","US-DFC",
#                "US-Dmg","US-DPW","US-EDN","US-EML","US-HB1",
#                "US-Ho1","US-HRA","US-HRC","US-Hsm","US-ICh",
#                "US-ICs","US-ICt","US-Ivo","US-JRn","US-KPL",
#                "US-LA1","US-LA2","US-Los","US-MBP","US-Myb",
#                "US-NC2","US-NC4","US-NGB","US-NGC","US-ONA",
#                "US-ORv","US-OWC","US-PAS","US-PFa","US-Sag",
#                "US-Sne","US-Snf","US-Srr","US-StJ","US-TrB",
#                "US-TrS","US-Tw1","US-Tw3","US-Tw4","US-Tw5",
#                "US-Twt","US-Uaf","US-WPT")

# drop site from drop.ls
target.site <- target.site[!target.site %in% drop.ls]

# drop sites with less than min.year.to.run
target.site <- target.site[target.site %in% sgi.amf$SITE_ID[!is.na(sgi.amf$DATA_YEAR) &
                                                              sgi.amf$DATA_YEAR >= min.year.to.run]]

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
    "US-Cwt"
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

src.list <- src.list[-which(substr(
  src.list,
  start = nchar(src.list) - 3,
  stop = nchar(src.list)
) == "None")]

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
full.ls$SC_source <- NA
full.ls$NEE_source <- NA
full.ls.col <- ncol(full.ls)
full.ls <- data.frame(full.ls,
                      matrix(0, nrow = nrow(full.ls), 
                             ncol = length(sel.var) * 2))
colnames(full.ls)[c((full.ls.col + 1):ncol(full.ls))] <-
  paste(rep(sel.var, each = 2), c("", "_F"), sep = "")


#
#fores.ls <- full.ls$SITE_ID[full.ls$IGBP %in% forest.igbp.ls]
#fores.ls[!fores.ls %in% unique(c(var.info$Site_ID, zm.list$site))]


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
    latest.two(grab.version(base.list.tmp)[[1]],
               grab.version(base.list.tmp)[[2]])[1]
  base.list <-
    paste("AMF_", target.site[k], "_BASE-BADM_", target.base.ver, sep = "")
  
  ## specify the BASE files under zip to grab
  case.ls <- NULL
  for (j in 1:length(target.base.ver)) {
    file.grab <-
      unzip(paste(base.in, base.list[j], ".zip", sep = ""), list = TRUE)[, "Name"]
    case.ls <-
      c(case.ls, file.grab[which(
        substr(file.grab, start = 12, stop = 15) == "BASE" &
          substr(file.grab, start = 17, stop =
                   18) == target.res[k]
      )])
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
      skip = 2,
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
  
  ## work on exceptional cases
  # only filled TA provided for US-PHM
  if(target.site[k] %in% c("US-PHM")){
    colnames(data.in)[which(colnames(data.in) == "TA_PI_F")] <- "TA"
  }
  
  var.name <- names(data.in)
  
  # convert VPD from kPa to hPa
  if (substr(
    target.base.ver,
    start = nchar(target.base.ver) - 1,
    stop = nchar(target.base.ver)
  ) == "-1" |
  target.site[k] %in% VPD.need.correct) {
    data.in[, which(substr(var.name, start = 1, stop = 3) == "VPD")] <-
      data.in[, which(substr(var.name, start = 1, stop = 3) == "VPD")] *
      10
  }
  
  #################################################################################################################
  ## Parses variable names/qualifiers
  
  ## decode the variable name, info of aggregation level, gap-filling,...
  basename_decode <- basename_parse(var.name = var.name,
                                    FP.ls = FP.ls$Name,
                                    echo = T)

  # output variable name parsed results
  write.csv(
    basename_decode,
    paste0(path.out, case.name.short, "_basename_decode.csv"),
    quote = F,
    row.names = F
  )
  
  #################################################################################################################
  ## Filter data based on the physical range
  data.work <- filter_physical_range(
    data.in = data.in,
    basename_decode = basename_decode,
    limit.ls = FP.ls,
    echo = T,
    loose.filter = 0.1
  )
  
  ## drop LE data from Us-Cav as many large negative values
  if(target.site[k] == "US-CaV"){
    data.work$LE <- NA
  }
  
  #################################################################################################################
  ## Handle the aggregation, ONLY on non-filled variables
  
  #### I. handle the equal/naive aggregation, i.e., only renaming var
  #   equal: unique variable, no aggregation needed,
  #             -> use base name
  #             e.g., ignore any _1_1_1 or _1 if provided
  #   naive: already layer-aggregated & not unique,
  #             -> use provided layer index
  
  data.work.holder <-
    rename_aggregate(basename_decode = basename_decode,
                     data.in = data.work)
  
  #### II. handle the 1step/2step aggregation
  #   1step: already replicate-averaged & not unique,
  #             -> need to aggregate to layer
  #   2step: quadruplet provided & not unique,
  #             -> need 1 or 2 steps aggregation
  
  data.work.holder <-
    simple_aggregate(basename_decode.work = data.work.holder[[2]],
                     data.work = data.work.holder[[1]],
                     FP.ls = FP.ls$Name)
  
  #### III. handle manual aggregation (not implemented)
  #    manual: mix of aggregation levels
  #              -> need manual inspection before planning for aggregation
  
  data.work <- data.work.holder[[1]]
  #basename_decode.work=data.work.holder[[2]]
  
  #################################################################################################################
  ## Deal with time stamps
  # Working on TIMESTAMP index
  data.work <- data.frame(
    TIMESTAMP = NA,
    DATE_ID = NA,
    HOUR_ID = NA,
    data.work
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
  
  full.ls$data_year_length[k] <- round(nrow(data.work) / d.hr.org / 365, digits = 1)
  
  if(full.ls$data_year_length[k] < min.year.to.run &
     !full.ls$SITE_ID[k] %in% border.year.site){
    
    print(paste("[Warning] Skip the site as less than", min.year.to.run , "year record"))
    
  } else {
    
    #############################################################################################################
    ## calculate potential radiation based on site's location & time zone,
    #  this is mainly used fo day/night separation in ustar threshold determination
    SW_IN_POT <- ext_radiation_v2(
      latitude = full.ls$LOCATION_LAT[k],
      longitude = full.ls$LOCATION_LONG[k],
      year = c(data.work$TIMESTAMP$year[1]+1900,
               data.work$TIMESTAMP$year[nrow(data.work)]+1900),
      utc.offset = full.ls$UTC_OFFSET[k],
      res = target.res[k]
    )
    
    SW_IN_POT <-
      SW_IN_POT[c(
        which(SW_IN_POT$TIMESTAMP_END == data.work$TIMESTAMP_END[1]):
          which(SW_IN_POT$TIMESTAMP_END == data.work$TIMESTAMP_END[nrow(data.work)])
      ), ]
    
    data.work <- merge.data.frame(data.work,
                                  SW_IN_POT[, -which(colnames(SW_IN_POT) %in% c("TIMESTAMP_START", "TIMESTAMP"))],
                                  by = "TIMESTAMP_END",
                                  all.y = T)
    
    
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
    
    basename_decode.work <-
      basename_parse(var.name = colnames(data.work),
                     FP.ls = FP.ls$Name,
                     echo = F)
    
    ############################################################################################################################
    ## Deal with RH and/or VPD
    vpd.loc <- which(basename_decode.work$basename == "VPD" &
                       !basename_decode.work$is_gapfill)
    rh.loc <- which(basename_decode.work$basename == "RH" &
                      !basename_decode.work$is_gapfill)
    ta.loc <- c(
      which(
        basename_decode.work$basename == "TA" &
          !basename_decode.work$is_gapfill
      ),
      which(
        basename_decode.work$basename == "T_SONIC" &
          !basename_decode.work$is_gapfill
      )
    )
    
    if (length(vpd.loc) > 0) {
      print(paste("[Info] using VPD"))
      
    } else if (length(rh.loc) > 0 & length(ta.loc) > 0) {
      data.work$VPD <-
        (0.612 * exp(17.27 * data.work[, ta.loc[1]] / (data.work[, ta.loc[1]] +
                                                         237.3))) *
        (1 - data.work[, rh.loc[1]] / 100) * 10
      
      print(paste("[Info] using RH"))
    } else{
      #data.work$e<-1 ## dummy
      print("[Warning] no moisture variables")
    }
    
    ############################################################################################################################
    ## Deal with P
    if("P" %in% sel.var){
      p.loc <- which(basename_decode.work$basename == "P" &
                       !basename_decode.work$is_gapfill)
      
      if (length(p.loc) > 0) {
        print(paste("[Info] calculate cumulative P"))
        
        YEAR_ID <- data.work$TIMESTAMP$year - min(data.work$TIMESTAMP$year) + 1
        
        data.work[is.na(data.work[, p.loc]), p.loc] <- 0
        for(yy in 1:max(YEAR_ID)){
          data.work[YEAR_ID == yy, p.loc] <- cumsum(data.work[YEAR_ID == yy, p.loc])
        }
      }
    }
    
    ############################################################################################################################
    ## Deal with PPFD_IN or SW_IN
    sw.loc <- which(basename_decode.work$basename == "SW_IN" &
                      !basename_decode.work$is_gapfill)
    ppfd.loc <- which(basename_decode.work$basename == "PPFD_IN" &
                        !basename_decode.work$is_gapfill)
    
    if (length(sw.loc) > 0) {
      print(paste("[Info] using SW_IN"))
      
    } else if (length(ppfd.loc) > 0) {
      data.work$SW_IN <- data.work[, ppfd.loc[1]] * 0.45
      sw.loc <- which(colnames(data.work) == "SW_IN")
      
      print(paste("[Info] using PPFD_IN"))
    } else{
      #data.work$e<-1 ## dummy
      print("[Warning] no radiation variables")
    }
    
    #################################################################################################################
    ## Handle ustar-filtering FC, using REddyProc function 
    
    fc.loc <- which(basename_decode.work$basename == "FC" &
                      !basename_decode.work$is_gapfill)[1]
    ustar.loc <- which(basename_decode.work$basename == "USTAR" &
                         !basename_decode.work$is_gapfill)[1]
    sw_pot.loc <- which(basename_decode.work$variable_names == "SW_IN_POT")
    
    if (sum(!is.na(fc.loc)) > 0 & sum(!is.na(ustar.loc)) > 0 &
        sum(!is.na(sw_pot.loc)) > 0 & sum(!is.na(ta.loc)) > 0
    ) {
      uStarTh <-
        suppressWarnings(suppressMessages(
          get_ustar_threshold(
            data.work = data.work,
            fc.loc = fc.loc,
            sw_pot.loc = sw_pot.loc,
            ta.loc = ta.loc[1],
            ustar.loc = ustar.loc,
            lat = as.numeric(full.ls$LOCATION_LAT[k]),
            long  = as.numeric(full.ls$LOCATION_LONG[k]),
            d.hr.org = d.hr.org
          )
        ))
      ## filtering using yearly thresholds
      year.ls <-
        uStarTh$seasonYear[uStarTh$aggregationMode == "year"]
      
      # catch exception for failed get_ustar_threshold()
      if (is.null(year.ls)) {
        data.work[!is.na(data.work[, ustar.loc]) &
                    data.work[, ustar.loc] < 0.2,
                  fc.loc] <- NA
        
      } else if(is.na(uStarTh$uStar[which(uStarTh$aggregationMode == "single")])){ 
        data.work[!is.na(data.work[, ustar.loc]) &
                    data.work[, ustar.loc] < 0.2,
                  fc.loc] <- NA
        
      } else{
        for (uu in seq_len(length(year.ls))) {
          ## use yearly threshold
          uStarTh.use <-
            uStarTh$uStar[which(uStarTh$aggregationMode == "year" &
                                  uStarTh$seasonYear == year.ls[uu])]
          
          ## if yearly threshold unavailable, use cross-year threshold
          if (is.na(uStarTh.use))
            uStarTh.use <-
              uStarTh$uStar[which(uStarTh$aggregationMode == "single")]
          
          data.work[data.work$TIMESTAMP$year + 1900 == year.ls[uu] &
                      !is.na(data.work[, ustar.loc]) &
                      data.work[, ustar.loc] < uStarTh.use,
                    fc.loc] <- NA
          
        }
      }
    }
    
    #################################################################################################################
    ## Work on the storage correction for FC -> NEE
    sc.loc <- which(basename_decode.work$basename == "SC" &
                      !basename_decode.work$is_gapfill)
    nee.loc <- which(basename_decode.work$basename == "NEE" &
                       !basename_decode.work$is_gapfill)
    co2.loc<-which(basename_decode.work$basename=="CO2"&
                     !basename_decode.work$is_gapfill)
    # pressure, used for SC calculation
    pa.loc <- which(basename_decode.work$basename == "PA" &
                      !basename_decode.work$is_gapfill)
    
    if (full.ls$IGBP[k] %in% forest.igbp.ls &
        (length(nee.loc) == 0 | 
         target.site[k] %in% skip.nee.ls)) {
      
      ## check if SC provided
      if(length(sc.loc) == 0){
        
        ## check top level CO2
        if (length(co2.loc) > 0 &
            sum(
              basename_decode.work$aggregate_method[co2.loc] == "naive" |
              basename_decode.work$aggregate_method[co2.loc] == "equal"
            ) == length(co2.loc)) {
          if (length(co2.loc) > 1) {
            co2.top <-
              co2.loc[which(as.numeric(basename_decode.work$layer[co2.loc]) ==
                              min(as.numeric(basename_decode.work$layer[co2.loc])))]
          } else{
            co2.top <- co2.loc
          }
          
          ## get EC measurement height
          site.zm.list <-
            work_on_zm_list3(
              work.path = paste0(path, "R\\"),
              target.site = target.site[k],
              var.info = var.info,
              FP.ls = FP.ls$Name
            )
          
          ## work on single level SC
          if (sum.notna(site.zm.list$htower) > 0 &
              max(site.zm.list$htower) > 5) {
            
            # parse CO2/EC height
            if (nrow(site.zm.list) == 1) {
              zec <- rep(site.zm.list$htower, nrow(data.work))
            } else{
              zec.break <- 0
              zec <- NULL
              for (k2 in 1:(nrow(site.zm.list) - 1)) {
                zec.break <-
                  c(zec.break, which(
                    as.character(data.work$TIMESTAMP_END) == as.character(site.zm.list$end[k2])
                  ))
                zec <-
                  c(zec, rep(site.zm.list$htower[k2], (zec.break[k2 + 1] - zec.break[k2])))
              }
              zec <-
                c(zec, rep(site.zm.list$htower[nrow(site.zm.list)], nrow(data.work) - length(zec)))
            }
            
            # calculate SC from top CO2
            if(length(pa.loc) == 1) {
              pa.data <- data.work[, pa.loc]
              
            } else if (!is.na(full.ls$LOCATION_ELEV[k])) {
              pa.data <-
                101.3 * exp(-full.ls$LOCATION_ELEV[k] / 8200)
              
            } else{
              pa.data <- 101.3
            }
            
            if(length(ta.loc) >= 1){
              ta.data <- data.work[, ta.loc[1]]
              
            } else{
              ta.data <- 15
              
            }
            
            data.work$SC <-
              c(NA, 
                0.5 * (data.work[c(3:(nrow(data.work))), co2.top] - 
                         data.work[c(1:(nrow(data.work) - 2)), co2.top]),
                NA) * 
              pa.data/ 
              (ta.data + 273.15) / 
              8.3143 * 1000 * zec / (hr * 60)
            
            # filter out-of-range SC
            data.work$SC[!is.na(data.work$SC) & abs(data.work$SC) > 50] <- NA
            
            print("[Info] forest/no SC, calculate SC from top CO2")
            full.ls$SC_source[k] <- "CO2"
            
          } else if(sum.notna(site.zm.list$htower) > 0 &
                    max(site.zm.list$htower) <= 5){
            data.work$SC <- 0
            print("[Info] short tower/no SC, assume SC = 0")
            full.ls$SC_source[k] <- "negligible"
            
          } else{
            print("[error] forest tower/no SC, missing EC height")
            full.ls$SC_source[k] <- "failed"
            
          }
        }
        
      } else{
        print("[Info] forest with SC")
        full.ls$SC_source[k] <- "SC"
        
      }
      
    } else if(!full.ls$IGBP[k] %in% forest.igbp.ls & 
              (length(nee.loc) == 0 | 
               target.site[k] %in% skip.nee.ls)){
      ## short vegetation, no NEE
      
      ## check if SC provided
      if (length(sc.loc) == 0 | 
          target.site[k] %in% skip.sc.ls) {
        data.work$SC <- 0
        print("[Info] short tower/no SC, assume SC = 0")
        full.ls$SC_source[k] <- "negligible"
        
      } else{
        print("[Info] short tower with SC")
        full.ls$SC_source[k] <- "SC"
        
      }
    }
    
    basename_decode.work <-
      basename_parse(var.name = colnames(data.work),
                     FP.ls = FP.ls$Name,
                     echo = F)
    
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
      sel.var.loc <- which(basename_decode.work$basename == sel.var[j1])
      
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
    
    ### handle NEE
    if (length(which(basename_decode.work$basename == "NEE")) == 0 | 
        target.site[k] %in% skip.nee.ls) {
      full.ls[k, c("NEE")] <- min(full.ls[k, c("FC")], full.ls[k, c("SC")])
      full.ls[k, c("NEE_F")] <- min(full.ls[k, c("FC_F")], full.ls[k, c("SC_F")])
      
      nee.loc2 <- which(colnames(data.out3) == "NEE")
      fc.loc2 <- which(colnames(data.out3) == "FC")
      sc.loc2 <- which(colnames(data.out3) == "SC")
      
      #data.out1[,nee.loc2]<-data.out1[,fc.loc2]+data.out1[,sc.loc2]
      data.out2[, nee.loc2] <- data.out2[, fc.loc2] + data.out2[, sc.loc2]
      data.out3[, nee.loc2] <- data.out3[, fc.loc2] + data.out3[, sc.loc2]
      data.out4[, nee.loc2] <- data.out4[, fc.loc2] + data.out4[, sc.loc2]
      #data.out5[,nee.loc2]<-data.out5[,fc.loc2]+data.out5[,sc.loc2]
      
      print(paste("[Info] generate NEE diurnal-seasonal from FC+SC"))
      full.ls$NEE_source[k] <- "FC+SC"
      
    } else{
      print(paste("[Info] generate NEE diurnal-seasonal from NEE"))
      full.ls$NEE_source[k] <- "NEE"
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
  paste(path.out, "ALL_BASE_site_list.csv", sep = ""),
  quote = F,
  na = "-9999",
  row.names = F
)

## distribution of data availability
hist(apply(full.ls[, c("NEE_F", "LE_F", "H_F", "USTAR_F", 
                       "NETRAD_F", "VPD_F", "TA_F", "SWC_F")], 1, function(x) sum(x == 1)))

## select sites with 2+ flux variables, 4+ met variables that have no gaps
select.ls <-
  full.ls[apply(full.ls[, c("NETRAD_F", "VPD_F", "TA_F", "SWC_F", 
                            "NEE_F", "LE_F", "H_F", "USTAR_F")],
                  1, function(x) sum(x == 1)) > 0, ]

write.csv(
  select.ls,
  paste(path.out, "ALL_BASE_site_short_list.csv", sep = ""),
  quote = F,
  na = "-9999",
  row.names = F
)

if (sink.log.to.file) {
  sink(NULL)
}
