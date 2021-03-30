################################################################################################################
####  Extraterrestrial (Potential incoming) radiation for hourly or half-hourly periods (SW_IN_POT)  ###########

# By Housen Chu, 5/10/2018
# [Input Argument]
# latitude: e.g.,41.5545 # latitude of the measurement site [degrees north of equator]  
# longitude: e.g., -83.8438 # longitude of the measurement site 
# year: e.g., 2010, c(2010:2015) # a number of a vector  
# utc.offset: e.g., -8 for Pacific standard time
# res: HH (30 min) or HR (60 min)   

# [Output Value]
# A R data.frame consist of three columns:
# TIMESTAMP_START: (YYYYMMDDHHMM) starting time the row
# TIMESTAMP_END: (YYYYMMDDHHMM) ending time of the row
# SW_IN_POT: (W/m-2) potential incoming radiation

ext_radiation_v2<-function(latitude,  
                           longitude,
                           year,
                           utc.offset=0,
                           res="HH"){
  
  hr<-ifelse(res=="HH",30,60)
  
  ## prepare full time stamps
  #  TIMESTAMP for mid-point time
  #  TIMESTAMP_START for starting time
  #  TIMESTAMP_END for ending time
  if(res=="HH"){
    TIMESTAMP<-seq.POSIXt(from=strptime(paste(year[1],"01010000",sep=""),
                                        format="%Y%m%d%H%M",tz="UTC"),
                          to=strptime(paste(year[length(year)]+1,"01010000",sep=""),
                                      format="%Y%m%d%H%M",tz="UTC"),
                          by="30 min")
  }else{
    TIMESTAMP<-seq.POSIXt(from=strptime(paste(year[1],"01010000",sep=""),
                                        format="%Y%m%d%H%M",tz="UTC"),
                          to=strptime(paste(year[length(year)]+1,"01010000",sep=""),
                                      format="%Y%m%d%H%M",tz="UTC"),
                          by="60 min")
  }
  TIMESTAMP<-strptime(TIMESTAMP[-length(TIMESTAMP)]+0.5*hr*60,
                      format="%Y-%m-%d %H:%M:%S",tz="UTC")
  
  start.time<-strptime(TIMESTAMP-0.5*hr*60,
                       format="%Y-%m-%d %H:%M:%S",tz="UTC")
  TIMESTAMP_START<-((start.time$year+1900)*10^8+
                      (start.time$mon+1)*10^6+
                      start.time$mday*10^4+
                      start.time$hour*10^2+
                      start.time$min)
  
  end.time<-strptime(TIMESTAMP+0.5*hr*60,
                     format="%Y-%m-%d %H:%M:%S",tz="UTC")
  TIMESTAMP_END<-((end.time$year+1900)*10^8+
                    (end.time$mon+1)*10^6+
                    (end.time$mday)*10^4+
                    end.time$hour*10^2+
                    end.time$min)
  
  DOY<-TIMESTAMP$yday+1
  mid.time<-(TIMESTAMP$hour+TIMESTAMP$min/60)
  
  # longitude of the centre of the local time zone [degrees west of Greenwich] 
  # t.zone = 75, 90, 105 and 120 for Eastern (-5), Central (-6), Rocky Mountain (-7), Pacific (-8) time zones
  # t.zone = 0 for Greenwich, 330 (+2) for Cairo (Egypt), and 255 (+7) for Bangkok (Thailand)
  if(utc.offset>0){
    t.zone<-360-15*utc.offset
  }else{
    t.zone<-15*(-utc.offset)  
  }
  
  l.time<-ifelse(res=="HH",0.5,1)
  latitude<-latitude/180*pi # latitudeitude (radians), 
  longitude<-ifelse(longitude>0,360-longitude,-longitude)  # longitudeitude of the measurement site [degrees west of Greenwich] 0-360 
  
  Gsc<-0.0820 # solar constant, MJ m-2 min-1
  delta<-0.409*sin(2*pi*(DOY+284)/365) # solar declination (radians)
  dr2<-1+0.033*cos(2*pi/365*DOY)  # inverse relatitudeive distance Earth-Sun
  Sea.crr<-0.1645*sin(2*2*pi*(DOY-81)/364)-0.1255*cos(2*pi*(DOY-81)/364)-0.025*sin(2*pi*(DOY-81)/364) # seasonal correction for solar time [hour].
  
  omega.r<-acos(-tan(latitude)*tan(delta))  #sunrise hour angle (radians) 
  omega.s<-(-omega.r)  #sunset hour angle (radians) 
  omega<-pi/12*((mid.time+0.06667*(t.zone-longitude)+Sea.crr)-12)  # solar time angle at midpoint of hourly or shorter period [rad],
  omega.1<-omega-pi*l.time/24 # solar time angle at beginning of period  (radians) 
  omega.2<-omega+pi*l.time/24 # solar time angle at end of period 
  
  ## clear sky radiation, # W m-2 or J m-2 s-1
  SW_IN_POT<-12*60/pi*Gsc*dr2*((omega.2-omega.1)*sin(latitude)*sin(delta)+cos(latitude)*cos(delta)*(sin(omega.2)-sin(omega.1)))*1000000/1800 
  SW_IN_POT[SW_IN_POT<0]<-0
  
  return(data.frame(TIMESTAMP_START,TIMESTAMP_END,SW_IN_POT=SW_IN_POT))  
}
