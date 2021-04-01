work_on_zm_list3<-function(work.path,    ## where a previous version of zm list exist
                          target.site,  ## which site to look for
                          var.info,      ## var info
                          FP.ls){
  
  # EC level height // FLUXNET2015 version // focus on US & CA sites
  badm.tmp.ls <-
    dir(work.path)[which(grepl("FLX_AA-Flx_BIF_HH_20200205.xlsx", dir(work.path)))]
  
  data.in <- read_excel(
    paste0(work.path, badm.tmp.ls),
    sheet = 1,
    col_names = T,
    na = "-9999"
  )
  
  data.in$GROUP_ID <- as.character(data.in$GROUP_ID)
  data.in$VARIABLE_GROUP <- as.character(data.in$VARIABLE_GROUP)
  data.in$DATAVALUE <- as.character(data.in$DATAVALUE)
  data.in$SITE_ID <- as.character(data.in$SITE_ID)
  
  zm.list <- badm.extract(data.in = data.in, sel.grp = "GRP_VAR_INFO")
  
  zm.list<-zm.list[substr(zm.list$SITE_ID,1,2) %in% c("CA","US", "MX", "BR", "AR", "CR"),]
  zm.list<-zm.list[zm.list$VAR_INFO_VARNAME %in% c("USTAR"),]
  
  zm.list$end<-NA
  zm.list<-zm.list[,c("SITE_ID","VAR_INFO_DATE","end","VAR_INFO_HEIGHT")]
  colnames(zm.list) <- c("site","start","end","htower")
  
  ## pad all to 12 digits
  zm.list$start[which(nchar(zm.list$start)==8)]<-paste0(zm.list$start[which(nchar(zm.list$start)==8)],"0000")
  
  ## fill end time stamps / if single entry
  single.zm.site<-names(table(zm.list$site))[which(table(zm.list$site)==1)]
  zm.list$end[which(zm.list$site %in% single.zm.site)]<-"all"
  
  multiple.zm.site<-data.frame(table(zm.list$site)[which(table(zm.list$site)>1)],
                               stringsAsFactors = F)
  
  for(i in 1:nrow(multiple.zm.site)){
    zm.list$end[which(zm.list$site == paste(multiple.zm.site$Var1[i]))][-multiple.zm.site$Freq[i]]<-
      zm.list$start[which(zm.list$site == paste(multiple.zm.site$Var1[i]))][-1]
    
    zm.list$end[which(zm.list$site == paste(multiple.zm.site$Var1[i]))][multiple.zm.site$Freq[i]] <- "all"
  }
  
  #### work on EC height
  ## manually adding new sites
  zm.list.new <- data.frame(site=c("US-A32","US-A74","US-NGB","US-Srr","US-StJ",
                                   "US-Rls","US-Rws","US-Rms","US-RC1","US-RC2",
                                   "US-RC3","US-RC4","US-RC5","CA-ARB","CA-ARF",
                                   "CA-SCC","US-ADR","CA-SCB",
                                   "US-CF1","US-CF2","US-CF3","US-CF4",
                                   "US-Hn3","US-Hn2","US-Hn1",
                                   "US-LA1","US-LA2",
                                   "US-PHM","US-Vcs",
                                   "US-A03","US-A10",  
                                   "US-EDN","US-NGC",
                                   "US-KS3","US-DPW",
                                   "US-xHA","US-xKZ","US-xBR","US-xSR",
                                   "US-EML","US-SdH","CA-Na1",
                                   "CA-MR3","CA-MR5","US-HRC",
                                   "CA-DBB","US-Uaf","US-Fcr",
                                   "US-ALQ","US-MtB","US-Snf",
                                   "US-xCP","US-xDL","US-xKA","US-xRM","US-xWD",
                                   "US-Wpp","CA-Ca3","US-Vcp","CA-SJ2","US-Bar"),  
                            start="all",
                            end="all",
                            htower=c(3.77,3.77,4.15,3.4,3.6,
                                     2,2,2,1.91,2.21,
                                     2.13,2.07,1.96,7.5,7.5,
                                     15.2,2,1.9,
                                     2.2,2.2,2.2,2.2,
                                     2.35,2.35,2.35,
                                     3.4,3.6,
                                     15,30,
                                     3.6,3.6,
                                     4.34,2.75,
                                     2.85,1.9,
                                     39.16,8.4,35.68,8.33,
                                     3.6,3.8,18.5,
                                     2.5,2.5,2.28,
                                     1.8,6,1.85,
                                     2.4,29.8,3.9,
                                     9.32,42.43,8.41,25.35,8.58,
                                     21,12,23.8,2,24.5))
  
  zm.list.new2<-data.frame(site=c("US-NC4","US-NC4",
                                  "US-NC3","US-NC3","US-NC3",
                                  "US-NC2","US-NC2","US-NC2","US-NC2"),
                           start=c("all",201301010000,
                                   "all",201604040000,201705010000,
                                   "all",201005200000,201204170000,201505140000),
                           end=c(201301010000,"all",
                                 201604040000,201705010000,"all",
                                 201005200000,201204170000,201505140000,"all"),
                           htower=c(28.2,33.2,
                                    3.53,6.55,9,
                                    22.5,25.6,22.4,27))
  
  zm.list.new10 <- data.frame(site = c("US-xAB", "US-xBL", "US-xBN", "US-xDJ", "US-xGR",
                                      "US-xJE", "US-xML", "US-xNW", "US-xSB", "US-xSC",
                                      "US-xSE", "US-xSJ", "US-xSP", "US-xST", "US-xTA", 
                                      "US-xTR", "US-xUK", "US-xUN", "US-xWR", "US-xYE",
                                      "PR-xGU", "PR-xLA", "US-xLE", "US-xPU", "US-xRN"),
                             start = rep("all", 20),
                             end = rep("all", 20),
                             htower = c(18.9, 8.55, 19.78, 22.5, 45.89,
                                        42.46, 29.06, 8.46, 35.81, 51.88,
                                        61.61, 39.35, 52.49, 22.36, 35.7,
                                        35.97, 35.77, 39.05, 74.18, 18.43,
                                        23.27, 8, 47.13, NA, 39.19))
  
  zm.list.new11 <- data.frame(site = c("BR-Npw", "CA-LP1", "CR-Lse", "MX-Aog", "MX-Tes",
                                       "US-BZS", "US-CMW", "US-CS2", "US-Cwt", "US-HBK",
                                       "US-LPH", "US-LS2", "US-SRS", "US-BZB"),
                             start = rep("all", 13),
                             end = rep("all", 13),
                             htower = c(20, NA, 42, 14.7, NA,
                                        2.5, 14, NA, NA, 33.2,
                                        20.5, 8, 7, 2.5))
  
  zm.list.new12 <- data.frame(site = c("US-Ha2"),
                              start = c("all", 200611010000, 201406040000),
                              end = c(200611010000, 201406040000, "all"),
                              htower = c(28, 29, 33.5))

  # zm.list.new2<-data.frame(site=c("US-Tw5","US-Tw5","US-OWC","US-OWC"),
  #                          start=c("all",201806281200,"all",201708011200),
  #                          end=c(201806281200,"all",201708011200,"all"),
  #                          htower=c(5.82,6.08,3.2,5.2))
  # ## US-RO2 heights
  # zm.list.new3<-data.frame(site=rep("US-Ro2",11),
  #                          start=c("all",200807181200,200810011200,200907281200,200908310000,201001010000,
  #                                  201107211200,201110281200,201307151200,201310151200,201510301200),
  #                          end=c(200807181200,200810011200,200907281200,200908310000,201001010000,201107211200,
  #                                201110281200,201307151200,201310151200,201510301200,201701010000),
  #                          htower=c(1.93,2.79,1.93,2.79,1.93,1.75,
  #                                   2.42,1.75,2.79,1.75,2.00))
  # 
  # ## add 2013-2015 heights
  # zm.list.new4<-data.frame(site=rep("US-Ro1",5),
  #                          start=c(201301010000,201507091200,201507211200,201508051200,201511041200),
  #                          end=c(201507091200,201507211200,201508051200,201511041200,201701010000),
  #                          htower=c(2.05,2.95,3.48,3.75,2.05))
  # 
  # zm.list.new5<-data.frame(site=rep("US-Ro5",7),
  #                          start=c("all",201606241200,201608041200,201610191200,201807091200,201807181200,201811021200),
  #                          end=c(201606241200,201608041200,201610191200,201807091200,201807181200,201811021200,201901010000),
  #                          htower=c(2.06,2.2,3.85,2.45,3.25,3.8,2.78))
  # 
  # zm.list.new6<-data.frame(site=rep("US-Br1",5),
  #                          start=c("all",201004220000,201010020000,201105040000,201110070000),
  #                          end=c(201004220000,201010020000,201105040000,201110070000,"all"),
  #                          htower=c(2.25,2.25,2.25,5.20,2.25))
  # 
  # zm.list.new7<-data.frame(site=rep("US-Br3",5),
  #                          start=c("all",201004220000,201010020000,201105040000,201110070000),
  #                          end=c(201004220000,201010020000,201105040000,201110070000,"all"),
  #                          htower=c(2.00,5.00,2.00,2.00,2.00))
  # 
  # zm.list.new8<-data.frame(site=rep("CA-ER1",10),
  #                          start=c(201501010000,201506221200,201507211200,201511021200,201601010000,201607111200,201610071200,201706291200,201707251200,201711101200),
  #                          end=c(201506221200,201507211200,201601010000,201607111200,201511021200,201610071200,201706291200,201707251200,201711101200,"all"),
  #                          htower=c(1.93,2.42,3.66,2.03,1.93,3.67,2,2.85,4.5,2.25))
  
  zm.list.new9<-data.frame(site=rep("US-Rpf",3),
                           start=c(200801010000,201307010000,201809010000),
                           end=c(201307010000,201809010000,"all"),
                           htower=c(2.6,3.9,5.6134))
  
  ## overwrite sites by above inputs
  # replace.ls <- which(
  #   zm.list$site == "US-Ro2" | 
  #     zm.list$site == "US-Ro5" |
  #     zm.list$site == "US-NC3" |
  #     zm.list$site == "US-NC4"
  # )
   
  # if (length(replace.ls) > 0)
  #   zm.list <- zm.list[-replace.ls, ]
  
  zm.list <- rbind.data.frame(
    zm.list,
    zm.list.new,
    zm.list.new2,
    #zm.list.new3,
    #zm.list.new4,
    #zm.list.new5,
    #zm.list.new6,
    #zm.list.new7,
    #zm.list.new8,
    zm.list.new9,
    zm.list.new10,
    zm.list.new11,
    zm.list.new12
  )
  
  ## use var info if not found above
  if (!target.site %in% zm.list$site) {
    ## get EC height from VI
    var.info.grab <- var.info[var.info$Site_ID == target.site, ]
    var.info.grab <- var.info.grab[!is.na(var.info.grab$Height), ]
    
    if (nrow(var.info.grab) > 0) {
      
      ec.loc <-
        which(
          grepl("USTAR", var.info.grab$Variable) |
            grepl("FC", var.info.grab$Variable) |
            grepl("LE", var.info.grab$Variable)
        )
      ### issue: don't look for H, as it also returns H2O, RH......
      
      if (length(ec.loc) == 1) {
        zm.list <- data.frame(zm.list,
                              c(target.site),
                              "all",
                              "all",
                              var.info.grab$Height[ec.loc])
        
      } else if (length(ec.loc) > 1 &
                 sd(var.info.grab$Height[ec.loc]) == 0) {
        zm.list <- rbind.data.frame(zm.list,
                                    data.frame(site = target.site,
                                               start = "all",
                                               end = "all",
                                               htower = var.info.grab$Height[ec.loc][1]))
        
      } else if (length(ec.loc) == 0) {
        #print(paste("[warning] cannot find EC height"))
        stop("[warning] cannot find EC height")
        
      } else{
        #print(paste("[warning] may have multiple EC height"))
        print(var.info.grab[ec.loc, c("Variable", "Start_Date", "Height")])
        stop(paste("[warning] may have multiple EC height"))
        
      }
      
    } else{
      #print(paste("[warning] cannot find valid VI"))
      stop("[warning] cannot find valid VI")
    }
  }
  
  zm.list.out <- zm.list[zm.list$site == target.site,]
  zm.list.out$htower <- as.numeric(zm.list.out$htower)

  return(zm.list.out)
}





