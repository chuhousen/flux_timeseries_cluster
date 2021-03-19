work_on_zm_list<-function(work.path,## where a previous version of zm list exist
                          target.site  ## which site to look for
){
  
  source("get_latest_VI.R")
  
  ########################################################################
  ## get the latest measurement height from ftp
  var.info.tmp<-get_latest_VI(path.out=work.path)
  var.info.ver<-paste("_",var.info.tmp[[2]],sep="")
  var.info<-var.info.tmp[[1]]
  
  # EC level height
  zm.list<-read.csv(paste(work.path,"FLX_BADM_EC_HEIGHT-20180722.csv",sep=""),header=T,stringsAsFactors = F)   
  zm.list$start<-as.vector(zm.list$start)
  zm.list$end<-as.vector(zm.list$end)
  
  #### work on EC height
  ## overwrite existing ones
  zm.list$htower[zm.list$site=="US-NR1"]<-21.5  ## over-write, provided by Sean Burns
  zm.list$htower[zm.list$site=="US-Syv"]<-36    ## provided by Ankur
  zm.list$htower[zm.list$site=="US-GLE"]<-22.65 ## provided by Ankur
  
  ## manually adding new sites
  zm.list.new<-data.frame(site=c("US-A32","US-A74","US-NGB","US-Srr","US-StJ",
                                 "US-Rls","US-Rws","US-Rms","US-RC1","US-RC2",
                                 "US-RC3","US-RC4","US-RC5","CA-ARB","CA-ARF",
                                 "CA-SCC","US-ADR","CA-SCB",
                                 "US-CF1","US-CF2","US-CF3","US-CF4",
                                 "US-Hn3","US-Hn2","US-Hn1",
                                 "US-LA1","US-LA2",
                                 "US-PHM","US-Vcs",
                                 "US-A03","US-A10",  ## place holding for these two sites, need update
                                 "US-EDN","US-NGC",
                                 "US-KS3","US-DPW",
                                 "US-xHA","US-xKZ","US-xBR","US-xSR",
                                 "US-EML","US-SdH","CA-Na1",
                                 "CA-MR3","CA-MR5","US-HRC",
                                 "CA-DBB","US-Uaf","US-Fcr",
                                 "US-ALQ","US-MtB","US-Snf",
                                 "US-xCP","US-xDL","US-xKA","US-xRM","US-xWD",
                                 "US-Wpp"),  
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
                                   21))
  
  zm.list.new2<-data.frame(site=c("US-NC4","US-NC4","US-NC3","US-NC3","US-NC3","US-Tw5","US-Tw5","US-OWC","US-OWC"),
                           start=c("all",201301010000,"all",201604040000,201705010000,"all",201806281200,"all",201708011200),
                           end=c(201301010000,"all",201604040000,201705010000,"all",201806281200,"all",201708011200,"all"),
                           htower=c(28.2,33.2,3.53,6.55,9,5.82,6.08,3.2,5.2))
  
  ## US-RO2 heights
  zm.list.new3<-data.frame(site=rep("US-Ro2",11),
                           start=c("all",200807181200,200810011200,200907281200,200908310000,201001010000,
                                   201107211200,201110281200,201307151200,201310151200,201510301200),
                           end=c(200807181200,200810011200,200907281200,200908310000,201001010000,201107211200,
                                 201110281200,201307151200,201310151200,201510301200,201701010000),
                           htower=c(1.93,2.79,1.93,2.79,1.93,1.75,
                                    2.42,1.75,2.79,1.75,2.00))
  
  ## add 2013-2015 heights
  zm.list.new4<-data.frame(site=rep("US-Ro1",5),
                           start=c(201301010000,201507091200,201507211200,201508051200,201511041200),
                           end=c(201507091200,201507211200,201508051200,201511041200,201701010000),
                           htower=c(2.05,2.95,3.48,3.75,2.05))
  
  zm.list.new5<-data.frame(site=rep("US-Ro5",7),
                           start=c("all",201606241200,201608041200,201610191200,201807091200,201807181200,201811021200),
                           end=c(201606241200,201608041200,201610191200,201807091200,201807181200,201811021200,201901010000),
                           htower=c(2.06,2.2,3.85,2.45,3.25,3.8,2.78))
  
  zm.list.new6<-data.frame(site=rep("US-Br1",5),
                           start=c("all",201004220000,201010020000,201105040000,201110070000),
                           end=c(201004220000,201010020000,201105040000,201110070000,"all"),
                           htower=c(2.25,2.25,2.25,5.20,2.25))
  
  zm.list.new7<-data.frame(site=rep("US-Br3",5),
                           start=c("all",201004220000,201010020000,201105040000,201110070000),
                           end=c(201004220000,201010020000,201105040000,201110070000,"all"),
                           htower=c(2.00,5.00,2.00,2.00,2.00))
  
  zm.list.new8<-data.frame(site=rep("CA-ER1",10),
                           start=c(201501010000,201506221200,201507211200,201511021200,201601010000,201607111200,201610071200,201706291200,201707251200,201711101200),
                           end=c(201506221200,201507211200,201601010000,201607111200,201511021200,201610071200,201706291200,201707251200,201711101200,"all"),
                           htower=c(1.93,2.42,3.66,2.03,1.93,3.67,2,2.85,4.5,2.25))
  
  zm.list.new9<-data.frame(site=rep("US-Rpf",3),
                           start=c(200801010000,201307010000,201809010000),
                           end=c(201307010000,201809010000,"all"),
                           htower=c(2.6,3.9,5.6134))
  ## overwrite sites by above inputs
  replace.ls<-which(zm.list$site=="US-Ro2"|zm.list$site=="US-Ro5"|
                      zm.list$site=="US-NC3"|zm.list$site=="US-NC4")
  
  if(length(replace.ls)>0) zm.list<-zm.list[-replace.ls,]
  
  zm.list<-rbind.data.frame(zm.list,
                            zm.list.new,
                            zm.list.new2,
                            zm.list.new3,
                            zm.list.new4,
                            zm.list.new5,
                            zm.list.new6,
                            zm.list.new7,
                            zm.list.new8,
                            zm.list.new9)
  
  #target.site<-"AR-TF1"
  if(target.site %in% zm.list$site){
    
    
  }else{
    
    ## get EC height from VI
    var.info.grab<-var.info[var.info$Site_ID==target.site,]
    var.info.grab<-var.info.grab[!is.na(var.info.grab$Height),]
    
    if(nrow(var.info.grab)>0){
      ## decode the variable name, info of aggregation level, gap-filling,...
      basename_decode<-basename_pharse(var.name=var.info.grab$Variable,
                                       FP.ls=FP.ls,echo=F)
      
      ec.loc<-as.numeric(rownames(basename_decode)[which(basename_decode$basename=="USTAR"|
                                                           basename_decode$basename=="H"|
                                                           basename_decode$basename=="LE"|
                                                           basename_decode$basename=="FC")])
      if(length(ec.loc)==1){
        
        zm.list<-data.frame(zm.list,c(target.site),
                            "all",
                            "all",
                            var.info.grab$Height[ec.loc])
        
      }else if(length(ec.loc)>1&sd(var.info.grab$Height[ec.loc])==0){
        
        zm.list<-rbind(zm.list,c(target.site),
                       "all",
                       "all",
                       var.info.grab$Height[ec.loc][1])
        
      }else if(length(ec.loc)==0){
        
        print(paste("[warning] cannot find EC height"))
        
      }else{
        
        print(paste("[warning] may have multiple EC height"))
        print(var.info.grab[ec.loc,c("Variable","Start_Date","Height")])
        
      }
      
    }else{
      
      print(paste("[warning] cannot find valid VI"))
      
    }
  }
  zm.list$htower<-as.numeric(zm.list$htower)  
  
  zm.list.out<-zm.list[zm.list$site==target.site,]
  
  return(zm.list.out)
  
}





