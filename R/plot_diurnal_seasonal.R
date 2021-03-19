plot_diurnal_seasonal<-function(data.median,
                                data.upp1=NULL,  # optional, upper bound 1
                                data.upp2=NULL,  # optional, upper bound 2
                                data.low1=NULL,  # optional, lower bound 1
                                data.low2=NULL,  # optional, lower bound 2
                                path.out,        # output path
                                case.name,       # case name for putput files
                                d.hr,            # number of points per day in the diurnal-seasonal file  
                                n.wd,            # number of window per year in the diurnal-seasonal file  
                                l.wd,            # length of wondow (days)
                                skip.var.loc=c(1:2),  # skip plotting for such variables
                                n.data.yr=NA          # optional, number of years in original data
                                ){
  
  var.loc.ls<-c(1:ncol(data.median))
  for(i in 1:length(skip.var.loc)){
    var.loc.ls<-var.loc.ls[-which(var.loc.ls==skip.var.loc[i])]
  }
  
  png(filename=paste(path.out,case.name.short,"_diurnalseasonal-%d.png",sep=""),width=10,height=9,units="in",pointsize=10,res=300)
  par(mfrow=c(4,2),oma=c(2.5,0.5,2.5,0.5),mar=c(3.5,4.5,3.5,0.5))
  
  
  
  for(j2 in var.loc.ls){
    
    ## only plot not all empty variables
    if(sum(!is.na(data.median[,j2]))>0){
      
      rng<-range(c(data.upp2[,j2],data.upp1[,j2],data.median[,j2],data.low1[,j2],data.low2[,j2]),na.rm=T)
      
      plot(data.median[,j2],ylim=rng,ylab="",pch=16,cex=0.9,xaxt="n",xlab="")
      
      lines(data.upp1[,j2],col="darkgrey")
      lines(data.low1[,j2],col="darkgrey")
      points(data.upp1[,j2],pch=16,cex=0.5,col="darkgrey")
      points(data.low1[,j2],pch=16,cex=0.5,col="darkgrey")
      
      lines(data.upp2[,j2],col="grey")
      lines(data.low2[,j2],col="grey")
      points(data.upp2[,j2],pch=16,cex=0.5,col="grey")
      points(data.low2[,j2],pch=16,cex=0.5,col="grey")
      
      axis(side=3,at=c(0.5*d.hr+d.hr*c(0:(n.wd-1))),labels=c(0.5*l.wd+l.wd*c(0:(n.wd-1))))
      
      mtext(side=2,paste(colnames(data.median)[j2]),line=2.5,font=2)
      mtext(side=3,"DOY (central date)",line=0,outer=T)
      
      abline(v=d.hr*c(0:n.wd)+0.5,col="grey",lty=4)
      axis(side=1,at=seq(1,nrow(data.median),by=d.hr/4)-0.5,labels=rep(c(0,6,12,18),times=nrow(data.median)/d.hr),cex.axis=0.5)
      
      mtext(side=1,"HOUR",line=0,outer=T)
      
      if(!is.na(n.data.yr)) mtext(side=3,paste("N=",n.data.yr,"year"),adj=0,line=0,outer=T)
      
    }
  }
  dev.off()
  
}



