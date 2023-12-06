rm(list = ls()) # clean working memory

require(ggplot2)
require(hexbin)
require(lmodel2)
require(MASS)
require(reshape2)
library(multcompView)
require(vioplot)

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\flux_timeseries_cluster\\"
path.in.root <- paste0(path, "cluster-output\\")
path.out.root <- paste0(path, "summary\\")
RDir <- paste0(path, "R\\")

out.version <- "20231019"

ver.ls <- c("20231019-8h5d", "20231019-12h7d", "20231019-24h15d")
ver.name.ls <- c("5 day - 3 hour",
                 "7 day - 2 hour",
                 "15 day - 1 hour")
dist.ls <- c("dtw")#, "euclidean")
target.var.ls <- c("NEE", "LE", "H", "USTAR", "NETRAD", "TA", "VPD", "SWC", "ALL_FLUX")

full.ls1 <- read.csv(paste0(path.in.root, ver.ls[1], "\\",
                            "ALL_BASE_site_short_list3.csv"),
                     header = T)
full.ls2 <- read.csv(paste0(path.in.root, ver.ls[2], "\\",
                            "ALL_BASE_site_short_list3.csv"),
                     header = T)
full.ls3 <- read.csv(paste0(path.in.root, ver.ls[3], "\\",
                            "ALL_BASE_site_short_list3.csv"),
                     header = T)

################################################################################
## dist by groups
full.ls2$eco_L1_name_short <- NA
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 1] <- "ARC"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 2] <- "TRA"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 3] <- "TGA"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 4] <- "HDP"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 5] <- "NTF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 6] <- "NWF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 7] <- "MCF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 8] <- "ETF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 9] <- "GPL"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 10] <- "NAD"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 11] <- "MCA"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 12] <- "SSA"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 13] <- "TSE"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 14] <- "TDF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 15] <- "THF"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 16] <- "CAI"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 17] <- "NAE"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 18] <- "CAE"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 19] <- "SAE"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 20] <- "AOL"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 21] <- "EHL"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 22] <- "GCH"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 23] <- "PMP"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 24] <- "MPA"
full.ls2$eco_L1_name_short[full.ls2$eco_L1 == 25] <- "HAW"

pair <- NULL
igbp.pair <- NULL
eco.pair <- NULL
for(l in 1:(nrow(full.ls2)-1)) {
  pair <- c(pair,
            paste0(full.ls2$SITE_ID[l],
                   "-",
                   full.ls2$SITE_ID[(l+1):nrow(full.ls2)]))
  igbp.pair <- c(igbp.pair,
                 paste0(full.ls2$IGBP[l],
                        "-",
                        full.ls2$IGBP[(l+1):nrow(full.ls2)]))
  eco.pair <- c(eco.pair,
                paste0(full.ls2$eco_L1_name_short[l],
                       "-",
                       full.ls2$eco_L1_name_short[(l+1):nrow(full.ls2)]))
}
dist.pile <- data.frame(pair,
                        igbp.pair,
                        eco.pair)

for(ii in 1:8){
  dist.tmp <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                                     "AMF-diurnal-seasonal-cluster-", 
                                     target.var.ls[ii], 
                                     "-distMatrix-",
                                     dist.ls,
                                     ".csv"),
                              header = T)
  dist.tmp <- dist.tmp[, -1]  
  colnames(dist.tmp) <- sub("\\.", "-", colnames(dist.tmp))
  for(i in 1:nrow(dist.tmp)){
    dist.tmp[i, i] <- NA
  }
  dist.tmp <-
    (dist.tmp - min(dist.tmp, na.rm = T)) / (max(dist.tmp, na.rm = T) - min(dist.tmp, na.rm = T))
  
  dist.tmp.flat <- data.frame(NULL)
    for(ll in 1:ncol(dist.tmp)) {
      dist.tmp.flat <- rbind.data.frame(dist.tmp.flat,
                                        data.frame(
                                          pair = paste0(colnames(dist.tmp)[ll],
                                                        "-",
                                                        colnames(dist.tmp[,ll:nrow(dist.tmp)])),
                                          dist = dist.tmp[ll:nrow(dist.tmp), ll]
                                        ))
  }
  colnames(dist.tmp.flat)[2] <- paste0(target.var.ls[ii])
  dist.tmp.flat<-na.omit(dist.tmp.flat)

  dist.pile <- merge.data.frame(dist.pile,
                                dist.tmp.flat,
                                all.x = T,
                                by = "pair")
  
}

dist.pile$igbp.group <- "mixed"
for(i in 1: nrow(dist.pile)){
  if(substr(dist.pile$igbp.pair[i], 1, 3) == substr(dist.pile$igbp.pair[i], 5, 7)){
    dist.pile$igbp.group[i] <- substr(dist.pile$igbp.pair[i], 1, 3)
  }
}
dist.pile$igbp.group[dist.pile$igbp.group %in% c("CSH", "OSH")] <- "SHB"

dist.pile$eco.group <- "mixed"
for(i in 1: nrow(dist.pile)){
  if(substr(dist.pile$eco.pair[i], 1, 3) == substr(dist.pile$eco.pair[i], 5, 7)){
    dist.pile$eco.group[i] <- substr(dist.pile$eco.pair[i], 1, 3)
  }
}

##############################################################################
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable) {
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][, 4]
  Tukey.labels <-
    data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment) ,]
  return(Tukey.labels)
}

######## work on Tukey Test
eco.order <- c("mixed", "NTF", "NWF", "ETF", "GPL", "NAD", "MCA")
igbp.order <- c("mixed", "ENF", "DBF", "CRO", "GRA", "SHB", "WET")

dist.pile.igbp <- dist.pile[dist.pile$igbp.group %in% igbp.order, ]
dist.pile.igbp$igbp.group <- factor(dist.pile.igbp$igbp.group,
                                    levels = igbp.order)
igbp.lab <- c("ALL", "ENF", "DBF", "CRO", "GRA", "SHB", "WET")

dist.pile.eco <- dist.pile[dist.pile$eco.group %in% eco.order,]
dist.pile.eco$eco.group <- factor(dist.pile.eco$eco.group,
                                    levels = eco.order)
eco.lab <- c("ALL", "NTF", "NWF", "ETF", "GPL", "NAD", "MCA")

######## Flux DTW
jpeg(paste0(path.out.root, out.version, "\\DTW_flux_byGroup.jpg"),
     width = 6.5, height = 4.5, units = "in", res = 300, pointsize = 11)
par(mfrow = c(2, 4), mar = c(2.5, 0, 0.5, 0), oma = c(0.5, 5.5, 0.5, 0.5))
for(vv in 1:4){
  igbp.var <- lm(dist.pile.igbp[, target.var.ls[vv]] ~ dist.pile.igbp$igbp.group)
  igbp.var.aov <- aov(igbp.var)
  
  # Tukey test to study each pair of treatment :
  igbp.var.tukey <- TukeyHSD(x = igbp.var.aov, 'dist.pile.igbp$igbp.group', conf.level = 0.95)
  print(igbp.var.tukey[[1]][c(1:6),])
  #plot(igbp.var.tukey, las = 1)
  
  # Apply the function on my dataset
  igbp.var.label <- generate_label_df(igbp.var.tukey , 'dist.pile.igbp$igbp.group')
  
  # Draw the basic boxplot
  a <-
    boxplot(
      dist.pile.igbp[, target.var.ls[vv]] ~ dist.pile.igbp$igbp.group,
      ylim = c(0, 1),
      col = "gray",
      border = "gray",
      ylab = "",
      xlab = "",
      main = "",
      las = 1,
      xaxt = "n",
      yaxt = "n",
      pch = 1,
      lty = 1,
      boxlty = 0,
      outcex = 0.2,
      medcol = "black",
      outcol = "lightgray",
      whisklwd = 2
    )
  points(tapply(dist.pile.igbp[, target.var.ls[vv]], 
                dist.pile.igbp$igbp.group, 
                function(x) mean(x, na.rm = T)),
         pch = 16)
  text(0.5, 1, paste0("(", letters[vv], ") ", target.var.ls[vv]), adj = c(0, 1), font = 2)
  text(c(1:length(igbp.lab))- 0.5, par("usr")[3]-0.02, labels = igbp.lab, srt = 60, pos = 1, xpd = TRUE)
  
  # I want to write the letter over each box. Over is how high I want to write it.
  over <- 0.1*max( a$stats[nrow(a$stats),] )
  #over <- 0.1*max( a$stats[nrow(a$stats),] )
  
  #Add the labels
  text(c(1:length(unique(dist.pile.igbp$igbp.group))),
       a$stats[nrow(a$stats), ] + over,
       igbp.var.label[igbp.order, 1],
       cex = 0.8)
  if(vv == 1){
    axis(2, 
         seq(0, 1, 0.25), 
         paste(c("(Similar)", "", "", "", "(Unique)"),
               seq(0, 1, 0.25)),
         las = 1)
  }
}

for(vv in 1:4){
  eco.var <- lm(dist.pile.eco[, target.var.ls[vv]] ~ dist.pile.eco$eco.group)
  eco.var.aov <- aov(eco.var)
  
  # Tukey test to study each pair of treatment :
  eco.var.tukey <- TukeyHSD(x = eco.var.aov, 'dist.pile.eco$eco.group', conf.level = 0.95)
  print(eco.var.tukey[[1]][c(1:6),])
  #plot(eco.var.tukey,
  #     las = 1)
  
  # Apply the function on my dataset
  eco.var.label <- generate_label_df(eco.var.tukey , 'dist.pile.eco$eco.group')
  
  # Draw the basic boxplot
  a <-
    boxplot(
      dist.pile.eco[, target.var.ls[vv]] ~ dist.pile.eco$eco.group,
      ylim = c(0, 1),
      col = "gray",
      border = "gray",
      ylab = "",
      xlab = "",
      main = "",
      las = 1,
      xaxt = "n",
      yaxt = "n",
      pch = 1,
      lty = 1,
      boxlty = 0,
      outcex = 0.2,
      medcol = "black",
      outcol = "lightgray",
      whisklwd = 2
    )
  points(tapply(dist.pile.eco[, target.var.ls[vv]], 
                dist.pile.eco$eco.group, 
                function(x) mean(x, na.rm = T)),
         pch = 16)
  text(0.5, 1, paste0("(", letters[vv + 4], ") ", target.var.ls[vv]), adj = c(0, 1), font = 2)
  text(c(1:length(eco.lab)) - 0.5, par("usr")[3]-0.02, labels = eco.lab, srt = 60, pos = 1, xpd = TRUE)
  
  # I want to write the letter over each box. Over is how high I want to write it.
  over <- 0.1*max( a$stats[nrow(a$stats),] )
  
  #Add the labels
  text(c(1:length(unique(dist.pile.eco$eco.group))),
       a$stats[nrow(a$stats), ] + over,
       eco.var.label[eco.order, 1],
       cex = 0.8)
  if(vv == 1){
    axis(2, 
         seq(0, 1, 0.25), 
         paste(c("(Similar)", "", "", "", "(Unique)"),
               seq(0, 1, 0.25)),
         las = 1)
  }
}
dev.off()

###### Met DTW
jpeg(paste0(path.out.root, out.version, "\\DTW_met_byGroup.jpg"),
     width = 6.5, height = 4.5, units = "in", res = 300, pointsize = 11)
par(mfrow = c(2, 4), mar = c(2.5, 0, 0.5, 0), oma = c(0.5, 5.5, 0.5, 0.5))
for(vv in 5:8){
  igbp.var <- lm(dist.pile.igbp[, target.var.ls[vv]] ~ dist.pile.igbp$igbp.group)
  igbp.var.aov <- aov(igbp.var)
  
  # Tukey test to study each pair of treatment :
  igbp.var.tukey <- TukeyHSD(x = igbp.var.aov, 'dist.pile.igbp$igbp.group', conf.level = 0.95)
  print(igbp.var.tukey[[1]][c(1:6),])
  #plot(igbp.var.tukey,
  #     las = 1)
  
  # Apply the function on my dataset
  igbp.var.label <- generate_label_df(igbp.var.tukey , 'dist.pile.igbp$igbp.group')
  
  # Draw the basic boxplot
  a <-
    boxplot(
      dist.pile.igbp[, target.var.ls[vv]] ~ dist.pile.igbp$igbp.group,
      ylim = c(0, 1),
      col = "gray",
      border = "gray",
      ylab = "",
      xlab = "",
      main = "",
      las = 1,
      xaxt = "n",
      yaxt = "n",
      pch = 1,
      lty = 1,
      boxlty = 0,
      outcex = 0.2,
      medcol = "black",
      outcol = "lightgray",
      whisklwd = 2
    )
  points(tapply(dist.pile.igbp[, target.var.ls[vv]], 
                dist.pile.igbp$igbp.group, 
                function(x) mean(x, na.rm = T)),
         pch = 16)
  text(0.5, 1, paste0("(", letters[vv-4], ") ", target.var.ls[vv]), adj = c(0, 1), font = 2)
  text(c(1:length(igbp.lab))- 0.5, par("usr")[3]-0.02, labels = igbp.lab, srt = 60, pos = 1, xpd = TRUE)
  
  # I want to write the letter over each box. Over is how high I want to write it.
  over <- 0.1*max( a$stats[nrow(a$stats),] )
  
  #Add the labels
  text(c(1:length(unique(dist.pile.igbp$igbp.group))),
       a$stats[nrow(a$stats), ] + over,
       igbp.var.label[igbp.order, 1],
       cex = 0.8)
  if(vv == 5){
    axis(2, 
         seq(0, 1, 0.25), 
         paste(c("(Similar)", "", "", "", "(Unique)"),
               seq(0, 1, 0.25)),
         las = 1)
  }
}

for(vv in 5:8){
  eco.var <- lm(dist.pile.eco[, target.var.ls[vv]] ~ dist.pile.eco$eco.group)
  eco.var.aov <- aov(eco.var)
  
  # Tukey test to study each pair of treatment :
  eco.var.tukey <- TukeyHSD(x = eco.var.aov, 'dist.pile.eco$eco.group', conf.level = 0.95)
  print(eco.var.tukey[[1]][c(1:6),])
  #plot(eco.var.tukey,
  #     las = 1)
  
  # Apply the function on my dataset
  eco.var.label <- generate_label_df(eco.var.tukey , 'dist.pile.eco$eco.group')
  
  # Draw the basic boxplot
  a <-
    boxplot(
      dist.pile.eco[, target.var.ls[vv]] ~ dist.pile.eco$eco.group,
      ylim = c(0, 1),
      col = "gray",
      border = "gray",
      ylab = "",
      xlab = "",
      main = "",
      las = 1,
      xaxt = "n",
      yaxt = "n",
      pch = 1,
      lty = 1,
      boxlty = 0,
      outcex = 0.2,
      medcol = "black",
      outcol = "lightgray",
      whisklwd = 2
    )
  points(tapply(dist.pile.eco[, target.var.ls[vv]], 
                dist.pile.eco$eco.group, 
                function(x) mean(x, na.rm = T)),
         pch = 16)
  text(0.5, 1, paste0("(", letters[vv], ") ", target.var.ls[vv]), adj = c(0, 1), font = 2)
  text(c(1:length(eco.lab)) - 0.5, par("usr")[3]-0.02, labels = eco.lab, srt = 60, pos = 1, xpd = TRUE)
  
  # I want to write the letter over each box. Over is how high I want to write it.
  over <- 0.1*max( a$stats[nrow(a$stats),] )
  
  #Add the labels
  text(c(1:length(unique(dist.pile.eco$eco.group))),
       a$stats[nrow(a$stats), ] + over,
       eco.var.label[eco.order, 1],
       cex = 0.8)
  if(vv == 5){
    axis(2, 
         seq(0, 1, 0.25), 
         paste(c("(Similar)", "", "", "", "(Unique)"),
               seq(0, 1, 0.25)),
         las = 1)
  }
}
dev.off()

################################################################################
###### Cluster comparison
clust.sum <- data.frame(match = c(paste(ver.name.ls[1], "vs", ver.name.ls[2]),
                                  paste(ver.name.ls[2], "vs", ver.name.ls[3]),
                                  "total site"))
for(l1 in 1:length(target.var.ls)){
  clust1 <- full.ls1[, paste(target.var.ls[l1], 
                             "clust_group",
                             dist.ls[1],
                             sep = "_")]
  clust2 <- full.ls2[, paste(target.var.ls[l1], 
                             "clust_group",
                             dist.ls[1],
                             sep = "_")]
  clust3 <- full.ls3[, paste(target.var.ls[l1], 
                             "clust_group",
                             dist.ls[1],
                             sep = "_")]
  
  clust.comp12 <- data.frame(table(clust1, clust2))
  clust.comp12 <- clust.comp12[clust.comp12$Freq != 0,]
  clust.comp12$match <- FALSE
  
  clust.comp12$match[
    clust.comp12$clust1 %in% names(table(clust.comp12$clust1))[table(clust.comp12$clust1) == 1] &
      clust.comp12$clust2 %in% names(table(clust.comp12$clust2))[table(clust.comp12$clust2) == 1]] <- TRUE
  for(m1 in 1:nrow(clust.comp12)){
    if(!clust.comp12$match[m1]){
      clust12.seach <- unique(c(which(clust.comp12$clust1 == clust.comp12$clust1[m1]),
                                which(clust.comp12$clust2 == clust.comp12$clust2[m1])))
      
      clust.comp12$match[clust12.seach[which(clust.comp12$Freq[clust12.seach] == max(clust.comp12$Freq[clust12.seach]))[1]]] <- TRUE
    }
  }
  
  clust.comp23 <- data.frame(table(clust2, clust3))
  clust.comp23 <- clust.comp23[clust.comp23$Freq != 0,]
  clust.comp23$match <- FALSE
  
  clust.comp23$match[
    clust.comp23$clust2 %in% names(table(clust.comp23$clust2))[table(clust.comp23$clust2) == 1] &
      clust.comp23$clust3 %in% names(table(clust.comp23$clust3))[table(clust.comp23$clust3) == 1]] <- TRUE
  for(m2 in 1:nrow(clust.comp23)){
    if(!clust.comp23$match[m2]){
      clust23.seach <- unique(c(which(clust.comp23$clust2 == clust.comp23$clust2[m2]),
                                which(clust.comp23$clust3 == clust.comp23$clust3[m2])))
      
      clust.comp23$match[clust23.seach[which(clust.comp23$Freq[clust23.seach] == max(clust.comp23$Freq[clust23.seach]))[1]]] <- TRUE
    }
  }
  clust.sum$tmp <- c(sum(clust.comp12$Freq[clust.comp12$match]),
                     sum(clust.comp23$Freq[clust.comp23$match]),
                     min(c(sum(clust.comp12$Freq), sum(clust.comp23$Freq))))
  colnames(clust.sum)[which(colnames(clust.sum) == "tmp")] <- target.var.ls[l1]
}
write.csv(clust.sum,
          paste0(path.out.root, out.version,
                 "\\Cluser_match_summary.csv"),
          row.names = F,
          quote = F)

################################################################################
###### Distance matrix comparison  
dist.sum <- data.frame(match = c(paste(ver.name.ls[1], "vs", ver.name.ls[2]),
                                  paste(ver.name.ls[2], "vs", ver.name.ls[3])))

fig.dist12.ls <- list(NULL)
fig.dist23.ls <- list(NULL)
for(l1 in 1:length(target.var.ls)){
  if(!target.var.ls[l1] %in% c("ALL_MET", "ALL_FLUX", "ALL_FLUXMET")){
    dist1 <- read.csv(paste0(path.in.root, ver.ls[1], "\\", 
                             "AMF-diurnal-seasonal-cluster-", 
                             target.var.ls[l1], 
                             "-distMatrix-",
                             dist.ls[1],
                             ".csv"),
                      header = T)
    dist2 <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                             "AMF-diurnal-seasonal-cluster-", 
                             target.var.ls[l1], 
                             "-distMatrix-",
                             dist.ls[1],
                             ".csv"),
                      header = T)
    dist3 <- read.csv(paste0(path.in.root, ver.ls[3], "\\", 
                             "AMF-diurnal-seasonal-cluster-", 
                             target.var.ls[l1], 
                             "-distMatrix-",
                             dist.ls[1],
                             ".csv"),
                      header = T)
    
    
    colnames(dist1)[-1] <- sub("\\.", "-", colnames(dist1)[-1])
    colnames(dist2)[-1] <- sub("\\.", "-", colnames(dist2)[-1])
    colnames(dist3)[-1] <- sub("\\.", "-", colnames(dist3)[-1])
    
    dist1.flat <- data.frame(NULL) 
    for(l2 in 2:ncol(dist1)) {
      dist1.flat <- rbind.data.frame(dist1.flat,
                                     data.frame(pair = paste0(colnames(dist1)[l2],
                                                              "-",
                                                              dist1[l2:nrow(dist1), 1]),
                                                dist = dist1[l2:nrow(dist1), l2]))
    }
    dist2.flat <- data.frame(NULL)
    for (l2 in 2:ncol(dist2)) {
      dist2.flat <- rbind.data.frame(dist2.flat,
                                     data.frame(pair = paste0(colnames(dist2)[l2],
                                                              "-",
                                                              dist2[l2:nrow(dist2), 1]),
                                                dist = dist2[l2:nrow(dist2), l2]))
    }
    dist3.flat <- data.frame(NULL)
    for (l2 in 2:ncol(dist3)) {
      dist3.flat <- rbind.data.frame(dist3.flat,
                                     data.frame(pair = paste0(colnames(dist3)[l2],
                                                              "-",
                                                              dist3[l2:nrow(dist3), 1]),
                                                dist = dist3[l2:nrow(dist3), l2]))
    }
    
    dist.merge <- merge.data.frame(
      merge.data.frame(
        dist1.flat,
        dist2.flat,
        by = "pair",
        suffixes = c("1", "2")
      ),
      dist3.flat,
      by = "pair"
    )
    colnames(dist.merge)[ncol(dist.merge)] <- "dist3"
    
    dist.merge <- na.omit(dist.merge)
    
    ### normalize
    for(dd in 2:4){
      dist.merge[, dd] <- (dist.merge[, dd] - min(dist.merge[, dd])) /
                             (max(dist.merge[, dd]) - min(dist.merge[, dd]))
    }
    
    
    ### prepare plot dist 1 vs dist 2 
    dist.sum12 <- lm(dist1 ~ dist2, data = dist.merge)
    summary(dist.sum12)
    mean(abs(dist.sum12$residuals))
    
    # get the kde2d information:
    dist.sum12.kde <-
      kde2d(dist.merge$dist2, dist.merge$dist1, n = 400)
    dx <-
      diff(dist.sum12.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
    dy <- diff(dist.sum12.kde$y[1:2])
    sz <- sort(dist.sum12.kde$z)
    c1 <- cumsum(sz) * dx * dy
    
    # specify desired contour levels:
    prob <- c(0.9, 0.8, 0.7, 0.6, 0.5)
    
    # plot:
    dimnames(dist.sum12.kde$z) <-
      list(dist.sum12.kde$x, dist.sum12.kde$y)
    dc <- reshape2::melt(dist.sum12.kde$z)
    dc$prob <- approx(sz, 1 - c1, dc$value)$y
    
    fig.dist12 <- ggplot(dc,
                         aes(x = Var1, y = Var2)) +
      theme_grey(base_size = 11) +
      labs(x = ver.name.ls[2],
           y = ver.name.ls[1],
           tag = paste0("(",
                        ifelse(l1 < 5, letters[2 * l1 - 1], letters[2 * (l1 -
                                                                           4) - 1]) ,
                        ")")) +
      xlim(0, max(dist.merge[, 2:4]) * 0.9) +
      ylim(0, max(dist.merge[, 2:4]) * 0.9) +
      geom_point(
        data = dist.merge,
        aes(x = dist2, y = dist1),
        shape = '.',
        alpha = 0.6,
        col = "grey30"
      ) +
      geom_contour(aes(z = prob),
                   breaks = prob,
                   col = "yellow",
                   lwd = 0.6) +
      geom_abline(
        slope = 1,
        intercept = 0,
        show.legend = "1:1 reference line",
        col = "black",
        lwd = 0.5
      ) +
      geom_abline(
        slope = dist.sum12$coefficients[2],
        intercept = dist.sum12$coefficients[1],
        show.legend = "regression",
        col = "red",
        lwd = 0.8
      )
    
    ### prepare plot dist 2 vs dist 3
    dist.sum23 <- lm(dist3 ~ dist2, data = dist.merge)
    summary(dist.sum23)
    mean(abs(dist.sum23$residuals))
    
    # get the kde2d information:
    dist.sum23.kde <-
      kde2d(dist.merge$dist2, dist.merge$dist3, n = 400)
    dx <-
      diff(dist.sum23.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
    dy <- diff(dist.sum23.kde$y[1:2])
    sz <- sort(dist.sum23.kde$z)
    c1 <- cumsum(sz) * dx * dy
    
    # specify desired contour levels:
    prob <- c(0.9, 0.8, 0.7, 0.6, 0.5)
    
    # plot: 
    dimnames(dist.sum23.kde$z) <-
      list(dist.sum23.kde$x, dist.sum23.kde$y)
    dc <- reshape2::melt(dist.sum23.kde$z)
    dc$prob <- approx(sz, 1 - c1, dc$value)$y
    
    fig.dist23 <- ggplot(dc,
                         aes(x = Var1, y = Var2)) +
      theme_grey(base_size = 11) +
      labs(x = ver.name.ls[2],
           y = ver.name.ls[3],
           tag = paste0("(",
                        ifelse(l1 < 5, letters[2 * l1], letters[2 * (l1 - 4)]) ,
                        ")")) +
      xlim(0, max(dist.merge[, 2:4]) * 0.9) +
      ylim(0, max(dist.merge[, 2:4]) * 0.9) +
      geom_point(
        data = dist.merge,
        aes(x = dist2, y = dist3),
        shape = '.',
        alpha = 0.6,
        col = "grey30"
      ) +
      geom_contour(aes(z = prob),
                   breaks = prob,
                   col = "yellow",
                   lwd = 0.6) +
      geom_abline(
        slope = 1,
        intercept = 0,
        show.legend = "1:1 reference line",
        col = "black",
        lwd = 0.5
      ) +
      geom_abline(
        slope = dist.sum23$coefficients[2],
        intercept = dist.sum23$coefficients[1],
        show.legend = "regression",
        col = "red",
        lwd = 0.8
      )
    
    fig.dist12.ls[[l1]] <- fig.dist12
    fig.dist23.ls[[l1]] <- fig.dist23
    
    dist.sum$R2 <- c(round(summary(dist.sum12)$r.squared, digits = 2),
                     round(summary(dist.sum23)$r.squared, digits = 2))
    dist.sum$intercept <- c(round(dist.sum12$coefficients[1], digits = 2),
                            round(dist.sum23$coefficients[1], digits = 2))
    dist.sum$slope <- c(round(dist.sum12$coefficients[2], digits = 2),
                            round(dist.sum23$coefficients[2], digits = 2))
    dist.sum$rmse <- c(round(sqrt(sum(dist.sum12$residuals^2/dist.sum12$df.residual)), digits = 0),
                       round(sqrt(sum(dist.sum23$residuals^2/dist.sum23$df.residual)), digits = 0))
    colnames(dist.sum)[c((ncol(dist.sum) - 3):ncol(dist.sum))] <-
      paste(target.var.ls[l1], c("R2", "intercept", "slope", "rmse"))
  }
}
write.csv(dist.sum,
          paste0(path.out.root, out.version,
                 "\\Dist_matrix_", dist.ls[1], "summary.csv"),
          row.names = F,
          quote = F)

#### Distance matrix figure Flux
png(
  paste0(
    path.out.root,
    out.version,
    "\\Dist_matrix_", 
    dist.ls[1],
    "_window_compare_flux.png"
  ),
  height = 11,
  width = 6,
  res = 300,
  units = "in",
  pointsize = 11
)
print(gridExtra::grid.arrange(fig.dist12.ls[[1]],
                              fig.dist23.ls[[1]],
                              fig.dist12.ls[[2]],
                              fig.dist23.ls[[2]],
                              fig.dist12.ls[[3]],
                              fig.dist23.ls[[3]],
                              fig.dist12.ls[[4]],
                              fig.dist23.ls[[4]],
                              ncol = 2))
dev.off()

#### Distance matrix figure Met
png(
  paste0(
    path.out.root,
    out.version,
    "\\Dist_matrix_", 
    dist.ls[1],
    "_window_compare_met.png"
  ),
  height = 11,
  width = 6,
  res = 300,
  units = "in",
  pointsize = 11
)
print(gridExtra::grid.arrange(fig.dist12.ls[[5]],
                              fig.dist23.ls[[5]],
                              fig.dist12.ls[[6]],
                              fig.dist23.ls[[6]],
                              fig.dist12.ls[[7]],
                              fig.dist23.ls[[7]],
                              fig.dist12.ls[[8]],
                              fig.dist23.ls[[8]],
                              ncol = 2))
dev.off()

################################################################################
### harmonic Map
dist1 <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                         "AMF-diurnal-seasonal-cluster-", 
                         target.var.ls[1], 
                         "-distMatrix-",
                         dist.ls,
                         ".csv"),
                  header = T)
dist2 <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                         "AMF-diurnal-seasonal-cluster-", 
                         target.var.ls[2], 
                         "-distMatrix-",
                         dist.ls,
                         ".csv"),
                  header = T)
dist3 <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                         "AMF-diurnal-seasonal-cluster-", 
                         target.var.ls[3], 
                         "-distMatrix-",
                         dist.ls,
                         ".csv"),
                  header = T)
dist4 <- read.csv(paste0(path.in.root, ver.ls[2], "\\", 
                         "AMF-diurnal-seasonal-cluster-", 
                         target.var.ls[4], 
                         "-distMatrix-",
                         dist.ls,
                         ".csv"),
                  header = T)
dist1 <- dist1[, -1]
dist2 <- dist2[, -1]
dist3 <- dist3[, -1]
dist4 <- dist4[, -1]

colnames(dist1) <- sub("\\.", "-", colnames(dist1))
colnames(dist2) <- sub("\\.", "-", colnames(dist2))
colnames(dist3) <- sub("\\.", "-", colnames(dist3))
colnames(dist4) <- sub("\\.", "-", colnames(dist4))

for(i in 1:nrow(dist1)){
  dist4[i, i] <- dist3[i, i] <- dist2[i, i] <- dist1[i, i] <- NA
}

all.site <- full.ls2$SITE_ID[!is.na(full.ls2$ALL_FLUX_clust_group_dtw)]

dist1.get <- which(colnames(dist1) %in% all.site)
dist1 <- dist1[dist1.get, dist1.get]
dist2.get <- which(colnames(dist2) %in% all.site)
dist2 <- dist2[dist2.get, dist2.get]
dist3.get <- which(colnames(dist3) %in% all.site)
dist3 <- dist3[dist3.get, dist3.get]
dist4.get <- which(colnames(dist4) %in% all.site)
dist4 <- dist4[dist4.get, dist4.get]

dist1 <-
  (dist1 - min(dist1, na.rm = T)) / (max(dist1, na.rm = T) - min(dist1, na.rm = T))
dist2 <-
  (dist2 - min(dist2, na.rm = T)) / (max(dist2, na.rm = T) - min(dist2, na.rm = T))
dist3 <-
  (dist3 - min(dist3, na.rm = T)) / (max(dist3, na.rm = T) - min(dist3, na.rm = T))
dist4 <-
  (dist4 - min(dist4, na.rm = T)) / (max(dist4, na.rm = T) - min(dist4, na.rm = T))

dist.all <- (dist1 + dist2 + dist3 + dist4) / 4
HUP <- data.frame(SITE_ID = colnames(dist.all),
                  HUP = apply(dist.all, 2, function(x) mean(x, na.rm = T)))

full.ls2 <- merge.data.frame(full.ls2,
                             HUP,
                             by = "SITE_ID",
                             all.x = T)

full.ls2$HUP.col <- NA
full.ls2$HUP.col[which(full.ls2$HUP < 0.18)] <- 1
full.ls2$HUP.col[which(full.ls2$HUP >= 0.18 & full.ls2$HUP < 0.20)] <- 2
full.ls2$HUP.col[which(full.ls2$HUP >= 0.20 & full.ls2$HUP < 0.24)] <- 3
full.ls2$HUP.col[which(full.ls2$HUP >= 0.24 & full.ls2$HUP < 0.29)] <- 4
full.ls2$HUP.col[which(full.ls2$HUP >= 0.29)] <- 5
#full.ls2$HUP.col[which(full.ls2$HUP >= 0.35)] <- 6

## Map 
png(
  paste0(path.out.root,
         out.version,
         "\\HUP-ALLFLUX-map.png"),
  height = 5,
  width = 4.6,
  res = 400,
  units = "in",
  pointsize = 10
)
par(
  fig = c(0, 1, 0, 1),
  oma = c(0, 0, 0, 0),
  mar = c(0, 0, 0, 0)
)
plot(
  0,
  0,
  type = "n",
  xlim = c(-170, -30),
  ylim = c(-60, 75),
  yaxt = "n",
  xaxt = "n",
  xlab = "",
  ylab = "",
  xaxs = "i",
  yaxs = "i",
  bty = "n"
)
maps::map(
  "world",
  fill = TRUE,
  border = "white",
  col = "grey80",
  bg = "white",
  add = T,
  xlim = c(-170, -30),
  ylim = c(-60, 75),
  mar = c(0, 0, 0, 0)
)
polygon(c(-60, -30, -30, -60),
        c(60, 60, 75, 75),
        border = "white",
        col= "white")
polygon(c(-70, -30, -30, -70),
        c(75, 75, 90, 90),
        border = "white",
        col= "white")
points(
  full.ls2$LOCATION_LONG[!is.na(full.ls2[,paste0("ALL_FLUX_clust_group_dtw")])],
  full.ls2$LOCATION_LAT[!is.na(full.ls2[,paste0("ALL_FLUX_clust_group_dtw")])],
  bg = rev(hcl.colors(5, palette = "Dark 3"))[full.ls2$HUP.col[!is.na(full.ls2[,paste0("ALL_FLUX_clust_group_dtw")])]],
  col = "white",
  cex = 0.75,
  pch = 21,
  lwd = 0.4
)
legend(
  -165,
  10,
  title = "Group",
  title.adj = 0,
  legend = c(#"> 0.35",
             "> 0.29",
             "0.24 - 0.29",
             "0.20 - 0.24",
             "0.18 - 0.20",
             "0.15 - 0.18"),
  fill = hcl.colors(5, palette = "Dark 3"),
  border = NA,
  bty = "n",
  cex = 0.75,
  ncol = 1
)
dev.off()

write.csv(full.ls2,
          paste0(path.out.root, out.version, 
                 "\\ALL_BASE_site_short_list4.csv"),
          quote = F,
          row.names = F)
