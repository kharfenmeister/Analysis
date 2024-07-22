### Smooth Sen1-Timeseries with LOESS and calculate minima, maxima and breakpoints for single fields ###
library(plyr)
library(strucchange)

# Data direction
data_pol_path <- "Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Data_Field_Mean_pol_spk/"
data_bs_path <- "Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Data_Field_Mean/"
bbch_path <- "Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/bbch/csv/"
bbch_dwd_path <- "Z:/Katharina/AgriFusion/Analyse/Phänologie/"

# Field Info
field_combis <- read.csv("Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/field_combis.csv", header=TRUE, sep=";", stringsAsFactors = FALSE)

bbch_stages <- c("31", "55", "75", "87", "99")

## Write all Sen1 combinations
crop <- c("wheat", "barley")
bs_param <- c("Alpha", "Entropy", "Anisotropy", "Sigma_VH", "Sigma_VV", "Sigma_VHVV")
bs_pass <- c("ASC_44", "ASC_146", "DESC_95", "DESC_168")

sets <- expand.grid(crop, bs_param, bs_pass)
names(sets) <- c("crop", "param", "pass")


for (s in 1:nrow(sets)){
  crop <- as.character(sets$crop[s])
  param <- as.character(sets$param[s])
  pass <- strsplit(as.character(sets$pass[s]), "_")[[1]][1]
  orbit <- strsplit(as.character(sets$pass[s]), "_")[[1]][2]
  
  if (param == "Alpha" | param == "Entropy" | param == "Anisotropy"){
    data_list <- list.files(data_pol_path, pattern=paste(param, pass, orbit, sep="_"), full.names=TRUE, recursive=TRUE)
  } else if (param == "Sigma_VH" | param == "Sigma_VV" | param == "Sigma_VHVV"){
    data_list <- list.files(data_bs_path, pattern=paste(param, pass, orbit, sep="_"), full.names=TRUE, recursive=TRUE)
  }
  
  field <- unique(unlist(lapply(strsplit(basename(data_list), "_"), "[[", 1)))
  field_crop <- field_combis$field[field_combis$crop == crop]
  field <- intersect(field, field_crop)
  
  for (f in 1:length(field)){
    fieldname <- field[f]
    year <- field_combis$year[field_combis$field == fieldname]
    site <- field_combis$site[field_combis$field == fieldname]
    data_field <- read.csv(data_list[grepl(paste(field[f], "_", param, sep=""), data_list)], header=TRUE)
    
    if (fieldname == "KarBio" & as.character(sets$pass[s]) == "ASC_44") next
    
    # BBCH
    bbch_data <- read.csv(list.files(bbch_path, pattern=paste(crop, "_", fieldname, ".csv", sep=""), full.names=TRUE), header=TRUE)
    bbch_data <- bbch_data[!is.na(bbch_data$bbch_median),]
    bbch_data$dates <- as.character(as.Date(bbch_data$dates, format="%d.%m.%Y"))
    bbch_data <- merge(bbch_data, data_field[,c(1,2)], by.x="dates", by.y="dates_all")
    
    #BBCH DWD
    bbch_dwd_annual <- read.csv(paste(bbch_dwd_path, "annual/", year, "_", crop, "_", site, "_pheno_annual.csv", sep=""), header=TRUE, sep=",")
    stations_col_annual <- grep("^X", names(bbch_dwd_annual))
    names(bbch_dwd_annual)[stations_col_annual] <- paste(names(bbch_dwd_annual)[stations_col_annual], "a", sep="_")
    bbch_dwd_annual <- bbch_dwd_annual[,1:(length(stations_col_annual)+1)]
    
    bbch_dwd_immediate <- read.csv(paste(bbch_dwd_path, year, "_", crop, "_", site, "_pheno.csv", sep=""), header=TRUE, sep=",")
    stations_col_immediate <- grep("^X", names(bbch_dwd_immediate))
    names(bbch_dwd_immediate)[stations_col_immediate] <- paste(names(bbch_dwd_immediate)[stations_col_immediate], "i", sep="_")
    bbch_dwd_immediate <- bbch_dwd_immediate[,1:(length(stations_col_immediate)+1)]
    
    bbch_dwd_all <-  merge(bbch_dwd_annual, bbch_dwd_immediate, by="bbch_stages", all.x=TRUE)
    stations_col <- grep("^X", names(bbch_dwd_all))
    
    bbch_dwd <-  as.data.frame(bbch_dwd_all$bbch_stages)
    names(bbch_dwd) <- c("bbch_stages")
    
    # 1. station
    bbch_dwd_all[,2] <- as.character(as.Date(bbch_dwd_all[,2], format="%d-%m-%Y"))
    bbch_dwd_x <- merge(bbch_dwd_all[,c(1,2)], data_field[,c(1,2)], by.x=colnames(bbch_dwd_all)[2], by.y="dates_all")
    bbch_dwd_x$station <- names(bbch_dwd_x[1])
    
    bbch_dwd <- merge(bbch_dwd, bbch_dwd_x, by="bbch_stages", all.x=TRUE)
    names(bbch_dwd)[2] <- "date"
    
    # Following stations
    for (x in 3:tail(stations_col, n=1)){
      bbch_dwd_all[,x] <- as.character(as.Date(bbch_dwd_all[,x], format="%d-%m-%Y"))
      bbch_dwd_x <- merge(bbch_dwd_all[,c(1,x)], data_field[,c(1,2)], by.x=colnames(bbch_dwd_all)[x], by.y="dates_all")
      bbch_dwd_x$station <- names(bbch_dwd_x[1])
      
      bbch_dwd_x <- bbch_dwd_x[,c(2,1,3,4)]
      names(bbch_dwd_x)[2] <- "date"
      bbch_dwd <- rbind(bbch_dwd, bbch_dwd_x)
    }
    
    bbch_dwd <- na.omit(bbch_dwd)
    
    bbch_dwd <- bbch_dwd[order(bbch_dwd$day),]
    
    # LOESS Smoothing
    loessMod <- loess(data_field[,4] ~ data_field$day, span=0.3, degree=1, na.action=na.exclude) 
    data_field$loess <- predict(loessMod) 
    
    # Local Minima and Maxima
    loess_data <- data_field$loess[!is.na(data_field$loess)]
    
    loc_min <- c()
    loc_max <- c()
    for (i in 2:(length(loess_data)-1)){
      if ((loess_data[i] >= loess_data[i+1]) & (loess_data[i] >= loess_data[i-1])){
        loc_max <- c(loc_max, loess_data[i])
      } else if ((loess_data[i] <= loess_data[i+1]) & (loess_data[i] <= loess_data[i-1])){
        loc_min <- c(loc_min, loess_data[i])
      }
    }
    
    loc_min_days <- c()
    for (i in 1:length(loc_min)){
      loc_min_days <- c(loc_min_days, data_field$day[which(data_field$loess == loc_min[i])])
    }
    
    loc_max_days <- c()
    for (i in 1:length(loc_max)){
      loc_max_days <- c(loc_max_days, data_field$day[which(data_field$loess == loc_max[i])])
    }
    
    loess_min <- as.data.frame(cbind(loc_min_days, loc_min))
    loess_max <- as.data.frame(cbind(loc_max_days, loc_max))
    
    # Breakpoints
    ts_days <- na.omit(cbind(data_field$day,data_field[,4]))[,1]
    ts_value <- na.omit(cbind(data_field$day,data_field[,c("loess")]))[,2]
    
    break_points <- breakpoints(ts_value~1)
    break_dates <- ts_days[break_points$breakpoints]
    
    ### Plot
    output <- paste("Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Plots_Mean_combis_LOESS_maxminbreak_bbch_all/", crop, "/", param, "_", pass, orbit, "_", fieldname, "_loess_bbch_all.png", sep="")
    windowsFonts(Calibri=windowsFont("Calibri"))
    png(output, width=3200, height=2200, res=300)
    par(mar = c(6, 5, 3, 3), yaxs="i", xpd=NA, family="Calibri")
    
    if (param == "Entropy" | param == "Anisotropy"){
      ylim <- c(0,1)
      bbch_text_pos <- 0.05
    } else if (param == "Alpha"){
      ylim <- c(-20,50)
      bbch_text_pos <- 1
    } else if (param == "Sigma_VH" | param == "Sigma_VV" | param == "Sigma_VHVV"){
      ylim <- c(-30,0)
      bbch_text_pos <- 1
    }
    
    plot(1, type="n", xlab="", ylim=ylim, xlim=c(0,153), xaxt="n", yaxt="n", ylab=paste(param, pass, orbit, sep=" "), cex.main=2, cex.axis=1.5, cex.lab=1.5)
    axis(2, cex.axis=1.5)
    axis(1, at=c(1, 32, 62, 93, 123), c("01.03.", "01.04.", "01.05.", "01.06.", "01.07."), cex.axis=1.5)
    abline(v=c(1, 32, 62, 93, 123), col="gray", xpd=FALSE)
    mtext(fieldname, 1, line=2.5, cex=2)
    
    # Harvest
    if (crop == "wheat" & year == "2018"){
      rect(139, ylim[1], 149, ylim[2], col=adjustcolor("gray40", alpha.f = 0.3), border=NA) # wheat 2018
    } else if (crop == "barley" & year == "2017"){
      rect(127, ylim[1], 139, ylim[2], col=adjustcolor("gray", alpha.f = 0.5), border=NA) # barley 2017
    } else if (crop == "barley" & year == "2018"){
      rect(119, ylim[1], 131, ylim[2], col=adjustcolor("gray40", alpha.f = 0.3), border=NA) # barley 2018
    } else if (crop == "wheat" & year == "2017"){
      rect(153, ylim[1], 163, ylim[2], col=adjustcolor("gray", alpha.f = 0.5), border=NA) # wheat 2017
    }
    
    # BBCH
    abline(v=bbch_data$day, col="black", lwd=2, xpd=FALSE)
    
    for (b in 1:nrow(bbch_data)){
      text(bbch_data$day[b], ylim[2]+bbch_text_pos, bbch_data$bbch_median[b])
    }
    
    # BBCH DWD
    bbch_dwd_cols <- c("green", "blue", "purple", "green", "blue")
    for (b in 1:length(bbch_stages)){
      bbch <- bbch_stages[b]
      abline(v=bbch_dwd$day[bbch_dwd$bbch_stages == bbch], col=bbch_dwd_cols[b], lwd=2, xpd=FALSE)
    }
    
    # Data
    lines(na.omit(cbind(data_field$day,data_field[,4])), type="b", col=field_combis$col[field_combis$field == field[f]], lwd=3, pch=16, xaxt="n", xlab="")
    lines(na.omit(cbind(data_field$day,data_field[,c("loess")])), type="b", col="green", lwd=3, pch=16, xaxt="n", xlab="")
    
    # Min / Max / Break
    points(loess_min, col="black", pch=25, cex=2, bg="black")
    points(loess_max, col="black", pch=24, cex=2, bg="black")
    
    abline(v=break_dates, col="red", lty=2, xpd=FALSE)
    
    dev.off()
    
  } #fields
} # sets
