### Plot mean values of all barley / all wheat fields ###

crop <- "wheat"
year <- "2017"

field_info <- read.csv("Z:/Katharina/AgriFusion/Feldmessungen/fieldnames_years_sites_crop.csv", header=TRUE, sep=";")
fields <- as.character(field_info$field[field_info$year == year & field_info$crop == crop])

## Precipitation
prec_data_demmin <- read.csv(paste("Z:/Katharina/AgriFusion/Daten/Wetterdaten/precipitation_heydenhof_", year, "_daily.csv", sep=""), header=TRUE)
prec_data_demmin$dates <- as.Date(prec_data_demmin$dates)
prec_data_demmin$prec_daily[prec_data_demmin$prec_daily == 0]<- NA

prec_data_bloen <- read.csv(paste("Z:/Katharina/AgriFusion/Daten/Wetterdaten/precipitation_naundorf_", year, "_daily.csv", sep=""), header=TRUE)
prec_data_bloen$dates <- as.Date(prec_data_bloen$dates)
prec_data_bloen$prec_daily[prec_data_bloen$prec_daily == 0]<- NA

## Temperature
temp_demmin <- read.csv(paste("Z:/Katharina/AgriFusion/Daten/Wetterdaten/Temperaturen/temp_gdd_heydenhof_", year, "_daily.csv", sep=""), header=TRUE)
temp_demmin$dates <- as.Date(temp_demmin$dates)
temp_demmin <- temp_demmin[60:nrow(temp_demmin),]

temp_bloen <- read.csv(paste("Z:/Katharina/AgriFusion/Daten/Wetterdaten/Temperaturen/temp_gdd_langenlipsdorf_", year, "_daily.csv", sep=""), header=TRUE)
temp_bloen$dates <- as.Date(temp_bloen$dates)
temp_bloen<- temp_bloen[60:nrow(temp_bloen),]

sen1_dir <- "Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Data_Field_Mean_all_pol_Spk/" 

bs_param <- c("Alpha", "Entropy", "Anisotropy", "m1", "m2", "m3", "m4")
bs_pass <- c("ASC_44", "ASC_146", "DESC_95", "DESC_168")


for (p in 1:length(bs_param)){
  sen1_data <- list.files(paste(sen1_dir, year, "/", crop, "/", sep=""), pattern=paste("_", bs_param[p], "_", sep=""), full.names=TRUE)
  
  ## Plot Time Series ASC and DESC
  output <- paste("Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Plots_Field_Mean_all_pol_Spk/", year, "/", crop, "/", year, "_", crop, "_", bs_param[p], "_fields_all_mean_gdd_sm.png", sep="")
  windowsFonts(Calibri=windowsFont("Calibri"))
  png(output, width=3200, height=2200, res=300)
  par(mar = c(6, 5, 3, 9), yaxs="i", family="Calibri")
  labels <- c("ASC 44", "ASC 146", "DESC 95",  "DESC 168")
  ifelse(bs_param[p] == "Alpha", ylim <- c(-20, 50), ylim <- c(0, 1))
  
  # 1) ASC 44
  sen1_asc44 <- read.csv(sen1_data[grepl("ASC_44", sen1_data)], header=TRUE)
  
  plot(1, type="n", xlab="", ylim=ylim, xlim=c(0,153), xaxt="n", yaxt="n", ylab=paste(bs_param[p], sep=" "), cex.main=2, cex.axis=2, cex.lab=2)
  axis(2, cex.axis=1.5)
  
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
  
  # Soil Moisture
  if (crop == "wheat" & year == "2017"){
    soilmoisture <- read.csv("Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/soil_moisture/2017_wheat_soilmoisture.csv", header=TRUE)
  } else if (crop == "wheat" & year == "2018"){
    soilmoisture <- read.csv("Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/soil_moisture/2018_wheat_soilmoisture.csv", header=TRUE)
  } else if (crop == "barley" & year == "2017"){
    soilmoisture <- read.csv("Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/soil_moisture/2017_barley_soilmoisture.csv", header=TRUE)
  } else if (crop == "barley" & year == "2018"){
    soilmoisture <- read.csv("Z:/Katharina/AgriFusion/Analyse/Feldmessungen_Timeseries/soil_moisture/2018_barley_soilmoisture.csv", header=TRUE)
  }
  
  poly.x <- c(na.omit(cbind(sen1_asc44$day.x, sen1_asc44$all_mean+sen1_asc44$all_sd))[,1], rev(na.omit(cbind(sen1_asc44$day.x, sen1_asc44$all_mean+sen1_asc44$all_sd))[,1]))
  poly.y <- c(na.omit(cbind(sen1_asc44$day.x, sen1_asc44$all_mean+sen1_asc44$all_sd))[,2], rev(na.omit(cbind(sen1_asc44$day.x, sen1_asc44$all_mean-sen1_asc44$all_sd))[,2]))
  polygon(poly.x, poly.y, col=adjustcolor("blue", alpha.f = 0.2), border=NA)
  lines(na.omit(cbind(sen1_asc44$day.x, sen1_asc44$all_mean)), type="b", col="blue", lwd=3, pch=1, xaxt="n", xlab="")
  #axis(1, at=sen1_asc44$day.x[!is.na(sen1_asc44$all_mean)], sen1_asc44$label.x[!is.na(sen1_asc44$all_mean)], cex.axis=1.5, las=2)
  
  # 2) ASC 146
  sen1_asc146 <- read.csv(sen1_data[grepl("ASC_146", sen1_data)], header=TRUE)
  
  par(new=T)
  plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5) 
  poly.x <- c(na.omit(cbind(sen1_asc146$day.x, sen1_asc146$all_mean+sen1_asc146$all_sd))[,1], rev(na.omit(cbind(sen1_asc146$day.x, sen1_asc146$all_mean+sen1_asc146$all_sd))[,1]))
  poly.y <- c(na.omit(cbind(sen1_asc146$day.x, sen1_asc146$all_mean+sen1_asc146$all_sd))[,2], rev(na.omit(cbind(sen1_asc146$day.x, sen1_asc146$all_mean-sen1_asc146$all_sd))[,2]))
  polygon(poly.x, poly.y, col=adjustcolor("red", alpha.f = 0.2), border=NA)
  lines(na.omit(cbind(sen1_asc146$day.x, sen1_asc146$all_mean)), type="b", col="red", lwd=3, pch=1, xaxt="n", xlab="")
  #axis(1, at=sen1_asc146$day.x[!is.na(sen1_asc146$all_mean)], sen1_asc146$label.x[!is.na(sen1_asc146$all_mean)], cex.axis=1.5, las=2)
  
  # 3) DESC 95
  sen1_desc95 <- read.csv(sen1_data[grepl("DESC_95", sen1_data)], header=TRUE)
  
  par(new=T)
  plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5) 
  poly.x <- c(na.omit(cbind(sen1_desc95$day.x, sen1_desc95$all_mean+sen1_desc95$all_sd))[,1], rev(na.omit(cbind(sen1_desc95$day.x, sen1_desc95$all_mean+sen1_desc95$all_sd))[,1]))
  poly.y <- c(na.omit(cbind(sen1_desc95$day.x, sen1_desc95$all_mean+sen1_desc95$all_sd))[,2], rev(na.omit(cbind(sen1_desc95$day.x, sen1_desc95$all_mean-sen1_desc95$all_sd))[,2]))
  polygon(poly.x, poly.y, col=adjustcolor("orange", alpha.f = 0.2), border=NA)
  lines(na.omit(cbind(sen1_desc95$day.x, sen1_desc95$all_mean)), type="b", col="orange", lwd=3, pch=1, xaxt="n", xlab="")
  #axis(1, at=sen1_desc95$day.x[!is.na(sen1_desc95$all_mean)], sen1_desc95$label.x[!is.na(sen1_desc95$all_mean)], cex.axis=1.5, las=2)
  
  # 4) DESC 168
  sen1_des168 <- read.csv(sen1_data[grepl("DESC_168", sen1_data)], header=TRUE)
  
  par(new=T)
  plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5) 
  poly.x <- c(na.omit(cbind(sen1_des168$day.x, sen1_des168$all_mean+sen1_des168$all_sd))[,1], rev(na.omit(cbind(sen1_des168$day.x, sen1_des168$all_mean+sen1_des168$all_sd))[,1]))
  poly.y <- c(na.omit(cbind(sen1_des168$day.x, sen1_des168$all_mean+sen1_des168$all_sd))[,2], rev(na.omit(cbind(sen1_des168$day.x, sen1_des168$all_mean-sen1_des168$all_sd))[,2]))
  polygon(poly.x, poly.y, col=adjustcolor("green", alpha.f = 0.2), border=NA)
  lines(na.omit(cbind(sen1_des168$day.x, sen1_des168$all_mean)), type="b", col="green", lwd=3, pch=1, xaxt="n", xlab="")
  axis(1, at=c(1, 32, 62, 93, 123), c("01.03.", "01.04.", "01.05.", "01.06.", "01.07."), cex.axis=1.5)
  abline(v=c(1, 32, 62, 93, 123), col="gray")
  
  # Precipitation
  par(new=TRUE)
  plot(1, type="n", xlab="", ylab="", ylim=c(0,100), xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5)
  lines(sen1_asc146$day.x, prec_data_demmin$prec_daily, type="h", col="blue", lwd=3, xlab="", ylab="")
  lines(sen1_asc146$day.x, prec_data_bloen$prec_daily, type="h", col="deepskyblue", lwd=3, xlab="", ylab="")
  axis(side=4, cex.axis=1.2, line=1)
  mtext(side = 4, line=3, "Precipitation in mm / Soil Moisture in Vol%", cex=1.2)
  
  # Temperature
  par(new=TRUE)
  plot(1, type="n", xlab="", ylab="", ylim=c(0,2000), xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5)
  lines(sen1_asc146$day.x, temp_demmin$gdd_kumm, type="l", col="brown3", lwd=3, xlab="", ylab="")
  lines(sen1_asc146$day.x, temp_bloen$gdd_kumm, type="l", col="coral1", lwd=3, xlab="", ylab="")
  axis(side=4, cex.axis=1.2, line=4.5)
  mtext(side = 4, line=6.5, "Growing Degree Days", cex=1.2)
  
  # Soil Moisture
  par(new=TRUE)
  plot(1, type="n", xlab="", ylab="", ylim=c(0,100), xlim=c(0,153), axes=FALSE, xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5)
  lines(na.omit(cbind(sen1_asc146$day.x, soilmoisture$demmin)), type="l", col="seagreen", lwd=3, xlab="", ylab="")
  lines(na.omit(cbind(sen1_asc146$day.x, soilmoisture$bloen)), type="l", col="springgreen3", lwd=3, xlab="", ylab="")

  # Legend
  #output <- paste("Z:/Katharina/AgriFusion/Analyse/Sen1/TimeSeries_Sen1/Plots_Field_Mean_all/legend.png", sep="")
  #legend("topleft", labels, lty=1, lwd=3, col=c("blue", "red", "orange", "green"), bty='n', bg="white", cex=1.5, seg.len=1, x.intersp=0.5)
  #legend("topright", c("Precipitation Demmin", "Precipitation Blönsdorf"), lty=1, lwd=3, col=c("blue", "deepskyblue"), bty='n', bg="white", cex=1.5, seg.len=1, x.intersp=0.5)
  #legend("topleft", c("Growing Degree Days Demmin", "Growing Degree Days Blönsdorf"), lty=1, lwd=3, col=c("brown3", "coral1"), cex=1.5, seg.len=1, x.intersp=0.5)
  #legend("topleft", c("Soil Moisture Demmin", "Soil Moisture Blönsdorf"), lty=1, lwd=3, col=c("seagreen", "springgreen3"), cex=1.5, seg.len=1, x.intersp=0.5)
  dev.off()
  
}