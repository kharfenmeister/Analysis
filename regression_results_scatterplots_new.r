### Sentinel-1 vs- Fielddata Plot - ASC/DESC ###

croptypes <- c("wheat", "barley")

for (ct in 1:length(croptypes)){
  
  crop <- croptypes[ct]
  if(crop == "wheat"){
    fields <- c("BuBio", "A", "D", "S", "Hey", "E", "G")
  } else {
    fields <- c("KarBio", "B", "C", "F", "H", "Sie")
  }
  
  sen1_pass <- c("ASC_44", "ASC_146", "DESC_95", "DESC_168")
  
  sen1_params <- c("SigmaVH", "SigmaVV", "SigmaVHVV")
  veg_params <- c("biomass_yield", "dry_mass_yield", "vwc", "vwc_rel", "LAI", "plant_height")
  
  combis <- expand.grid(sen1_params, veg_params)
  names(combis) <- c("sen1", "veg_params")
  
  
  for (p in 1:length(sen1_pass)){
    pass <- strsplit(sen1_pass[p], "_")[[1]][1]
    orbit <- strsplit(sen1_pass[p], "_")[[1]][2] 
    
    output_path <- paste("Z:/Katharina/AgriFusion/Analyse/Sen1/Sen1_vs_Fielddata/plots/FINAL/linear/all/", gsub("_", "", sen1_pass[p]), "/all/", sep="")
    results_path <- "Z:/Katharina/AgriFusion/Analyse/Sen1/Sen1_vs_Fielddata/Results_FINAL/linear/"
    
    data_csv <- list.files("Z:/Katharina/AgriFusion/Analyse/Sen1/Sen1_vs_Fielddata/data/linear", pattern=paste("_", fields, "_", sep="", collapse="|"), full.names=TRUE)
    
    ## chose only images with fitting pass and orbit
    sen1_info <- list.files("Z:/Katharina/AgriFusion/Daten/Sen1/", pattern=".csv", full.names=TRUE)
    
    sen1_dates_x <- c()
    for (d in 1:length(sen1_info)){
      sen1_info_d <- read.csv(sen1_info[d], header=TRUE, sep=";")
      sen1_dates_fit <- sen1_info_d$Date[sen1_info_d$Path == pass & sen1_info_d$Orbit == orbit]
      sen1_dates_fit <- gsub("-", "", as.character(as.Date(sen1_dates_fit, format="%d.%m.%Y")))
      sen1_dates_x <- c(sen1_dates_x, sen1_dates_fit)
    }
    
    sen1_dates_x <- unique(sen1_dates_x)
    
    data_csv <- data_csv[grepl(paste(sen1_dates_x, collapse="|"), data_csv)]
    
    data_all <- read.csv(data_csv[1], header=TRUE)
    for (i in 2:length(data_csv)){
      data_i <- read.csv(data_csv[i], header=TRUE)
      data_all <- rbind(data_all, data_i)
    }
    
    
    data_all_2017 <- data_all[grepl("2017", data_all$sen1_date),]
    data_all_2018 <- data_all[grepl("2018", data_all$sen1_date),]
    
    outplot <- paste(results_path, "/plots/Scatterplots_", croptypes[ct], "_", gsub("_", "", sen1_pass[p]), "_all_new.png", sep="")
    
    plot.new()
    windowsFonts(Calibri=windowsFont("Calibri"))
    png(outplot, width=2500, height=4000, res=350)
    par(mfrow=c(6,6), oma = c(3, 4, 3, 0), mar = c(1, 1, 1, 1), yaxs="i", family="Calibri")
    
    for (c in 1:nrow(combis)){

      
      ### Splitted in Years
      
      sen1_param <- as.character(combis$sen1[c])
      veg_param <- as.character(combis$veg_params[c])
      
      xy_2017 <- cbind(data_all_2017[,sen1_param], data_all_2017[,veg_param])
      xy_2017 <- as.data.frame(na.omit(xy_2017))
      xy_2017 <- xy_2017[is.finite(rowSums(xy_2017)),]
      xy_2017 <- xy_2017[xy_2017[,2]!=0,]
      
      x_2017 <- xy_2017[,1]
      y_2017 <- xy_2017[,2]
      
      xy_2018 <- cbind(data_all_2018[,sen1_param], data_all_2018[,veg_param])
      xy_2018 <- as.data.frame(na.omit(xy_2018))
      xy_2018 <- xy_2018[is.finite(rowSums(xy_2018)),]
      xy_2018 <- xy_2018[xy_2018[,2]!=0,]
      
      x_2018 <- xy_2018[,1]
      y_2018 <- xy_2018[,2]
      
      if(nrow(xy) > 0){
        ## 1) Linear Model
        data_lm_2017 <- lm(y_2017~x_2017)
        data_lm_log_2017 <- lm(log(y_2017)~x_2017)
        #data_nls <- nls(y~a*exp(b*x)+d, start=c(a=3000, b=0.1, d=0.1))
        #data_glm <- glm(y~x, family=gaussian(link="identity"))
        
        # R²
        data_cor <- round(cor(x_2017, y_2017),2)
        rsquared_lin <- round(data_cor^2,3)
        data_cor_log <- round(cor(x_2017, log(y_2017)),2)
        rsquared_exp <- round(data_cor_log^2,3)
        
        
        # RMSE
        values_pred_2017 <- data_lm_2017$coefficients[1] + data_lm_2017$coefficients[2]*x_2017
        rmse_lin_2017 <- round(sqrt(mean((values_pred_2017 - y_2017)^2, na.rm = TRUE)),3)
        
        values_pred_log_2017 <- data_lm_log_2017$coefficients[1] + data_lm_log_2017$coefficients[2]*x_2017
        rmse_exp_2017 <- round(sqrt(mean((exp(values_pred_log_2017) - y_2017)^2, na.rm = TRUE)),3)
        
        data_lm_2018 <- lm(y_2018~x_2018)
        data_lm_log_2018 <- lm(log(y_2018)~x_2018)
        #data_nls <- nls(y~a*exp(b*x)+d, start=c(a=3000, b=0.1, d=0.1))
        #data_glm <- glm(y~x, family=gaussian(link="identity"))
        
        # R²
        data_cor <- round(cor(x_2018, y_2018),2)
        rsquared_lin <- round(data_cor^2,3)
        data_cor_log <- round(cor(x_2018, log(y_2018)),2)
        rsquared_exp <- round(data_cor_log^2,3)
        
        
        # RMSE
        values_pred_2018 <- data_lm_2018$coefficients[1] + data_lm_2018$coefficients[2]*x_2018
        rmse_lin_2018 <- round(sqrt(mean((values_pred_2018 - y_2018)^2, na.rm = TRUE)),3)
        
        values_pred_log_2018 <- data_lm_log_2018$coefficients[1] + data_lm_log_2018$coefficients[2]*x_2018
        rmse_exp_2018 <- round(sqrt(mean((exp(values_pred_log_2018) - y_2018)^2, na.rm = TRUE)),3)
        
        # Plot
        #output <- paste(output_path, veg_param, "/", crop, "_", sen1_param, "_", veg_param, "_", gsub("_", "", sen1_pass[p]), ".png", sep="")
        
        # Set ylab
        if (veg_param == "crop_coverage"){
          ylim <- c(0,100)
          ylab <- "Crop Coverage in %"
        } else if (veg_param == "LAI"){
          ylim <- c(0,9)
          ylab <- "LAI"
        } else if (veg_param == "dry_mass_yield"){
          ylim <- c(0,3000)
          ylab <- "Dry Biomass in g/m²"
        } else if (veg_param == "biomass_yield"){
          ylim <- c(0,10000)
          ylab <- "Biomass in g/m²"
        } else if (veg_param == "plant_height"){
          ylim <- c(0,100)
          ylab <- "Plant Height in cm"
        } else if (veg_param == "vwc"){
          ylim <- c(0,5000)
          ylab <- "VWC in g/m²"
        } else if (veg_param == "BBCH"){
          ylim <- c(0,100)
          ylab <- "BBCH"
        } else if (veg_param == "vwc_rel"){
          ylim <- c(0,100)
          ylab <- ylab <- "VWC in %"
        }
        
        # Set xlim and label
        if (sen1_param == "SigmaVHVV"){
          xlim <- c(-14, 1)
          sen1_param_label <- "VH/VV"
        } else if (sen1_param == "SigmaVH"){
          xlim <- c(-26,-10)
          sen1_param_label <- "VH"
        } else if (sen1_param == "SigmaVV"){
          xlim <- c(-23,-4)
          sen1_param_label <- "VV"
        } else if (sen1_param == "Entropy" | sen1_param == "Anisotropy" | sen1_param == "m1" | sen1_param == "m2" | sen1_param == "m3" | sen1_param == "m4"){
          xlim <- c(0,1)
        } else if (sen1_param == "Alpha"){
          xlim <- c(0,90)
        }
      
        
        print(paste("Min: ", veg_param, " ", sen1_param, ": ", round(min(x_2017),2), " ", round(min(x_2018),2), sep=""))
        print(paste("Max: ", veg_param, " ", sen1_param, ": ", round(max(x_2017),2), " ", round(max(x_2018),2), sep=""))
        
        if (veg_param == "biomass_yield" && sen1_param == "SigmaVH"){
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          mtext(bquote(bold(.(paste(sen1_param_label, "2017", sep=" ")))), side=3, line=1, font.lab=2 )
          mtext(bquote(bold(.(ylab))), side=2, line=3, font.lab=2 )
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          mtext(bquote(bold(.(paste(sen1_param_label, "2018", sep=" ")))), side=3, line=1, font.lab=2 )
          #points(x_2018, y_2018, pch=19, col="gray")
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
        } else if (veg_param == "plant_height" && sen1_param == "SigmaVH"){
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          mtext(bquote(bold(.(ylab))), side=2, line=3, font.lab=2 )
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2) 
        } else if (sen1_param == "SigmaVH"){
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          mtext(bquote(bold(.(ylab))), side=2, line=3)
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)   
        } else if (veg_param == "plant_height"){
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
        } else if (veg_param == "biomass_yield"){
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          mtext(bquote(bold(.(paste(sen1_param_label, "2017", sep=" ")))), side=3, line=1, font.lab=2 )
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          mtext(bquote(bold(.(paste(sen1_param_label, "2018", sep=" ")))), side=3, line=1, font.lab=2 )
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
        } else {
          par(mar=c(1,1,1,0))
          plot(x_2017, y_2017, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19)
          #abline(data_lm_2017, col="red", lwd=2)
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2017, list(x_2017=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
          #mtext(paste("R² = ", rsquared_lin, "\n", "RMSE = ", rmse_lin, sep=""), side=3, line=1, adj=0, cex=1.3, col="red")
          #mtext(paste("R² = ", rsquared_exp, "\n", "RMSE = ", rmse_exp, sep=""), side=3, line=1, adj=1, cex=1.3, col="green")
          par(mar=c(1,0,1,1))
          plot(x_2018, y_2018, xlab=sen1_param, ylab=ylab, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", cex.axis=1.3, cex.lab=1.3, pch=19, col="gray")
          #abline(data_lm_2018, col="red", lwd=2)  
          #lines(seq(xlim[1],xlim[2], 1), exp(predict(data_lm_log_2018, list(x_2018=seq(xlim[1],xlim[2], 1)))), col="green", lwd=2)
        }
        
        
      } else {
        print(paste(fields[f], " ", sen1_param, " ", veg_param, ": No data", sep=""))
      }
      
      
    } # combis
    
    dev.off()
  } # pass/orbit
  
} # croptypes
