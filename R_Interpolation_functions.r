## interpolation functions ##
#########################################################################################################################################
### list of functions in this file ####
# 	1.	getTimeSteps				get vector with time steps for loading data
#	2.	loadData.sums				load data from time intervals and calculate sums (e.g. for Precipitation)
#	3.	loadData.averages			load data from time intervals and calculate averages (e.g. for Barometric Pressure)
#	4.	standardElevation			standardize signal with elevation data to 18m
#	5.	kriging					execute kriging algorithmus for loaded data in specific time intervals
#	6.	correctElevation			correct interpolated signal for every pixel in kriging output with elevation data 
#	7.	layoutInterpolation			add layout with description and logos to kriging output 
#	8.	stackInterpolation			layerstack for all files from kriging output
#	9.	extractInterpolationValues		extract interpolated values from kriging output
#	10.	extractData.sums			extract data from array and calculate sums (e.g. for Precipitation)
#	11.	extractData.averages			extract data from array and calculate averages (e.g. for Barometric Pressure)
#	12.	pasteErrorRaster			paste error raster with values -999 if data is missing or errors occur
#	13.	createGrid				create grid for kriging function
#########################################################################################################################################
## constraints
# additional input files are required. Therefore relative paths are used. Data is stored in:
#	../../97_projectdata_tereno/
#	../../98_logos/

getTimeSteps <- function(start, end, interval){

# author: Katharina Heupel
# date: 2014-04-25
# update: 2014-05-14; Christian Hohmann
# optimierung
# input:
#	start		first date as string in format "yyyy-mm-dd hh:mm:ss"
#	end		last date as string in format "yyyy-mm-dd hh:mm:ss"
#	interval	time interval in seconds as numeric
# output:
#	zeiten		vector with all time steps from start to end 

start <- strptime(start, "%Y-%m-%d %H:%M:%S", tz="UTC")
end <- strptime(end, "%Y-%m-%d %H:%M:%S", tz="UTC")
s <- difftime(end, start, unit="secs")

zeiten <- as.POSIXlt(seq.int(start,end,interval))

datetime <- zeiten + ((interval/2)-1)
zeiten <- factor(zeiten,levels=zeiten)      

return(list(zeiten=zeiten,datetime=datetime))

}

loadData.sums <- function(coords, dir, sep, dec, err, skip, start, end, prefix, signal, gen_fun_file){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	zeiten		vector with time steps generated from function getTimeSteps
#	coords		file with coordinates of climate stations / BF-stations
#	dir		file directory
#	sep		column separator
#	dec		decimal separator
#	err		error code, to be replaced by NA's
#	skip		nb lines to skip
#	start		first date of data set
#	end		last date of data set
#	prefix		prefix in filename before station name
#	signal		signal to be loaded
#	gen_fun_file	path to generic functions file
# output:
#	data		matrix with coordinates and signals 



raw <- loadQualitycontrolledSignal(dir, sep, dec, err, skip, start, end, prefix, signal, gen_fun_file)
data <- raw[[1]]
error <- raw[[2]]
data[error<0] <- NA
sums <- colSums(data, na.rm=TRUE) + ifelse(colSums(is.na(data)) > nrow(data)-10, NA, 0)
sums <- as.matrix(sums)
colnames(sums)[1] <- signal

coords[,ncol(coords)+1] <- NA
k <- 1
while (k <= nrow(coords)){
  coords[k,ncol(coords)] <- sums[match(coords[k,1],rownames(sums)),1]  
  k <- k+1
  }
data <- coords[!is.na(coords[,ncol(coords)]),]

return(data)

}

loadData.averages <- function(zeiten, coords, dir, sep, dec, err, skip, start, end, prefix, signal, gen_fun_file){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	zeiten		vector with time steps generated from function getTimeSteps
#	coords		file with coordinates of climate stations / BF-stations
#	dir		file directory
#	sep		column separator
#	dec		decimal separator
#	err		error code, to be replaced by NA's
#	skip		nb lines to skip
#	start		first date of data set
#	end		last date of data set
#	prefix		prefix in filename before station name
#	signal		signal to be loaded
#	gen_fun_file	path to generic functions file
# output:
#	data		matrix with coordinates and signals 

start <- zeiten[i]
end <- zeiten[i+1]

raw <- loadQualitycontrolledSignal(dir, sep, dec, err, skip, start, end, prefix, signal, gen_fun_file)
data <- raw[[1]]
error <- raw[[2]]
data[error<0] <- NA
averages <- colMeans(data[,2:ncol(data)],na.rm=TRUE)
averages <- as.matrix(averages)
colnames(averages)[1] <- signal

coords[,ncol(coords)+1] <- NA
k <- 1
while (k <= nrow(coords)){
  coords[k,ncol(coords)] <- averages[match(coords[k,1],rownames(averages)),1]
  k <- k+1
  }
data <- coords[!is.na(coords[,ncol(coords)]),]

return(data)

}

standardElevation <- function(data, elev){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	data		output matrix from loadData
#	elev		vector with elevation of climate stations / BF-stations
# output:
#	data.sE		data + standardized signal in last column

bar <- data[,ncol(data)]

bar_norm <- ((elev-18)*0.125)+bar
data.sE <- cbind(data, bar_norm)

return(data.sE)

}

kriging <- function(data, sig, grid, output.kriging){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	data		output matrix from loadData
#	sig 		column of data with signal values (usually data[,5], after standardElevation data[,6])
#	grid		grid as base for kriging
#	output.kriging	path to output raster
# output:
#	krigout		kriging output as raster

for (j in 2:length(data)){
  data[,j] <- as.numeric(as.character(data[,j]))
  }

data <- na.omit(data)

coordinates(data)<-c("X","Y")
proj4string(data) <- CRS("+proj=utm +zone=33 ellps=WGS84")

vg <- variogram(sig~1,data,cutoff=28000,width=4000) 
vg.fit <- fit.variogram(vg, model=vgm(model="Sph", range=14000, nugget=TRUE))

okrige <- krige(sig~1, data, grid, model=vg.fit)

krigout_raster <- raster(okrige["var1.pred"])

krigout_values <- getValues(krigout_raster)
krigout_values[krigout_values<0] <- 0 

krigout <- setValues(krigout_raster, krigout_values)
krigout <- writeRaster(krigout, output.kriging, "GTiff", overwrite=TRUE)

return(krigout)

}

correctElevation <- function(krigout, dgm, output.elev){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	krigout		output from kriging function (raster)
#	dgm	 	Digital Elevation Model with same extent and cellsize as krigout
#	output.elev	path to output raster
# output:
#	krigout.elev	kriging output with consideration of elevation data as raster

dgm_values <- getValues(dgm)
krigout_values <- getValues(krigout)

bar_abs <- ((dgm_values-18)*(-0.125))+krigout_values

krigout.elev <- setValues(krigout, bar_abs)
proj4string(krigout.elev) <- CRS("+proj=utm +zone=33 ellps=WGS84")
krigout.elev <- writeRaster(krigout.elev, output.elev, "GTiff", overwrite=TRUE)

return(krigout.elev)

}

layoutInterpolation <- function(data, krigout, coords, coords.names, datetime, signal, output.layout, extent.layout){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	data		output matrix from loadData with column "Name" 
#	krigout 	output from kriging function (raster)
#	coords		file with coordinates of climate stations
#	coords.names	file with names(=labels) of climate stations and coordinates of point labels (columns: X, Y, Name)
#	datetime	date and time of kriging output for plot description text as string
#	signal		signal and unit for plot description as string
#	output.layout	path to output directory
#	extent.layout		matrix with border coordinates (c(x1, x2, y1, y2))
# output:
#	PNG-file in output directory

  coordinates(data)<-c("X","Y")
  proj4string(data) <- CRS("+proj=utm +zone=33 ellps=WGS84")

  coordinates(coords)<-c("X","Y")
  proj4string(coords) <- CRS("+proj=utm +zone=33 ellps=WGS84")

  coordinates(coords.names)<-c("X","Y")
  proj4string(coords.names) <- CRS("+proj=utm +zone=33 ellps=WGS84")

# Extent Raster
  b <- as(extent(extent.layout), 'SpatialPolygons')
  proj4string(b) <- CRS("+proj=utm +zone=33 ellps=WGS84")
  krigout_b <- crop(krigout, b)
  #Layout Settings
  sp.stationen <- list("sp.points", coords, pch=1, col="black")
  sp.stationen.use <- list("sp.points", data, pch=16, col="black")
  #sp.stationen.label <- list("sp.pointLabel", coords.names, label=coords.names$Name, cex=1)
  north_arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(extent.layout[2]-2000,extent.layout[3]+1000), scale = 3000)
  scale <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(extent.layout[1]+2500, extent.layout[3]+1000), scale=5000, fill=c("transparent","black"))
  scale.text2 <- list("sp.text", c(extent.layout[1]+4900,extent.layout[3]+2000), "5 km")
  layout <- list(sp.stationen, sp.stationen.use, north_arrow, scale, scale.text2)
  for (i in 1:nrow(coords.names)){
    sp.stationen.label <- list("sp.text", c(coords.names$X[i], coords.names$Y[i]), coords.names$Name[i])
    layout [[length(layout)+1]]<- sp.stationen.label
  }

  # Set color depending on signal
  if (signal == "Precipitation"){
    col <- colorRampPalette(brewer.pal(9,"Blues"))(256)
  } else if (signal == "Temperature" | signal == "CMP3oben"){
    col <- rev(heat.colors(256))
  } else {
    col <- rev(rainbow(256))
  }

  # Set scale depending on signal
  if (signal == "Precipitation"){
    scale.start <- 0
    scale.end <- 50
    scale.seq <- (log(scale.end)-log(scale.start+0.1))
    at <- exp(seq(log(scale.start+0.1), log(scale.end), scale.seq/255))
  } else if (signal == "BarometricPressure"){
    scale.start <- 965
    scale.end <- 1025
    scale.seq <- scale.end - scale.start
    at <- seq(scale.start, scale.end, scale.seq/255)
  } else if (signal == "CMP3oben") {
    scale.start <- 0
    scale.end <- 1000
    scale.seq <- scale.end - scale.start
    at <- seq(scale.start, scale.end, scale.seq/255)
  }

  # spplots von Kriging-Output
  #scale.seq <- scale.end - scale.start
  #at <- seq(scale.start, scale.end, scale.seq/256)
  spplot.krigout <- spplot(krigout_b, at=at, scales=list(draw=T), col.regions=col, cuts=255, sp.layout=layout)

  # Logos laden
  tereno <- readJPEG("../../98_logos/TERENO_logo.jpg")
  gfz <- readPNG("../../98_logos/GFZ_Logo_png.png")
  dlr <- readJPEG("../../98_logos/DLR_logo.jpg")
  dtl <- readPNG("../../97_projectdata_tereno/Projektgebiet_Dtl_lin.png")

  # Beschreibung abhängig von Signal (Input Parameters!)
  if (signal == "BarometricPressure"){
    plot_text <- paste("Integrated Stations:", nrow(data), "\nInput Parameters: coordinates, measurements, \n			      DEM100 \nMap Projection: UTM, WGS84 \nSoftware: R \nInterpolation Algorithm: Kriging \n \n     data available and integrated \n     no data available ")
  } else {
    plot_text <- paste("Integrated Stations:", nrow(data), "\nInput Parameters: coordinates, measurements \n \nMap Projection: UTM, WGS84 \nSoftware: R \nInterpolation Algorithm: Kriging \n \n      data available and integrated \n      no data available ")
  }
  # Rand + Text hinzufügen (width/height = Umrechnung der Weite/Länge des Plots in km und Pixel (1km = 22 Pixel) und Hinzufügen eines festen Randes (350 bzw. 450 px) für Plot-Umgebung)
  png(file=output.layout, width=(((extent.layout[2]-extent.layout[1])/1000)*22)+450, height=(((extent.layout[4]-extent.layout[3])/1000)*22)+350)
  plot.new()
  par(xpd=NA, lheight=1.2, usr=c(0,1,0,1))
  print(spplot.krigout, position=c(0,0,0.7,1))
  mtext("Published by: GeoForschungsZentrum Potsdam \nAuthors: Katharina Heupel, Christian Hohmann \nDevelopment Date: Jan-Jun 2014", side=1, outer=F, line=0.5, at=0.15, adj=0, font=11, cex=1.2)
  text(0.85,0.65, paste(signal, "\n", datetime_i, sep=""), font=11, cex=1.8)
  text(0.7,0.14, plot_text, font=11, pos=4, cex=1.2)
  points(0.72, 0.098, pch=16, col="black", cex=1.3)
  points(0.72, 0.073, pch=1, col="black", cex=1)
  rasterImage(tereno, 0.75, 0.75, 0.95, 0.825)
  rasterImage(gfz, 0, -0.075, 0.15, 0.075)
  rasterImage(dlr, 0.5, -0.04, 0.6, 0.06)
  rasterImage(dtl, 0.75, 0.3, 0.95, 0.6)
  segments(0.7, 0.6, 1, 0.6)
  segments(0.7, 0.7, 1, 0.7)
  segments(0.7, 0.3, 1, 0.3)
  rect(-0.05,-0.07,1.02,0.93)
  dev.off()
}

stackInterpolation <- function(krigout.dir, layer.names, output.stack){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	krigout.dir	directory to kriging output raster
#	layer.names	vector with names for layer
#	output.stack		output directory
# output:
#	layerstack in output directory

stack <- stack(krigout.dir)
names(stack) <- layer.names
writeRaster(stack, output.stack, "GTiff", overwrite=TRUE)

}

extractInterpolationValuesStack <- function(coords, input, start, end, output.extract){

# author: Katharina Heupel
# date: 2014-04-25
# input:
#	coords		SpatialPointsDataFrame with Extraction-Points and columns "Name", "X", and "Y"
#	input		raster or layer stack to extract values
#	start		first date of extraction period as string in form "yyyy-mm-dd hh:mm:ss"
#	end		last date of extraction period as string in form "yyyy-mm-dd hh:mm:ss"
#	output.extract		output filename	
# output:
#	CSV-file in output directory

coordinates(coords)<-c("X","Y")
proj4string(coords) <- CRS("+proj=utm +zone=33 ellps=WGS84")

start <- strptime(start, "%Y-%m-%d %H:%M:%S", tz="UTC")
end <- strptime(end, "%Y-%m-%d %H:%M:%S", tz="UTC")
doy_start <- as.numeric(strftime(start, format = "%j"))
doy_end <- as.numeric(strftime(end, format = "%j"))
year <- substr(as.character(start),1,4)

data <- extract(input, coords, df=TRUE)
rownames(data) <- coords$Name
data_t <- t(data)

data_t <- data_t[-1,]
data_t_sel <- data_t[doy_start:doy_end,]

write.csv(data_t_sel, output.extract)

}

extractInterpolationValuesList <- function(extract.coords, input, output.extract){

# author: Katharina Heupel
# date: 2014-05-08
# input:
#	extract.coords	Table with Extraction-Points and columns "Name", "X", and "Y"
#	input		list with raster files
#	output.extract	output filename	
# output:
#	CSV-file in output directory

coordinates(extract.coords) <- c("X","Y")
proj4string(extract.coords) <- CRS("+proj=utm +zone=33 ellps=WGS84")

extraction <- as.data.frame(extract.coords)

for (i in 1:length(input)){
  a <- extract(raster(input[i]), extract.coords, df=TRUE)
  extraction <- cbind(extraction, a[2])
}

extraction <- t(extraction)

write.csv(extraction, output.extract)

}

extractData.sums <- function(coords, data, start, end){

# author: Christian Hohmann
# date: 2014-05-14
# input:
#	coords		file with coordinates of climate stations / BF-stations
#	start		first date of data set
#	end		last date of data set
#	data		array with numeric data, first colum is timestamp
# output:
#	data		matrix with coordinates and signals 

data<-data[data[,1]>=as.numeric(as.POSIXct(start)),]
data<-data[data[,1]<as.numeric(as.POSIXct(end)),]
sums <- colSums(data, na.rm=TRUE) + ifelse(colSums(is.na(data)) > nrow(data)-10, NA, 0)
sums <- as.matrix(sums)
colnames(sums)[1] <- signal

coords[,ncol(coords)+1] <- NA
k <- 1
while (k <= nrow(coords)){
  coords[k,ncol(coords)] <- sums[match(coords[k,1],rownames(sums)),1]  
  k <- k+1
  }
data <- coords[!is.na(coords[,ncol(coords)]),]

return(data)

}

extractData.averages <- function(coords, data, start, end){

# author: Katharina Heupel
# date: 2014-05-16
# input:
#	coords		file with coordinates of climate stations / BF-stations
#	start		first date of data set
#	end		last date of data set
#	data		array with numeric data, first colum is timestamp
# output:
#	data		matrix with coordinates and signals 

data <- data[data[,1] >= as.numeric(as.POSIXct(start)),]
data <- data[data[,1] < as.numeric(as.POSIXct(end)),]
averages <- colMeans(data[,2:ncol(data)],na.rm=TRUE)
averages <- as.matrix(averages)
colnames(averages)[1] <- signal

coords[,ncol(coords)+1] <- NA
k <- 1
while (k <= nrow(coords)){
  coords[k,ncol(coords)] <- averages[match(coords[k,1],rownames(averages)),1]
  k <- k+1
  }
data <- coords[!is.na(coords[,ncol(coords)]),]

return(data)

}

pasteErrorRaster <- function(extent, cellsize, output.error){

# author: Katharina Heupel
# date: 2014-06-27
# input:
#	extent		matrix with border coordinates (c(x1, x2, y1, y2))
#	cellsize	size of gridcells
#	output.error	output path and filename
# output:
#	raster in output directory with error-value "-999" for each pixel

offset <- cellsize/2

error.raster <- raster(nrow=(extent[4]-extent[3])/cellsize, ncol=(extent[2]-extent[1])/cellsize, xmn=extent[1], xmx=extent[2], ymn=extent[3], ymx=extent[4])
error.raster <- setValues(error.raster, -999)
proj4string(error.raster) <- CRS("+proj=utm +zone=33 ellps=WGS84")

writeRaster(error.raster, output.error, "GTiff", overwrite=TRUE)

}

createGrid <- function(extent, cellsize){

# author: Katharina Heupel
# date: 2014-07-03
# input:
#	extent		matrix with border coordinates (c(x1, x2, y1, y2))
#	cellsize	size of gridcells
# output:
#	spatial grid 

offset <- cellsize/2

x <- seq(extent[1]+offset, extent[2]-offset, cellsize)
y <- seq(extent[3]+offset, extent[4]-offset, cellsize)

x <- rep(x, (extent[4]-extent[3])/cellsize)
y <- rep(y, each=(extent[2]-extent[1])/cellsize)

xy <- cbind(x, y)
xy.sp <- SpatialPoints(xy)
proj4string(xy.sp) <- CRS("+proj=utm +zone=33 ellps=WGS84")
gridded(xy.sp)=TRUE 

return(xy.sp)

}


