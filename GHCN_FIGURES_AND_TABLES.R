## The weather station locations to be analyzed
GHCN.stations <- c("USC00055531","USC00051886","USC00295084","USC00293031")
GHCN.stations.abbrev <- c("MVNP","CRTZ","LANL","ESPN")
GHCN.stations.long <- c("Mesa Verde National Park, Colorado","Cortez, Colorado","Los Alamos National Laboratory, New Mexico","Espanola, New Mexico")

## Get a SPDF of weather station metadata (for extracting locations from PRISM data)
GHCN.meta.sp <- getGHCNStationMetadataSPDF(GHCN.stations, data.dir=paste(MASTER.DATA,"GHCN/",sep=''))

## Local regression and output of only the GHCN locations ##
VEPIIN.GHCN.cellNums <- raster::extract(VEPIIN.PRISM.h20_year.prcp,GHCN.meta.sp, cellnumbers=T)[,"cells"]
names(VEPIIN.GHCN.cellNums) <- GHCN.meta.sp$ID
VEPIIS.GHCN.cellNums <- raster::extract(VEPIIS.PRISM.h20_year.prcp,GHCN.meta.sp, cellnumbers=T)[,"cells"]
names(VEPIIS.GHCN.cellNums) <- GHCN.meta.sp$ID
VEPIIN.GHCN.cellNums <- VEPIIN.GHCN.cellNums[!is.na(VEPIIN.GHCN.cellNums)]
VEPIIS.GHCN.cellNums <- VEPIIS.GHCN.cellNums[!is.na(VEPIIS.GHCN.cellNums)]

# Calculate the SLM for including all predictors over the calibration period
GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIN.PRISM.h20_year.prcp, predlist=predlist, name="GHCN_VEPIIN_PRISM_h20_year_prcp", cellNums=VEPIIN.GHCN.cellNums)
GHCN.VEPIIN.PRISM.h20_year.prcp.recons <- recon(name="GHCN_VEPIIN_PRISM_h20_year_prcp", coefficients=GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIN.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)

GHCN.VEPIIN.PRISM.grow.GDD.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIN.PRISM.grow.GDD, predlist=predlist, name="GHCN_VEPIIN_PRISM_grow_GDD", cellNums=VEPIIN.GHCN.cellNums)
GHCN.VEPIIN.PRISM.grow.GDD.recons <- recon(name="GHCN_VEPIIN_PRISM_grow_GDD", coefficients=GHCN.VEPIIN.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIN.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)

GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIS.PRISM.h20_year.prcp, predlist=predlist, name="GHCN_VEPIIS_PRISM_h20_year_prcp", cellNums=VEPIIS.GHCN.cellNums)
GHCN.VEPIIS.PRISM.h20_year.prcp.recons <- recon(name="GHCN_VEPIIS_PRISM_h20_year_prcp", coefficients=GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIS.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)

GHCN.VEPIIS.PRISM.grow.GDD.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIS.PRISM.grow.GDD, predlist=predlist, name="GHCN_VEPIIS_PRISM_grow_GDD", cellNums=VEPIIS.GHCN.cellNums)
GHCN.VEPIIS.PRISM.grow.GDD.recons <- recon(name="GHCN_VEPIIS_PRISM_grow_GDD", coefficients=GHCN.VEPIIS.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIS.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)


### TABLES SHOWING CAR REGRESSION RESULTS ###
dendros.study <- read.csv(paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_METADATA.csv", sep=''))

# A method to make a table of the results of CAR regression variable selection
table.segments.slm.model <- function(model, station, name){
  station <- which(GHCN.stations==station)
  site.abbrev <- GHCN.stations.abbrev[station] 
  site <- GHCN.stations.long[station]
  
  cars <- model$cars
  cars <- cars[[tail(which(as.numeric(names(cars)) < 1983), 1)]]
  cars <- cars^2
  cars <- cars[order(cars, decreasing=T)]
  cars <- cars[cars>0]
  
  betas <- model$betas
  betas <- betas[tail(which(as.numeric(rownames(betas)) < 1983), 1),]
  betas <- betas[names(cars)]
  
  if(grepl("prcp",name)){
    caption <- paste("\\textbf{CAR regression results for predicting net water-year precipitation at the ",site," weather station over the calibration period.}",
                     " $\\mathsf{n=",length(names(cars)),"}$ chronologies were selected by minimizing AIC.",
                                 " Summed CAR$\\mathsf{^2=",round(sum(cars),digits=3),"}$.",sep='')
  }else if(grepl("GDD",name)){
    caption <- paste("\\textbf{CAR regression results for predicting growing-season growing degree days (GDD) at the ",site," weather station over the calibration period.}",
                     " $\\mathsf{n=",length(names(cars)),"}$ chronologies were selected by minimizing AIC.",
                                 " Summed CAR$\\mathsf{^2=",round(sum(cars),digits=3),"}$.",sep='')
  }
  
  series.selected <- dendros.study[dendros.study$SERIES %in% names(cars),c("SERIES","NAME","SPECIES","ELEVATION")]
  
  regression.results <- data.frame(SERIES=names(betas),Car=cars,Beta=betas)
  final.table <- merge(series.selected,regression.results,by="SERIES", all=T)
  final.table <- final.table[order(-final.table$Car),]
  names(final.table) <- c("Series","Name","Species","Elevation (m)","CAR$^2$","\\multicolumn{1}{c}{$\\beta$}")
  ## Write mean absolute differences table to a text file
  fileConn <- file(paste("../TABLES/",name,"_",site.abbrev,".tex", sep=''))
  
  writeLines(print(
    xtable(
      final.table,
      caption = caption,
      label = paste("tab:",name,"_",site.abbrev, sep=''),
      digits=3,
      align="llllccr"),
    booktabs=TRUE,
    floating=TRUE,
    table.placement = "!h",
    caption.placement = "top",
    include.rownames=F,
    tabular.environment='tabular',
    sanitize.text.function=function(x){gsub("_","\\_",x,fixed=TRUE)},
    math.style.negative=T,
    hline.after = NULL,
    add.to.row = list(pos = list(-1,0,nrow(final.table)),
                      command = c("\\toprule \n",
                                  "\\midrule \n",
                                  "\\bottomrule \n"))),
    fileConn)
  close(fileConn)
}

tables <- lapply(1:length(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models),FUN=function(x){table.segments.slm.model(model=GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_h20_year_prcp")})
tables <- lapply(1:length(GHCN.VEPIIN.PRISM.grow.GDD.lm.models),FUN=function(x){table.segments.slm.model(model=GHCN.VEPIIN.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_grow_GDD")})
tables <- lapply(1:length(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models),FUN=function(x){table.segments.slm.model(model=GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_h20_year_prcp")})
tables <- lapply(1:length(GHCN.VEPIIS.PRISM.grow.GDD.lm.models),FUN=function(x){table.segments.slm.model(model=GHCN.VEPIIS.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_grow_GDD")})


### MAPS SHOWING CHRONLOGY SELECTION ###
WDB2 <- paste(MASTER.DATA,"WDB2",sep='')
NATIONAL_ATLAS <- paste(MASTER.DATA,"NATIONAL_ATLAS",sep='')

# A method to create a spatial polygon from the extent of a spatial* or raster* object
polygonFromExtent <- function(sp.object, buffer=0){
  x <- extent(sp.object)
  proj4string <- projection(sp.object)
  extent.matrix <- rbind(c(x@xmin-buffer,x@ymin-buffer), c(x@xmin-buffer,x@ymax+buffer), c(x@xmax+buffer,x@ymax+buffer), c(x@xmax+buffer,x@ymin-buffer), c(x@xmin-buffer,x@ymin-buffer) ) # clockwise, 5 points to close it
  extent.SP <- SpatialPolygons( list(Polygons(list(Polygon(extent.matrix)),"extent")), proj4string=CRS(proj4string) )
  return(extent.SP)
}

polygonUTM_NAD83 <- function(UTM_East,UTM_West,UTM_North,UTM_South,UTM_Zone){
  # Create matrix of coordinates
  datainUTM<-matrix(c(UTM_East, UTM_West, UTM_West, UTM_East,UTM_East, UTM_North,UTM_North,UTM_South,UTM_South,UTM_North),nrow=5)
  
  # Set universal projection
  master.proj <- CRS(paste("+proj=utm +datum=NAD83 +zone=",UTM_Zone,sep=''))
  
  # Create SpatialPolygon of simulation area
  sim.poly <- Polygons(list(Polygon(datainUTM, hole=FALSE)),ID='A')
  sim.poly <- SpatialPolygons(list(sim.poly), proj4string=master.proj)
  IDs <- sapply(slot(sim.poly, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  sim.poly <- SpatialPolygonsDataFrame(sim.poly,df)
  return(sim.poly)
}

# Lambert's conic conformal projection centered on the four corners
master.proj <- CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225")

states <- readOGR(paste(NATIONAL_ATLAS,"/statep010",sep=''), layer='statep010')
rivers <- readOGR(WDB2, layer="ne_10m_rivers_lake_centerlines")
lakes <- readOGR(WDB2, layer="ne_10m_lakes")

dendros.study <- read.csv(paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_METADATA.csv", sep=''))

dendros.study <- dendros.study[!is.na(dendros.study$LON),]
dendros.study <- SpatialPointsDataFrame(as.matrix(cbind(dendros.study$LON,dendros.study$LAT)),dendros.study, proj4string=CRS(projection(states)))

GHCN.stations.sp <- readOGR(paste(MASTER.DATA,"GHCN/ghcnd-stations.kml", sep=''),"ghcnd_stations")

states <- spTransform(states,master.proj)
rivers <- spTransform(rivers,master.proj)
lakes <- spTransform(lakes,master.proj)
dendros.study <- spTransform(dendros.study,master.proj)
GHCN.stations.sp <- spTransform(GHCN.stations.sp,master.proj)

states <- states[states$STATE %in% c("Arizona","Colorado","Utah","New Mexico"),]

sim.poly <- polygonFromExtent(states, buffer=20000)

rivers <- crop(rivers,sim.poly)
lakes <- crop(lakes,sim.poly)

VEPIIN <- spTransform(polygonUTM_NAD83(UTM_North = 4170000, UTM_South = 4102000, UTM_East = 740000, UTM_West = 672800, UTM_Zone = 12),master.proj)
VEPIIS <- spTransform(polygonUTM_NAD83(UTM_North = 4030400, UTM_South = 3939600, UTM_East = 435800, UTM_West = 359200, UTM_Zone = 13),master.proj)

# A method to plot the results of CAR regression variable selection
plot.segments.slm.model <- function(model, station, name){
  station.location <- GHCN.stations.sp[GHCN.stations.sp$ID==station,]
  station <- which(GHCN.stations==station)
  site.abbrev <- GHCN.stations.abbrev[station] 
  site <- GHCN.stations.long[station]
  
  cars <- model$cars
  cars <- cars[[tail(which(as.numeric(names(cars)) < 1983), 1)]]
  cars <- cars^2
  cars <- cars[order(cars, decreasing=T)]
  cars <- cars[cars>0]
  
  variable.locations <- dendros.study[dendros.study$SERIES %in% names(cars),]
  
  car <- as.vector(cars)[unlist(lapply(variable.locations$SERIES,FUN=function(x,...){which(names(cars) == x)}))]
  
  segments(x1=coordinates(station.location)[1,1],y1=coordinates(station.location)[1,2],x0=coordinates(variable.locations)[,1],y0=coordinates(variable.locations)[,2],lend=0, lwd=10*(car^1.5/sum(car^1.5)), col="gray50")
}



fig.ratio <- (ymax(extent(sim.poly))-ymin(extent(sim.poly)))/(xmax(extent(sim.poly))-xmin(extent(sim.poly)))

quartz(file="../FIGURES/CAR_SELECTION_PRCP.pdf", width=6, height = 6*fig.ratio, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=4, dpi=300)

par(mai=c(0.1,0.1,0.1,0.1), oma=c(0,0,0,0), lend=2, ljoin=1, xpd=F)

plot(1, type='n', xlab="", ylab="", xlim=c(xmin(extent(sim.poly)),xmax(extent(sim.poly))),ylim=c(ymin(extent(sim.poly)),ymax(extent(sim.poly))), xaxs="i", yaxs="i", axes=FALSE, main='')

par(ljoin=0, lend=2)

# plot(VEPIIN, add=T, col="gray80")
# plot(VEPIIS, add=T, col="gray80")

plot(rivers, add=T, col="black", lwd=0.75)
plot(lakes, add=T, col="white", lwd=0.75)
plot(states, add=T, border="black", lwd=1.5)

# lapply(1:length(prcp.slms),FUN=function(x){plot.segments.slm.model(model=prcp.slms[[x]],station=x)})

lapply(1:length(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_h20_year_prcp")})
lapply(1:length(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_h20_year_prcp")})

plot(dendros.study[dendros.study$SPECIES %in% c("PIPO","PIAR","PILO","PIFL","PICO", "PCEN","ABLO","ABCO"),], add=T, pch=19,cex=1.1, col="#e41a1c")
plot(dendros.study[dendros.study$SPECIES %in% c("PSME"),], add=T, pch=19,cex=1.1, col="#0052CC")
plot(dendros.study[dendros.study$SPECIES %in% c("PIED","PICM","JUSP"),], add=T, pch=19,cex=1.1, col="#4daf4a")

plot(GHCN.stations.sp, col="black", pch=17, cex=2, add=T)
text(GHCN.stations.sp[1,],"CRTZ", adj=c(1,-0.95), cex=1.2)
text(GHCN.stations.sp[2,],"MVNP", adj=c(1.18,1.5), cex=1.2)
text(GHCN.stations.sp[3,],"ESPN", adj=c(1,-1.5), cex=1.2)
text(GHCN.stations.sp[4,],"LANL", adj=c(0.55,-1.1), cex=1.2)
# plot(pollen, add=T, pch=1, cex=1.2, col="black")

# text(-91977.01,4386391,"VEPIIN", adj=c(-0.1,-0.5), cex=1.2)
# text(128269.4,4247709,"VEPIIS", adj=c(0,-0.5), cex=1.2)

label.coords <- coordinates(states)
label.coords[2,1] <- label.coords[2,1] + 175000
label.coords[3,1] <- label.coords[3,1] + 50000
text(label.coords, labels=states$STATE, cex=2.1)

legend("bottomleft", inset=-0.01,xpd=T, pch=c(19,19,19,17),col=c("#0052CC","#4daf4a","#e41a1c","black"), legend=c("Douglas fir","Pinyon and juniper","Spruce, pine, and true fir","GHCN stations"), bty = "n", cex=1.7)

dev.off()


quartz(file="../FIGURES/CAR_SELECTION_GDD.pdf", width=6, height = 6*fig.ratio, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=4, dpi=300)

par(mai=c(0.1,0.1,0.1,0.1), oma=c(0,0,0,0), lend=2, ljoin=1, xpd=F)

plot(1, type='n', xlab="", ylab="", xlim=c(xmin(extent(sim.poly)),xmax(extent(sim.poly))),ylim=c(ymin(extent(sim.poly)),ymax(extent(sim.poly))), xaxs="i", yaxs="i", axes=FALSE, main='')

par(ljoin=0, lend=2)

# plot(VEPIIN, add=T, col="gray80")
# plot(VEPIIS, add=T, col="gray80")

plot(rivers, add=T, col="black", lwd=0.75)
plot(lakes, add=T, col="white", lwd=0.75)
plot(states, add=T, border="black", lwd=1.5)

# lapply(1:length(prcp.slms),FUN=function(x){plot.segments.slm.model(model=prcp.slms[[x]],station=x)})

# lapply(1:length(GHCN.VEPIIN.PRISM.h20_year.prcp.slm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIN.PRISM.h20_year.prcp.slm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.h20_year.prcp.slm.models)[x]))), name="GHCN_VEPIIN_PRISM_h20_year_prcp")})
lapply(1:length(GHCN.VEPIIN.PRISM.grow.GDD.lm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIN.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_grow_GDD")})
# lapply(1:length(GHCN.VEPIIS.PRISM.h20_year.prcp.slm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIS.PRISM.h20_year.prcp.slm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.h20_year.prcp.slm.models)[x]))), name="GHCN_VEPIIS_PRISM_h20_year_prcp")})
lapply(1:length(GHCN.VEPIIS.PRISM.grow.GDD.lm.models),FUN=function(x){plot.segments.slm.model(model=GHCN.VEPIIS.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_grow_GDD")})


plot(dendros.study[dendros.study$SPECIES %in% c("PIPO","PIAR","PILO","PIFL","PICO", "PCEN","ABLO","ABCO"),], add=T, pch=19,cex=1.1, col="#e41a1c")
plot(dendros.study[dendros.study$SPECIES %in% c("PSME"),], add=T, pch=19,cex=1.1, col="#0052CC")
plot(dendros.study[dendros.study$SPECIES %in% c("PIED","PICM","JUSP"),], add=T, pch=19,cex=1.1, col="#4daf4a")

plot(GHCN.stations.sp, col="black", pch=17, cex=2, add=T)
text(GHCN.stations.sp[1,],"CRTZ", adj=c(1,-0.95), cex=1.2)
text(GHCN.stations.sp[2,],"MVNP", adj=c(-0.2,1.5), cex=1.2)
text(GHCN.stations.sp[3,],"ESPN", adj=c(-0.2,1.3), cex=1.2)
text(GHCN.stations.sp[4,],"LANL", adj=c(1.1,-0.7), cex=1.2)
# plot(pollen, add=T, pch=1, cex=1.2, col="black")

# text(-91977.01,4386391,"VEPIIN", adj=c(-0.1,-0.5), cex=1.2)
# text(128269.4,4247709,"VEPIIS", adj=c(0,-0.5), cex=1.2)

label.coords <- coordinates(states)
label.coords[2,1] <- label.coords[2,1] + 175000
label.coords[3,1] <- label.coords[3,1] + 50000
text(label.coords, labels=states$STATE, cex=2.1)

legend("bottomleft", inset=-0.01,xpd=T, pch=c(19,19,19,17),col=c("#0052CC","#4daf4a","#e41a1c","black"), legend=c("Douglas fir","Pinyon and juniper","Spruce, pine, and true fir","GHCN stations"), bty = "n", cex=1.7)

dev.off()


### TABLES SHOWING MODEL FIT STATISTICS ###
# This uses a 3-fold consecutive cross-validation over the calibration years (1924--1983).
# The same test is used for all other compared reconstruction methods
table.fits <- function(model, station, name, split=T){
  station <- which(GHCN.stations==station)
  site.abbrev <- GHCN.stations.abbrev[station] 
  site <- GHCN.stations.long[station]
  
  fits <- model$errors
  fits.mean <- do.call(rbind,lapply(fits,colMeans))
  fits.mean <- as.data.frame(fits.mean[,c("R2c","R2v","N_RMSE","RE","CE")])
  fits.mean$Year <- as.character(rownames(fits.mean))
  fits.mean <- fits.mean[,c("Year","R2c","R2v","N_RMSE","RE","CE")]
  
  
  if(grepl("prcp",name)){
    caption <- paste("\\textbf{Model fits for predicting net water-year precipitation at the ",site," weather station.}",
                     " Fits were generates using 3-fold consecutive cross-validation over the calibration period, 1924--1983.",sep='')
  }else if(grepl("GDD",name)){
    caption <- paste("\\textbf{Model fits for predicting growing-season growing degree days (GDD) at the ",site," weather station.}",
                     " Fits were generates using 3-fold consecutive cross-validation over the calibration period, 1924--1983.",sep='')
  }
  
  names(fits.mean) <- c("Start year (AD)","\\multicolumn{1}{c}{$R_c^2$}","\\multicolumn{1}{c}{$R_v^2$}","\\multicolumn{1}{c}{$RMSE_n$}","\\multicolumn{1}{c}{$RE$}","\\multicolumn{1}{c}{$CE$}")
  
  table <- print(
    xtable(
      fits.mean,
      caption = caption,
      label = paste("tab:MODEL_PERFORMANCE_",name,"_",site.abbrev, sep=''),
      digits=3,
      align="llrrrrr"),
#     size="\\footnotesize",
    print.results=F,
    booktabs=TRUE,
    floating=TRUE,
    table.placement = "!h",
    caption.placement = "top",
    include.rownames=F,
    tabular.environment='tabular',
    sanitize.text.function=function(x){gsub("_","_",x,fixed=TRUE)},
    math.style.negative=T,
    hline.after = NULL,
    add.to.row = list(pos = list(-1,0,nrow(fits.mean)),
                      command = c("\\toprule \n",
                                  "\\midrule \n",
                                  "\\bottomrule \n")))
  
#   split.year <- fits.mean[,"Start year"][floor(nrow(fits.mean)/2)+1]
#   table <- gsub("tabular}\\{","tabular}[t]\\{",table)
#   table <- gsub(paste("\n  ",split.year,sep=''),paste("\\\n \\\\bottomrule \\\n \\\\end{tabular} \\\n \\\\begin{tabular}[t]{lrrrrr} \\\n  \\\\toprule \\\n Start year & $R_c^2$ & $R_v^2$ & $NRMSE$ & $RE$ & $CE$ \\\\\\\\ \\\n  \\\\midrule \\\n  ",split.year, sep=''),table)
#   
  
  ## Write mean absolute differences table to a text file
  fileConn <- file(paste("../TABLES/MODEL_PERFORMANCE_",name,"_",site.abbrev,".tex", sep=''))
  
  writeLines(table,fileConn)
  close(fileConn)
  
  
}

tables <- lapply(1:length(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models),FUN=function(x){table.fits(model=GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_h20_year_prcp")})
tables <- lapply(1:length(GHCN.VEPIIN.PRISM.grow.GDD.lm.models),FUN=function(x){table.fits(model=GHCN.VEPIIN.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIN.GHCN.cellNums == as.numeric(names(GHCN.VEPIIN.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIN_PRISM_grow_GDD")})
tables <- lapply(1:length(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models),FUN=function(x){table.fits(model=GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.h20_year.prcp.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_h20_year_prcp")})
tables <- lapply(1:length(GHCN.VEPIIS.PRISM.grow.GDD.lm.models),FUN=function(x){table.fits(model=GHCN.VEPIIS.PRISM.grow.GDD.lm.models[[x]],station=names(which(VEPIIS.GHCN.cellNums == as.numeric(names(GHCN.VEPIIS.PRISM.grow.GDD.lm.models)[x]))), name="GHCN_VEPIIS_PRISM_grow_GDD")})


### COMPARE TO OTHER MODELS
compare.prcp.mvnp <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIN.PRISM.h20_year.prcp[5239]), regularize=regularize)
compare.gdd.mvnp <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIN.PRISM.grow.GDD[5239]), regularize=regularize)

compare.prcp.crtz <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIN.PRISM.h20_year.prcp[3629]), regularize=regularize)
compare.gdd.crtz <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIN.PRISM.grow.GDD[3629]), regularize=regularize)

compare.prcp.lanl <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIS.PRISM.h20_year.prcp[6894]), regularize=regularize)
compare.gdd.lanl <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIS.PRISM.grow.GDD[6894]), regularize=regularize)

compare.prcp.espn <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIS.PRISM.h20_year.prcp[5262]), regularize=regularize)
compare.gdd.espn <- cv.compare.methods(X=training.series, Y=as.numeric(VEPIIS.PRISM.grow.GDD[5262]), regularize=regularize)

prcp.compare <- list(compare.prcp.mvnp,compare.prcp.crtz,compare.prcp.lanl,compare.prcp.espn)
gdd.compare <- list(compare.gdd.mvnp,compare.gdd.crtz,compare.gdd.lanl,compare.gdd.espn)

prcp.values <- lapply(prcp.compare,FUN=function(x){paste(sprintf("%.3f", round(x[,"R2v"],digits=3))," (",sprintf("%.3f", round(x[,"N_RMSE"],digits=3)),")",sep='')})
gdd.values <- lapply(gdd.compare,FUN=function(x){paste(sprintf("%.3f", round(x[,"R2v"],digits=3))," (",sprintf("%.3f", round(x[,"N_RMSE"],digits=3)),")",sep='')})

all.values <- do.call('rbind',lapply(1:4,FUN=function(x){rbind(prcp.values[[x]],gdd.values[[x]])}))
rownames(all.values) <- c("PRCP","GDD","PRCP","GDD","PRCP","GDD","PRCP","GDD")
colnames(all.values) <- c("CAR$\\mathsf{^2}$ \\cite{Zuber2011,Zuber2013} \\\\ Shrinkage Regression \\\\ Min. RMSE","CAR$\\mathsf{^2}$ \\cite{Zuber2011,Zuber2013} \\\\ Shrinkage Regression \\\\ Min. AIC", "MCOR$\\mathsf{^2}$ \\\\ OLS Regression \\\\ Min. RMSE","MCOR$\\mathsf{^2}$ \\\\ OLS Regression \\\\ Min. AIC","PPR \\cite{Cook1999,Cook2004} \\\\ PC Regression \\\\ Min. RMSE","PPR \\cite{Cook1999,Cook2004} \\\\ PC Regression \\\\ Min. AIC")

all.values[1,order(prcp.compare[[1]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[1,order(prcp.compare[[1]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[1,order(prcp.compare[[1]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[1,order(prcp.compare[[1]][,"CV_RMSE"])[2]],"}}",sep='')
all.values[2,order(gdd.compare[[1]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[2,order(gdd.compare[[1]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[2,order(gdd.compare[[1]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[2,order(gdd.compare[[1]][,"CV_RMSE"])[2]],"}}",sep='')

all.values[3,order(prcp.compare[[2]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[3,order(prcp.compare[[2]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[3,order(prcp.compare[[2]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[3,order(prcp.compare[[2]][,"CV_RMSE"])[2]],"}}",sep='')
all.values[4,order(gdd.compare[[2]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[4,order(gdd.compare[[2]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[4,order(gdd.compare[[2]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[4,order(gdd.compare[[2]][,"CV_RMSE"])[2]],"}}",sep='')

all.values[5,order(prcp.compare[[3]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[5,order(prcp.compare[[3]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[5,order(prcp.compare[[3]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[5,order(prcp.compare[[3]][,"CV_RMSE"])[2]],"}}",sep='')
all.values[6,order(gdd.compare[[3]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[6,order(gdd.compare[[3]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[6,order(gdd.compare[[3]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[6,order(gdd.compare[[3]][,"CV_RMSE"])[2]],"}}",sep='')

all.values[7,order(prcp.compare[[4]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[7,order(prcp.compare[[4]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[7,order(prcp.compare[[4]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[7,order(prcp.compare[[4]][,"CV_RMSE"])[2]],"}}",sep='')
all.values[8,order(gdd.compare[[4]][,"N_RMSE"])[1]] <- paste("\\textcolor{BrickRed}{\\bf{",all.values[8,order(gdd.compare[[4]][,"CV_RMSE"])[1]],"}}",sep='')
all.values[8,order(gdd.compare[[4]][,"N_RMSE"])[2]] <- paste("\\textcolor{Cerulean}{\\bf{",all.values[8,order(gdd.compare[[4]][,"CV_RMSE"])[2]],"}}",sep='')


## Write mean absolute differences table to a text file
fileConn <- file("../TABLES/MODEL_COMPARE.tex")
writeLines(c("\\begin{sidewaystable}[!h]",
             "\\centering",
             "\\caption{\\textbf{Model comparison of different reconstruction techniques.} Values given are validation r-squared ($R_v^2$) and normalized root mean squared prediction error ($RMSE_n$; in parentheses). Red type indicates the technique that minimizes $RMSE_n$; blue type indicates the second-ranked model. Each of the six reconstruction techniques are characterized by a variable ranking method, a model classification, and a model selection method. Statistics reported are the mean of all cross-validated runs.}",
             "\\label{tab:MODEL_COMPARE}",
             "\\begin{tabular}{clcccccccc}",
             "\\toprule",
             "",
             "&",paste("& \\begin{sideways} \\parbox{4cm}{",colnames(all.values),"} \\end{sideways} ",sep=''),"\\bigstrut \\\\",
             "\\midrule",
             "&&$R_v^2$ ($RMSE_n$)&$R_v^2$ ($RMSE_n$)&$R_v^2$ ($RMSE_n$)&$R_v^2$ ($RMSE_n$)&$R_v^2$ ($RMSE_n$)&$R_v^2$ ($RMSE_n$)\\\\",
             "\\midrule",
             "",
             "\\multirow{2}[5]{*}{\\rotatebox{90}{MVNP}}","& PRCP", paste("& ",all.values[1,],sep=''),"\\bigstrut \\\\",
             "& GDD", paste("& ",all.values[2,],sep=''),"\\bigstrut \\\\",
             #             "\\\\",
             "\\multirow{2}[5]{*}{\\begin{sideways} CRTZ \\end{sideways}}","& PRCP", paste("& ",all.values[3,],sep=''),"\\bigstrut \\\\",
             "& GDD", paste("& ",all.values[4,],sep=''),"\\bigstrut \\\\",
             #             "\\\\",
             "\\multirow{2}[5]{*}{\\begin{sideways} LANL \\end{sideways}}","& PRCP", paste("& ",all.values[5,],sep=''),"\\bigstrut \\\\",
             "& GDD", paste("& ",all.values[6,],sep=''),"\\bigstrut \\\\",
             #             "\\\\",
             "\\multirow{2}[5]{*}{\\begin{sideways} ESPN \\end{sideways}}","& PRCP", paste("& ",all.values[7,],sep=''),"\\bigstrut \\\\",
             "& GDD", paste("& ",all.values[8,],sep=''),"\\bigstrut \\\\",
             "\\bottomrule",
             "\\end{tabular}",
             "\\end{sidewaystable}"),fileConn)
close(fileConn)


