######## BEGIN PREAMBLE ########

# Set the working directory to the location of R-scripts
# setwd("/Users/bocinsky/Desktop/Bocinsky_Kohler_2014/R/")
setwd("/Users/bocinsky/IMPORTANT/WSU/RESEARCH/MyPublications/SUBMITTED/BOCINSKY_KOHLER_2014/R")

# Load the functions for all analyses below
source('./UTILITY_FUNCTIONS.R')
source('./PACKAGES.R')
source('./GHCN_FUNCTIONS.R')
source('./PRISM_FUNCTIONS.R')
source('./ITRDB_FUNCTIONS.R')
source('./CROSS_VALIDATE_FUNCTIONS.R')
source('./COMPARE_METHODS_FUNCTIONS.R')
source('./COOK_PCR_FUNCTIONS.R')
source('./SLM_CUSTOM_FUNCTIONS.R')

# Suppress use of scientific notation
options(scipen=999)

# Create an output directory above the R script directory
dir.create("../OUTPUT/", showWarnings=F)

######## END PREAMBLE ########


## [INTENTIONALLY LEFT BLANK] ##


######## BEGIN PARAMETERS ########
## Set the master data directory
## This should be the location of PRISM data
# MASTER.DATA <- "/Users/bocinsky/Desktop/DATA/"
MASTER.DATA <- "/Volumes/BOCINSKY_DATA/DATA/"


## Set the calibration period
# Here, I use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983

# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)

## Define the VEPIIN and VEPIIS study areas.
VEPIIN.polygon <- createArea(North = 4170000, South = 4102000, East = 740000, West = 672800, projection.string = "+proj=utm +datum=NAD83 +zone=12")
VEPIIS.polygon <- createArea(North = 4030400, South = 3939600, East = 435800, West = 359200, projection.string = "+proj=utm +datum=NAD83 +zone=13")

######## END PARAMETERS ########


## [INTENTIONALLY LEFT BLANK] ##


######## BEGIN DATA IMPORT AND PREPROCESSING ########
## Load the PRISM (interpolated climate) data for eash of the study regions
# VEPIIN GDD
VEPIIN.PRISM.monthly.tmin <- getPRISM_MONTHLYData(template=VEPIIN.polygon, type='tmin', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIN", force.redo=F)
VEPIIN.PRISM.monthly.tmax <- getPRISM_MONTHLYData(template=VEPIIN.polygon, type='tmax', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIN", force.redo=F)
VEPIIN.PRISM.grow.GDD <- calcANNUALGDD_MONTHLY(extraction.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/VEPIIN/", sep=''), months=c(5,6,7,8,9), t.base=10, t.cap=30)
VEPIIN.PRISM.grow.GDD <- subset(VEPIIN.PRISM.grow.GDD,grep(paste(calibration.years,collapse="|"),names(VEPIIN.PRISM.grow.GDD)))

# VEPIIN water-year prcp
VEPIIN.PRISM.monthly.prcp <- getPRISM_MONTHLYData(template=VEPIIN.polygon, type='ppt', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIN", force.redo=F)
VEPIIN.PRISM.h20_year.prcp <- annualizePRISM_MONTHLY(prism.brick=VEPIIN.PRISM.monthly.prcp, months=c(-2,-1,0,1,2,3,4,5,6,7,8,9), fun="sum")
VEPIIN.PRISM.h20_year.prcp <- subset(VEPIIN.PRISM.h20_year.prcp,grep(paste(calibration.years,collapse="|"),names(VEPIIN.PRISM.h20_year.prcp)))

# Clean-up
rm(VEPIIN.PRISM.monthly.tmin,VEPIIN.PRISM.monthly.tmax,VEPIIN.PRISM.monthly.prcp)
gc()
gc()

# VEPIIS GDD
VEPIIS.PRISM.monthly.tmin <- getPRISM_MONTHLYData(template=VEPIIS.polygon, type='tmin', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIS", force.redo=F)
VEPIIS.PRISM.monthly.tmax <- getPRISM_MONTHLYData(template=VEPIIS.polygon, type='tmax', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIS", force.redo=F)
VEPIIS.PRISM.grow.GDD <- calcANNUALGDD_MONTHLY(extraction.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/VEPIIS/", sep=''), months=c(5,6,7,8,9), t.base=10, t.cap=30)
VEPIIS.PRISM.grow.GDD <- subset(VEPIIS.PRISM.grow.GDD,grep(paste(calibration.years,collapse="|"),names(VEPIIS.PRISM.grow.GDD)))

# VEPIIS water-year prcp
VEPIIS.PRISM.monthly.prcp <- getPRISM_MONTHLYData(template=VEPIIS.polygon, type='ppt', out.dir=paste(MASTER.DATA,"PRISM/EXTRACTIONS/", sep=''), monthly.dir=paste(MASTER.DATA,"PRISM/LT81_800M/",sep=''), label="VEPIIS", force.redo=F)
VEPIIS.PRISM.h20_year.prcp <- annualizePRISM_MONTHLY(prism.brick=VEPIIS.PRISM.monthly.prcp, months=c(-2,-1,0,1,2,3,4,5,6,7,8,9), fun="sum")
VEPIIS.PRISM.h20_year.prcp <- subset(VEPIIS.PRISM.h20_year.prcp,grep(paste(calibration.years,collapse="|"),names(VEPIIS.PRISM.h20_year.prcp)))

# Clean-up
rm(VEPIIS.PRISM.monthly.tmin,VEPIIS.PRISM.monthly.tmax,VEPIIS.PRISM.monthly.prcp)
gc()
gc()


## Load the ITRDB (tree ring) data, and crop to study period
if(file.exists(paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_DATA.csv", sep=''))){
  # Get a SPDF of ITRDB metadata and get ITRDB database
  ITRDB.data <- read.csv(paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_DATA.csv", sep=''))
  ITRDB.meta.sp <- getITRDBMetadataSPDF(names(ITRDB.data),data.dir=paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/", sep=''))
}else{
  ITRDB.data <- getITRDB(raw.dir = paste(MASTER.DATA,"DENDRO/ITRDB/",sep=''), output.dir=paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/", sep=''), type='standard', download=F, force.redo=T)
  ITRDB.meta <- read.csv(paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_METADATA.csv", sep=''),colClasses="character")
  
  # Get a SPDF of ITRDB metadata
  ITRDB.meta.sp <- getITRDBMetadataSPDF(names(ITRDB.data),data.dir=paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/", sep=''))
  
  # Get a shapefile of the 48 CONUS states
  states <- readOGR(paste(MASTER.DATA,"NATIONAL_ATLAS/statep010",sep=''), layer='statep010')
  # And select only the four corners
  states <- states[states$STATE %in% c("Arizona","Colorado","Utah","New Mexico"),]
  states <- spTransform(states,CRS(projection(ITRDB.meta.sp)))
  
  # Trim the IRTB database to the four corners states
  ITRDB.meta.sp <- ITRDB.meta.sp[as.vector(!is.na((ITRDB.meta.sp %over% states)[,1])),]
  # Select only the four corners series
  ITRDB.data <- ITRDB.data[,c("YEAR",as.character(ITRDB.meta.sp$SERIES))]
  # Remove years with no data
  ITRDB.data <- ITRDB.data[apply(ITRDB.data,1,FUN=function(i){any(!is.na(i[-1]))}),]
  # Clean up the chronology names
  ITRDB.meta.sp$NAME <- sanitizeITRDBnames(ITRDB.meta.sp$NAME)
  
  # Write the amended ITRDB data and metadata
  write.csv(ITRDB.data,paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_DATA.csv", sep=''),row.names=F)
  write.csv(ITRDB.meta.sp,paste(MASTER.DATA,"DENDRO/ITRDB/EXTRACTIONS/ITRDB_METADATA.csv", sep=''), row.names=F)
}

######## FINISH DATA IMPORT AND PREPROCESSING ########


## [INTENTIONALLY LEFT BLANK] ##


######## BEGIN ANALYSIS ########
#### Spatial retrodictions of all climate signals ####
## Isolating the calibration years for slm coefficient estimation
training.series <- ITRDB.data[ITRDB.data$YEAR %in% calibration.years,]
training.series <- training.series[,t(complete.cases(t(training.series)))]
training.series <- training.series[,-1]

## Get list of periods in retrodiction years with stable chronologies, and list chronologies
retro.series <- ITRDB.data[ITRDB.data$YEAR %in% 1:2000,]
retro.series <- retro.series[,c("YEAR",names(training.series))]
rownames(retro.series) <- retro.series$YEAR
retro.series <- retro.series[,-1]
retro.series[!is.na(retro.series)] <- TRUE
retro.series[is.na(retro.series)] <- FALSE
retro.series <- unique(retro.series)
predlist <- lapply(1:nrow(retro.series),function(i){which(retro.series[i,]==1)})
names(predlist) <- rownames(retro.series)
retro.series.count <- rowSums(retro.series)

recon.series <- ITRDB.data[ITRDB.data$YEAR %in% 1:2000,]
recon.series <- recon.series[,c("YEAR",names(training.series))]
rownames(recon.series) <- recon.series$YEAR
breaks <- c(as.numeric(rownames(retro.series)),2001)
breaks.lengths <- diff(breaks)
retro.series.list <- lapply(1:(length(breaks)-1),function(i){recon.series[as.character(breaks[i]):as.character(breaks[i+1]-1),-1,drop=F]})

VEPIIN.PRISM.h20_year.prcp.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIN.PRISM.h20_year.prcp, predlist=predlist, name="VEPIIN_PRISM_h20_year_prcp")
VEPIIN.PRISM.h20_year.prcp.recons <- recon(name="VEPIIN_PRISM_h20_year_prcp_unique", models=VEPIIN.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIN.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)
rm(VEPIIN.PRISM.h20_year.prcp.recons); gc(); gc()
VEPIIN.PRISM.h20_year.prcp.errors <- recon.errors(name="VEPIIN_PRISM_h20_year_prcp_unique", models=VEPIIN.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIN.PRISM.h20_year.prcp, recon.series=recon.series)
VEPIIN.PRISM.h20_year.prcp.lm.models.union <- slm.models.custom.brick.union(models=VEPIIN.PRISM.h20_year.prcp.lm.models, X=training.series, Y.brick=VEPIIN.PRISM.h20_year.prcp, name="VEPIIN_PRISM_h20_year_prcp")
rm(VEPIIN.PRISM.h20_year.prcp.lm.models); gc(); gc()
VEPIIN.PRISM.h20_year.prcp.recons.union <- recon(name="VEPIIN_PRISM_h20_year_prcp_union", models=VEPIIN.PRISM.h20_year.prcp.lm.models.union, Ytrain.brick=VEPIIN.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)
rm(VEPIIN.PRISM.h20_year.prcp.recons.union); gc(); gc()
VEPIIN.PRISM.h20_year.prcp.errors <- recon.errors(name="VEPIIN_PRISM_h20_year_prcp_union", models=VEPIIN.PRISM.h20_year.prcp.lm.models.union, Ytrain.brick=VEPIIN.PRISM.h20_year.prcp, recon.series=recon.series)
rm(VEPIIN.PRISM.h20_year.prcp.lm.models.union); gc(); gc()

VEPIIN.PRISM.grow.GDD.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIN.PRISM.grow.GDD, predlist=predlist, name="VEPIIN_PRISM_grow_GDD")
VEPIIN.PRISM.grow.GDD.recons <- recon(name="VEPIIN_PRISM_grow_GDD_unique", models=VEPIIN.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIN.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)
rm(VEPIIN.PRISM.grow.GDD.recons); gc(); gc()
VEPIIN.PRISM.grow.GDD.errors <- recon.errors(name="VEPIIN_PRISM_grow_GDD_unique", models=VEPIIN.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIN.PRISM.grow.GDD, recon.series=recon.series)
VEPIIN.PRISM.grow.GDD.lm.models.union <- slm.models.custom.brick.union(models=VEPIIN.PRISM.grow.GDD.lm.models, X=training.series, Y.brick=VEPIIN.PRISM.grow.GDD, name="VEPIIN_PRISM_grow_GDD")
rm(VEPIIN.PRISM.grow.GDD.lm.models); gc(); gc()
VEPIIN.PRISM.grow.GDD.recons.union <- recon(name="VEPIIN_PRISM_grow_GDD_union", models=VEPIIN.PRISM.grow.GDD.lm.models.union, Ytrain.brick=VEPIIN.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)
rm(VEPIIN.PRISM.grow.GDD.recons.union); gc(); gc()
VEPIIN.PRISM.grow.GDD.errors <- recon.errors(name="VEPIIN_PRISM_grow_GDD_union", models=VEPIIN.PRISM.grow.GDD.lm.models.union, Ytrain.brick=VEPIIN.PRISM.grow.GDD, recon.series=recon.series)
rm(VEPIIN.PRISM.grow.GDD.lm.models.union); gc(); gc()

VEPIIS.PRISM.h20_year.prcp.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIS.PRISM.h20_year.prcp, predlist=predlist, name="VEPIIS_PRISM_h20_year_prcp")
VEPIIS.PRISM.h20_year.prcp.recons <- recon(name="VEPIIS_PRISM_h20_year_prcp_unique", models=VEPIIS.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIS.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)
rm(VEPIIS.PRISM.h20_year.prcp.recons); gc(); gc()
VEPIIS.PRISM.h20_year.prcp.errors <- recon.errors(name="VEPIIS_PRISM_h20_year_prcp_unique", models=VEPIIS.PRISM.h20_year.prcp.lm.models, Ytrain.brick=VEPIIS.PRISM.h20_year.prcp, recon.series=recon.series)
VEPIIS.PRISM.h20_year.prcp.lm.models.union <- slm.models.custom.brick.union(models=VEPIIS.PRISM.h20_year.prcp.lm.models, X=training.series, Y.brick=VEPIIS.PRISM.h20_year.prcp, name="VEPIIS_PRISM_h20_year_prcp")
rm(VEPIIS.PRISM.h20_year.prcp.lm.models); gc(); gc()
VEPIIS.PRISM.h20_year.prcp.recons.union <- recon(name="VEPIIS_PRISM_h20_year_prcp_union", models=VEPIIS.PRISM.h20_year.prcp.lm.models.union, Ytrain.brick=VEPIIS.PRISM.h20_year.prcp, recon.series=recon.series, training.series=training.series)
rm(VEPIIS.PRISM.h20_year.prcp.recons.union); gc(); gc()
VEPIIS.PRISM.h20_year.prcp.errors <- recon.errors(name="VEPIIS_PRISM_h20_year_prcp_union", models=VEPIIS.PRISM.h20_year.prcp.lm.models.union, Ytrain.brick=VEPIIS.PRISM.h20_year.prcp, recon.series=recon.series)
rm(VEPIIS.PRISM.h20_year.prcp.lm.models.union); gc(); gc()

VEPIIS.PRISM.grow.GDD.lm.models <- slm.models.custom.brick.unique(X=training.series, Y.brick=VEPIIS.PRISM.grow.GDD, predlist=predlist, name="VEPIIS_PRISM_grow_GDD")
VEPIIS.PRISM.grow.GDD.recons <- recon(name="VEPIIS_PRISM_grow_GDD_unique", models=VEPIIS.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIS.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)
rm(VEPIIS.PRISM.grow.GDD.recons); gc(); gc()
VEPIIS.PRISM.grow.GDD.errors <- recon.errors(name="VEPIIS_PRISM_grow_GDD_unique", models=VEPIIS.PRISM.grow.GDD.lm.models, Ytrain.brick=VEPIIS.PRISM.grow.GDD, recon.series=recon.series)
VEPIIS.PRISM.grow.GDD.lm.models.union <- slm.models.custom.brick.union(models=VEPIIS.PRISM.grow.GDD.lm.models, X=training.series, Y.brick=VEPIIS.PRISM.grow.GDD, name="VEPIIS_PRISM_grow_GDD")
rm(VEPIIS.PRISM.grow.GDD.lm.models); gc(); gc()
VEPIIS.PRISM.grow.GDD.recons.union <- recon(name="VEPIIS_PRISM_grow_GDD_union", models=VEPIIS.PRISM.grow.GDD.lm.models.union, Ytrain.brick=VEPIIS.PRISM.grow.GDD, recon.series=recon.series, training.series=training.series)
rm(VEPIIS.PRISM.grow.GDD.recons.union); gc(); gc()
VEPIIS.PRISM.grow.GDD.errors <- recon.errors(name="VEPIIS_PRISM_grow_GDD_union", models=VEPIIS.PRISM.grow.GDD.lm.models.union, Ytrain.brick=VEPIIS.PRISM.grow.GDD, recon.series=recon.series)
rm(VEPIIS.PRISM.grow.GDD.lm.models.union); gc(); gc()

## Calculate the PRCP, GDD, and Total maize niche
source("./NICHE_FUNCTIONS.R")

## Do reconstructions and model comparison for the four GHCN station locations
source("./GHCN_FIGURES_AND_TABLES.R")

## Generate other figures for publication
source("../FIGURES/R/SW_MAP.R")
source("../FIGURES/R/TREE_AVAILABILITY.R")
