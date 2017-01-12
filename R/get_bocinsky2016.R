globalVariables(c("tile", "type"))
#' Download, mask, and calculate the maize growing niche from the Bocinsky 2016 data
#' available from the NOAA paleoclimate database:
#' 
#' Bocinsky, R. K., J. Rush, K. W. Kintigh, and T. A. Kohler. 2016. 
#' Exploration and exploitation in the macrohistory of the pre-Hispanic Pueblo Southwest. 
#' Science Advances 2:e1501532.
#'
#' @param template A Spatial* or Raster* object from which to define the area of interest.
#' If omitted, download entire reconstruction.
#' @param dir_out A directory in which to place the output.
#' @param prcp_threshold The minimum amount of water-year precipitation, in mm, required to be in the farming niche.
#' Defaults to 300 mm.
#' @param gdd_threshold The minimum number of Fahrenheit growing degree daysrequired to be in the farming niche.
#' Defaults to 1800 FGDD.
#' @param years An integer vector of years, between AD 1 and 2000, you wish to extract.
#' Defaults to 1:2000.

#' @return A logical RasterBrick of the niche, through time.
#' @importFrom magrittr %>% %<>%
#' @importFrom foreach %do%
#' @importFrom utils head tail
#' @export
get_bocinsky2016 <- function(template = NULL,
                               dir_out = "./OUTPUT/",
                               prcp_threshold = 300,
                               gdd_threshold = 1800,
                               years = 1:2000){
  
  dir.create(paste0(dir_out,"DATA/"), showWarnings = F, recursive = T)
  
  # Suppress scientific notation
  options(scipen=999999)
  # Force Raster to load large rasters into memory
  raster::rasterOptions(chunksize=2e+07,maxmemory=2e+08)
  
  url <- 'http://www1.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/usa/bocinsky2016/'
  req <- httr::GET(url)
  files <- XML::readHTMLTable(rawToChar(req$content),stringsAsFactors = FALSE)[[1]]$Name # Get the file listing
  files <- files[grep("nc4",files)] # Only download netcdf files
  
  if(!is.null(template)){
    template <- sp::spTransform(FedData::polygon_from_extent(template),sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
    extent.latlon <- raster::extent(template)
    
    # Open USGS NED download service.
    # NED tiles are labeled by their northwest corner. Thus, coordinate 36.42N, -105.71W is in grid n37w106
    wests <- seq(ceiling(abs(extent.latlon@xmax)),ceiling(abs(extent.latlon@xmin)))
    norths <- seq(floor(abs(extent.latlon@ymin)),floor(abs(extent.latlon@ymax)))
    
    tilesLocations <- as.matrix(expand.grid(norths,wests,stringsAsFactors = FALSE)) %>%
      apply(1,function(x){paste0(x[2],"W",x[1],"N")})
    
    
    message("Area of interest includes ",length(tilesLocations)," reconstruction tiles.")
    
    files <- sapply(tilesLocations, grep, x = files, value = T) %>% as.vector()
  }
  
  out <- foreach::foreach(tile = paste0(url,files)) %do% {
    FedData::download_data(url = tile,
                           destdir = paste0(dir_out,"DATA/"),
                           timestamping = FALSE)
  }
  
  foreach::foreach(type = c("PPT","GDD")) %do% {
    
    if(file.exists(paste0(dir_out,"/",type,"_niche_",head(years,1),"-",tail(years,1),".tif"))) return(NA)
    
    system(paste0("gdal_merge.py -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND -o ",paste0(dir_out,"/DATA/",type,"_merged.tif")," ", paste0("./OUTPUT/DATA/",grep(type,files,value = T),collapse=" ")))
    
    tile_brick <- raster::brick(paste0(dir_out,"/DATA/",type,"_merged.tif"))
    
    if(!is.null(template)){
      tile_brick %<>%
        raster::crop(template) %>%
        raster::mask(template)
    }
    
    raster::writeRaster(tile_brick, paste0(dir_out,"/",type,"_",head(years,1),"-",tail(years,1),".tif"),
                        datatype="INT2U",
                        options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                        overwrite=T,
                        setStatistics=FALSE)
    
    tile_niche <- tile_brick >= ifelse(type == "PPT", prcp_threshold, gdd_threshold)
    
    raster::writeRaster(tile_niche, paste0(dir_out,"/",type,"_niche_",head(years,1),"-",tail(years,1),".tif"),
                        datatype="INT1U",
                        options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                        overwrite=T,
                        setStatistics=FALSE)
    
    unlink(paste0(dir_out,"/DATA/",type,"_merged.tif"))
  }
  
  precip_niche <- raster::brick(paste0(dir_out,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"))
  gdd_niche <- raster::brick(paste0(dir_out,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"))
  
  out <- precip_niche * gdd_niche
  
  raster::writeRaster(out, paste0(dir_out,"/","niche_",head(years,1),"-",tail(years,1),".tif"),
                      datatype="INT1U",
                      options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                      overwrite=T,
                      setStatistics=FALSE)
  
  return(out)
}