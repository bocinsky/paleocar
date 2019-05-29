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
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NED/'.
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/NED/'.
#' @param prcp_threshold The minimum amount of water-year precipitation, in mm, required to be in the farming niche.
#' Defaults to 300 mm. Use 'NA' to suppress niche calculation.
#' @param gdd_threshold The minimum number of Fahrenheit growing degree daysrequired to be in the farming niche.
#' Defaults to 1800 FGDD. Use 'NA' to suppress niche calculation.
#' @param years An integer vector of years, between AD 1 and 2000, you wish to extract.
#' Defaults to 1:2000.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A logical RasterBrick of the niche, through time.
#' @importFrom magrittr %>% %<>%
#' @importFrom foreach %do%
#' @importFrom utils head tail
#' @importFrom purrr map
#' @export
get_bocinsky2016 <- function(template = NULL,
                             label = "paleocar",
                             raw.dir = "./RAW/PALEOCAR/",
                             extraction.dir = paste0("./EXTRACTIONS/", label, "/PALEOCAR/"),
                             prcp_threshold = 300,
                             gdd_threshold = 1800,
                             years = 1:2000,
                             force.redo = F){
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & file.exists(paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif")))
    return(raster::brick(paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif")))
  
  # Force Raster to load large rasters into memory
  raster::rasterOptions(chunksize=2e+07,maxmemory=2e+08)
  
  url <- 'http://www1.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/usa/bocinsky2016/'
  req <- httr::GET(url)
  files <- XML::readHTMLTable(rawToChar(req$content),stringsAsFactors = FALSE)[[1]]$Name # Get the file listing
  files <- files[grep("nc4",files)] # Only download netcdf files
  
  if(!is.null(template)){
    template <- sp::spTransform(FedData::polygon_from_extent(template),
                                sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
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
                           destdir = raw.dir)
  }
  
  # mosaick
  out_bricks <- foreach::foreach(type = c("PPT","GDD")) %do% {
    
    if(force.redo | !file.exists(paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"))){
      
      system(paste0("gdalbuildvrt temp.vrt ", paste0(raw.dir,"/",grep(type, files, value = T), collapse=" ")))
      system(paste0("gdal_translate -q -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND temp.vrt ",
                    paste0("'",extraction.dir,"/",type,"_merged.tif","'")))
      unlink("temp.vrt")
      
      tile_brick <- suppressWarnings(raster::brick(paste0(extraction.dir,"/",type,"_merged.tif")))
      
      if(!is.null(template)){
        tile_brick %<>%
          raster::crop(template) %>%
          raster::mask(template)
      }
      
      raster::writeRaster(tile_brick, paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"),
                          datatype="INT2U",
                          options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                          overwrite=T,
                          setStatistics=FALSE)
      
      unlink(paste0(extraction.dir,"/",type,"_merged.tif"))
      
    }
    
    return(
      raster::brick(paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"))
    )
  }
  
  names(out_bricks) <- c("PPT","GDD")
  
  out_bricks %<>% purrr::map(function(x){
    raster::projection(x) <- sp::CRS("+proj=longlat +datum=WGS84")
    return(x)
    }) 
  
  if(!is.na(prcp_threshold) & !is.na(gdd_threshold)){

    if(force.redo | !file.exists(paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"))){
      ppt_niche <- out_bricks$PPT >= prcp_threshold
        
        raster::writeRaster(ppt_niche, paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"),
                            datatype="INT1U",
                            options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                            overwrite=T,
                            setStatistics=FALSE)
      }
    ppt_niche <- raster::brick(paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"))
    

    if(force.redo | !file.exists(paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"))){
      gdd_niche <- out_bricks$GDD >= gdd_threshold
      
      raster::writeRaster(gdd_niche, paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"),
                          datatype="INT1U",
                          options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                          overwrite=T,
                          setStatistics=FALSE)
      
      
    }
    gdd_niche <- raster::brick(paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"))
    
    out <- ppt_niche * gdd_niche
    
    raster::writeRaster(out, paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif"),
                        datatype="INT1U",
                        options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                        overwrite=T,
                        setStatistics=FALSE)
    
    out_bricks$niche <- out
  }
  
  return(out_bricks)
}
