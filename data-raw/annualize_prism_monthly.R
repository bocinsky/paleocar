annualize_prism_monthly <- function(prism.brick, months=c(1:12), fun){
  brick.names <- gsub("X","",names(prism.brick))
  brick.years <- as.numeric(substring(brick.names,1,4))
  brick.months <- as.numeric(substring(brick.names,5,6))
  
  # If more months than a year, break
  if(length(months)>12){
    stop("ERROR! Too many months.")
  }
  
  # Process the months to get months from current, previous, and future years
  previous.year <- 12+months[months<1]
  current.year <- months[months>=1 & months<=12]
  next.year <- months[months>12]-12
  no.year <- (1:12)[!(1:12 %in% c(previous.year,current.year,next.year))]
  
  brick.years[brick.months %in% previous.year] <- brick.years[brick.months %in% previous.year]+1
  brick.years[brick.months %in% next.year] <- brick.years[brick.months %in% next.year]-1
  brick.years[brick.months %in% no.year] <- 0
  
  prism.brick <- prism.brick[[which(brick.years %in% as.numeric(names(table(brick.years))[table(brick.years)==length(months)]))]]
  brick.years <- brick.years[which(brick.years %in% as.numeric(names(table(brick.years))[table(brick.years)==length(months)]))]
  
  fun.type <- fun
  
  cl <- parallel::makeCluster(8)
  out <- parallel::parLapply(cl, unique(brick.years), function(prism.brick, brick.years, fun.type, year){
    
    if(fun.type=="sum"){
      out <- raster::calc(raster::subset(prism.brick, subset=which(brick.years %in% year)), fun=sum)
    }else if(fun.type=="mean"){
      out <- raster::calc(raster::subset(prism.brick, subset=which(brick.years %in% year)), fun=mean)
    }
    return(out)
  }, prism.brick=prism.brick, brick.years=brick.years, fun.type=fun.type)
  
  parallel::stopCluster(cl)
  
  out <- raster::stack(out, quick=T)
  names(out) <- unique(brick.years)
  
  return(out)
}
