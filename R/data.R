#' MVNP Spatial Polygon
#'
#' A dataset containing the spatial polygon defining the boundary 
#' of Mesa Verde National Park in southwestern Colorado.
#'
#' @format SpatialPolygonsDataFrame
#' @source \url{http://nrdata.nps.gov/programs/Lands/meve_tracts.zip}
"mvnp"

#' Data from the International Tree Ring Data Bank
#'
#' A dataset containing the tree-ring data from the ITRDB for 10-degree buffer around 
#' Mesa Verde National Park.
#' 
#' These data were extracted from the ITRDB using the \code{FedData::get_itrdb()} function.
#'
#' @format list
#' \describe{
#'   \item{metadata}{A data.table or SpatialPointsDataFrame (if makeSpatial==TRUE) of the 
#'   locations and names of extracted ITRDB chrononlogies}
#'   \item{widths}{A matrix of tree-ring widths/densities given user selection}
#'   \item{depths}{A matrix of tree-ring sample depths}
#' }
#' @source \url{https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring}
"itrdb"

#' PRISM water_year precipitation
#'
#' A dataset containing the water-year (October--September) gridded precipitation estimates for 
#' Mesa Verde National Park derived from the 800m PRISM dataset.
#'
#' @format RasterStack
#' @source \url{http://www.prism.oregonstate.edu/}
"mvnp_prism"