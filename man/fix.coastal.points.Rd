% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix.coastal.points.R
\name{fix.coastal.points}
\alias{fix.coastal.points}
\title{Shift points to the nearest non-NA cell of a base raster within a maximum distance}
\usage{
fix.coastal.points(data, xy.cols, r, ncell, clip.by.poly = F, poly,
  occ.desaggregation = F)
}
\arguments{
\item{data}{a matrix or data.frame containing the points coordinates}

\item{xy.cols}{a numeric vector indicating the IDs of the x-y columns}

\item{r}{a \code{RasterLayer} object to be used as base raster map}

\item{ncell}{an integer indicating the maximum number of cells away from the points where they can be shifted. If the nearest non-NA cell fall farther than this distance from a give point, this will be dropped}

\item{clip.by.poly}{a logical vector indicating if the base raster must be clipped along the boundaries of a polygon}

\item{poly}{a \code{SpatialPolygons} object used to clip the base raster map}

\item{occ.desaggregation}{a logical vector indicating if to apply the \code{occ.desaggregation.RASTER} function}
}
\value{
a matrix including the shifted points
}
\description{
The function moves points to the nearest non-NA cell of a base raster within a maximum distance.
}
\examples{
\dontrun{
library(raster)
raster(system.file("exdata/prediction_ground.gri", package="PaleoCore"))->prediction_ground
data(lgm)

data_clupus<-subset(lgm, spec=="Canis_lupus")

coordinates(data_clupus)<-~x+y
proj4string(data_clupus)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
data_clupus<-spTransform(data_clupus, proj4string(prediction_ground))
data_clupus<-as.data.frame(data_clupus)

data_mammuth_reduced<-fix.coastal.points(data_clupus, xy.cols = 4:5,
r = prediction_ground, ncell = 1, occ.desaggregation = F)
}
}
\author{
Mirko Di Febbraro
}
