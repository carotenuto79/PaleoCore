% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minosse.R
\name{minosse}
\alias{minosse}
\title{A ready to use function to compute species geographic range with default settings}
\usage{
minosse(dat,species.name=NULL,domain,coc_by="locality",min.occs=10,min.bkg=100,
sampling.by.distance=TRUE,n.sims=10,n.sims.clusters="automatic",prediction.ground=NULL,
abiotic.covs=NULL,combine.covs=FALSE,crop.by.mcp=FALSE,th_num=3,projection=NULL,
lon_0=NULL,lat_0=NULL,seed=NULL)
}
\arguments{
\item{dat}{A n x m dataframe where n are the single occurrences and m are the following columns: spec (the species name), x and y (longitude and latitude in decimal degrees, respectively) and loc_id (an id identifying the fossil locality).}

\item{species.name}{The character vector of the species names whose geographic ranges are to be estimated. If NULL, minosse runs for all the species in the fossil dataset.}

\item{domain}{Character. Eithter "land", for terrestrial species, or "sea", for marine species.}

\item{coc_by}{Character. This argument enables the cooccurrence analysis performed either at the locality level (then use "locality") or at cell of the prediction ground level (then use "cell"). See details below.}

\item{min.occs}{Numeric. For both target and predictor species. The number occurrences below which minosse.data discards a species from being a valid target or elegible predictors. Default 10.}

\item{min.bkg}{Numeric. minosse function by default simulates as many pseudo abseces as the presences, thereby this is the minimum number of pseudo absences to simulate if a species occurrence number is below this value.}

\item{sampling.by.distance}{Logical. If TRUE pseudo absences are simulated with an intensity proportional to the distance to the presence data. If FALSE a pure spatial random distribution is simulated.}

\item{n.sims}{Numeric. The number of pseudo absences simulations (see details).}

\item{n.sims.clusters}{Numeric. The number of machine cores to use when setting multiple pseudo absences simulations. Default "automatic".}

\item{prediction.ground}{Either a raster or a SpatialPolygons class object where to perform all the spatial interpolations. This will be the prediction ground used when running minosse.target.}

\item{abiotic.covs}{the raster or rasters' stack of additional environmental predictors.}

\item{combine.covs}{Logical. Should minosse.data collate species and abiotic predictors when performing variables' number reduction? Default False. See details.}

\item{crop.by.mcp}{Logical. If TRUE, the interpoalation of the predictors species data are limited to the prediction.grund area delimited by the MCP of ALL the fossil occurrences. Default FALSE}

\item{th_num}{Numeric. The Regression Kriging prediction map binarization threshold value index as reported in the optimal.thresholds function in the PresenceAbsence package. Default 3 = MaxSens+Spec,    maximizes (sensitivity+specificity)/2.}

\item{projection}{Character. This argument works only if prediction.ground is NULL. This is the euqual-area projection for spatial interpolations. A character string in the proj4 format or either "moll" (Mollweide) or "laea" (Lambert Azimuthal equal area) projections (see details in minosse.data function).}

\item{lon_0}{Numeric. Only if prediction.ground is NULL. The longitude of the projection centre used when setting either "moll" or "laea" projections.  If NULL the mean longitude of the whole fossil record is used. Default NULL.}

\item{lat_0}{Numeric. Only if prediction.ground is NULL. The latitude of the projection centre used when setting "laea" projection.  If NULL the mean latitude of whole fossil record is used. Default NULL.}

\item{seed}{Numeric. The seed number for experiment replication.}
}
\value{
a list of three objects where the first one is the polygon of the target species geographic range (a SpatialPolygons object), the second element is the output of minosse.target function (see minosse.target function for details) and the last one is the result (if available) of the cooccurrence analysis.
If minosse function is performed for multiple species all at once, then minosse output described above is replicated for each target species.
}
\description{
This function performs minosse.data and minosse.target all at once for one or multiple target species.
}
\examples{
  \donttest{
  library(raster)
  data(lgm)
  raster(system.file("exdata/prediction_ground.gri", package="PaleoCore"))->prediction_ground

  mam<-minosse(dat=lgm,species.name="Mammuthus_primigenius",domain="land",
  prediction.ground=prediction_ground,crop.by.mcp=FALSE,th_num=3,
  coc_by="locality",min.occs=3,min.bkg=100,sampling.by.distance=TRUE,
  n.sims=10,n.sims.clusters="automatic",projection=NULL,lon_0 = NULL,
  lat_0 = NULL,seed=625)

  }
}
\author{
Francesco Carotenuto, francesco.carotenuto@unina.it
}
