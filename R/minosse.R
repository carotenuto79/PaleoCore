#' @title A ready to use function to compute species geographic range with default settings
#' @description This function performs minosse.data and minosse.target all at once for one or multiple target species.
#' @usage minosse(dat,species.name=NULL,domain,coc_by="locality",min.occs=10,min.bkg=100,
#' sampling.by.distance=TRUE,n.sims=10,n.sims.clusters="automatic",prediction.ground=NULL,
#' abiotic.covs=NULL,combine.covs=FALSE,crop.by.mcp=FALSE,th_num=3,projection=NULL,
#' lon_0=NULL,lat_0=NULL,seed=NULL)
#' @param dat A n x m dataframe where n are the single occurrences and m are the following columns: spec (the species name), x and y (longitude and latitude in decimal degrees, respectively) and loc_id (an id identifying the fossil locality).
#' @param species.name The character vector of the species names whose geographic ranges are to be estimated. If NULL, minosse runs for all the species in the fossil dataset.
#' @param domain Character. Eithter "land", for terrestrial species, or "sea", for marine species.
#' @param coc_by Character. This argument enables the cooccurrence analysis performed either at the locality level (then use "locality") or at cell of the prediction ground level (then use "cell"). See details below.
#' @param min.occs Numeric. For both target and predictor species. The number occurrences below which minosse.data discards a species from being a valid target or elegible predictors. Default 10.
#' @param min.bkg Numeric. minosse function by default simulates as many pseudo abseces as the presences, thereby this is the minimum number of pseudo absences to simulate if a species occurrence number is below this value.
#' @param sampling.by.distance Logical. If TRUE pseudo absences are simulated with an intensity proportional to the distance to the presence data. If FALSE a pure spatial random distribution is simulated.
#' @param n.sims Numeric. The number of pseudo absences simulations (see details).
#' @param n.sims.clusters Numeric. The number of machine cores to use when setting multiple pseudo absences simulations. Default "automatic".
#' @param prediction.ground Either a raster or a SpatialPolygons class object where to perform all the spatial interpolations. This will be the prediction ground used when running minosse.target.
#' @param abiotic.covs the raster or rasters' stack of additional environmental predictors.
#' @param combine.covs Logical. Should minosse.data collate species and abiotic predictors when performing variables' number reduction? Default False. See details.
#' @param crop.by.mcp Logical. If TRUE, the interpoalation of the predictors species data are limited to the prediction.grund area delimited by the MCP of ALL the fossil occurrences. Default FALSE
#' @param th_num Numeric. The Regression Kriging prediction map binarization threshold value index as reported in the optimal.thresholds function in the PresenceAbsence package. Default 3 = MaxSens+Spec,	maximizes (sensitivity+specificity)/2.
#' @param projection Character. This argument works only if prediction.ground is NULL. This is the euqual-area projection for spatial interpolations. A character string in the proj4 format or either "moll" (Mollweide) or "laea" (Lambert Azimuthal equal area) projections (see details in minosse.data function).
#' @param lon_0 Numeric. Only if prediction.ground is NULL. The longitude of the projection centre used when setting either "moll" or "laea" projections.  If NULL the mean longitude of the whole fossil record is used. Default NULL.
#' @param lat_0 Numeric. Only if prediction.ground is NULL. The latitude of the projection centre used when setting "laea" projection.  If NULL the mean latitude of whole fossil record is used. Default NULL.
#' @param seed Numeric. The seed number for experiment replication.
#' @details NULL
#' @export
#' @return a list of three objects where the first one is the polygon of the target species geographic range (a SpatialPolygons object), the second element is the output of minosse.target function (see minosse.target function for details) and the last one is the result (if available) of the cooccurrence analysis.
#' If minosse function is performed for multiple species all at once, then minosse output described above is replicated for each target species.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \donttest{
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri", package="PaleoCore"))->prediction_ground
#'
#'   mam<-minosse(dat=lgm,species.name="Mammuthus_primigenius",domain="land",
#'   prediction.ground=prediction_ground,crop.by.mcp=FALSE,th_num=3,
#'   coc_by="locality",min.occs=3,min.bkg=100,sampling.by.distance=TRUE,
#'   n.sims=10,n.sims.clusters="automatic",projection=NULL,lon_0 = NULL,
#'   lat_0 = NULL,seed=625)
#'
#'   }


minosse<-function(dat,
                  species.name=NULL,
                  domain,
                  coc_by="locality",
                  min.occs=10,
                  min.bkg=100,
                  sampling.by.distance=TRUE,
                  n.sims=10,
                  n.sims.clusters="automatic",
                  prediction.ground=NULL,
                  abiotic.covs=NULL,
                  combine.covs=FALSE,
                  crop.by.mcp=FALSE,
                  th_num=3,
                  projection=NULL,
                  lon_0=NULL,
                  lat_0=NULL,
                  seed=NULL){

  if(is.null(species.name)) unique(as.character(dat$spec))->spec_list else strsplit(species.name," ")->spec_list
  if(is.null(projection)) stop("MInOSSE requires equal area-projected coordinates reference system to properly perform spatial analyses")
  minosse_dat<-list()
  minosse_res_list<-list()
  for(i in 1:length(spec_list)){
    if(class(try(minosse_dat[[i]]<-minosse.data(
      obj=dat,
      species_name=spec_list[[i]],
      domain=domain,
      coc.by=coc_by,
      min.occs=min.occs,
      abiotic.covs=abiotic.covs,
      combine.covs=combine.covs,
      reduce_covs_by="pca",
      covs_th=0.95,
      c.size="mean",
      bkg.predictors="presence",
      min.bkg=min.bkg,
      sampling.by.distance=sampling.by.distance,
      prediction.ground=prediction.ground,
      projection=projection,
      lon_0=lon_0,
      lat_0=lat_0,
      n.clusters="automatic",
      seed=seed)
    ))=="try-error") {next} else {minosse_dat[[i]]<-minosse_dat[[i]];

    if(class(try(minosse_res_list[[i]]<-minosse.target(
      resp=minosse_dat[[i]][[1]],
      predictors=minosse_dat[[i]][[2]],
      bkg="presence",
      min.bkg = min.bkg,
      n.sims=n.sims,
      sampling.by.distance=sampling.by.distance,
      n.folds=1,
      n.sims.clusters="automatic",
      seed=seed)))=="try-error") {next} else minosse_res_list[[i]]<-minosse_res_list[[i]]
    }
    print(i)
  }
  gc()

  minosse_res_list[which(!is.na(minosse_res_list))]->minosse_res_list
  spec_list[which(!is.na(minosse_res_list))]->spec_list
  minosse_dat[which(!is.na(minosse_res_list))]->minosse_dat
  lapply(minosse_res_list,function(x)minosse.poly(x,th_num))->minosse_poly_list
  minosse_res_list<-lapply(seq(1:length(spec_list)),function(x)list(minosse_poly=minosse_poly_list[[x]],minosse_res=minosse_res_list[[x]],sig_species=minosse_dat[[x]]$sig_species))
  names(minosse_res_list)<-spec_list
  if(length(minosse_res_list)==1) minosse_res_list[[1]]->minosse_res_list
  return(minosse_res_list)
}
