## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# this avoids console spam below
library(terra)
library(units)
library(sf)
library(sp)

## ----setup--------------------------------------------------------------------
library(snapKrig)

# used in examples
library(sp)
library(terra)
library(sf)
library(units)

## ----intro_id-----------------------------------------------------------------
# define a matrix and pass it to sk to get a grid object
mat_id = diag(10)
g = sk(mat_id)

# report info about the object on console
print(g)
class(g)

## ----intro_id_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'-------
# plot the grid with matrix theme
plot(g, ij=TRUE, main='identity matrix')

## ----intro_id_logi_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# make a grid of logical values
g0 = g == 0
print(g0)

# plot 
plot(g0, ij=TRUE, col_grid='white', main='matrix off-diagonals')

## ----intro_varmat_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# get a covariance matrix for 10 by 10 grid
vmat = sk_var(sk(10))

# plot the grid in matrix style
plot(sk(vmat), ij=TRUE, main='a covariance matrix')

## ----intro_vectorization------------------------------------------------------
# extract vectorization
vec = g[]

# compare to R's vectorization
all.equal(vec, c(mat_id))

## ----intro_sim_plot, out.width='75%', fig.dim=c(9, 5), fig.align='center'-----
# simulate data on a rectangular grid
gdim = c(100, 200)
g_sim = sk_sim(gdim)

# plot the grid in raster style
plot(g_sim, main='an example snapKrig simulation', cex=1.5)

## ----intro_sim_nonug_plot, out.width='75%', fig.dim=c(9, 5), fig.align='center'----
# get the default covariance parameters and modify nugget
pars = sk_pars(gdim)
pars[['eps']] = 1e-6

# simulate data on a rectangular grid
g_sim = sk_sim(gdim, pars)

# plot the result
plot(g_sim, main='an example snapKrig simulation (near-zero nugget)', cex=1.5)

## ----intro_sim_pars_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# plot the covariance footprint
sk_plot_pars(pars)

## ----intro_summary_method-----------------------------------------------------
summary(g_sim)

## ----intro_export-------------------------------------------------------------
sk_export(g_sim)

## ----intro_sim_upscaled_plot, out.width='75%', fig.dim=c(9, 5), fig.align='center'----
# upscale
g_sim_up = sk_rescale(g_sim, up=4)

# plot result
plot(g_sim_up, main='simulation data up-scaled by factor 4X', cex=1.5)

## ----intro_sim_downscaled_plot, out.width='75%', fig.dim=c(9, 5), fig.align='center'----
# downscale
g_sim_down = sk_rescale(g_sim_up, down=4)

# plot result
plot(g_sim_down, main='up-scaled by factor 4X then down-scaled by factor 4X', cex=1.5)

## ----intro_sim_downscaled_pred_plot, out.width='75%', fig.dim=c(9, 5), fig.align='center'----
# upscale
g_sim_down_pred = sk_cmean(g_sim_down, pars)

# plot result
plot(g_sim_down_pred, main='down-scaled values imputed by snapKrig', cex=1.5)

## ----intro_LL_bounds----------------------------------------------------------
# pick two model parameters for illustration
p_nm = stats::setNames(c('y.rho', 'x.rho'), c('y range', 'x range'))

# set bounds for two parameters and define test parameters
n_test = 25
bds = sk_bds(pars, g_sim_up)[p_nm, c('lower', 'upper')]
bds_test = list(y=seq(bds['y.rho', 1], bds['y.rho', 2], length.out=n_test),
                x=seq(bds['x.rho', 1], bds['x.rho', 2], length.out=n_test))

## ----intro_LL_surface---------------------------------------------------------
# make a grid of test parameters
g_test = sk(gyx=bds_test)
p_all = sk_coords(g_test)

# fill in the grid with log-likelihood values
for(i in seq_along(g_test))
{
  # modify the model parameters with test values
  p_test = sk_pars_update(pars)
  p_test[p_nm] = p_all[i,]
  
  # compute likelihood and copy to grid
  g_test[i] = sk_LL(sk_pars_update(pars, p_test), g_sim_up)
}

## ----intro_LL_surface_plot, out.width='50%', fig.dim=c(5, 5), fig.align='center'----
# plot the likelihood surface
plot(g_test, asp=2, main='log-likelihood surface', ylab=names(p_nm)[1], xlab=names(p_nm)[2], reset=FALSE)

# highlight the MLE
i_best = which.max(g_test[])
points(p_all[i_best,'x'], p_all[i_best,'y'], col='white', cex=1.5, lwd=1.5)

## ----intro_LL_surface_expected------------------------------------------------
# print the true values
print(c(x=pars[['x']][['kp']][['rho']], y=pars[['y']][['kp']][['rho']]))

## ----meuse_helper, include=FALSE----------------------------------------------

# load the Meuse data into a convenient format
get_meuse = function(dfMaxLength = units::set_units(50, m))
{
  # Note: dfMaxLength sets the interval used to sample line geometries of the river
  # using Voronoi tiles. This is a fussy and not well-tested algorithm for finding the
  # centre line of a river polygon, but it seems to work well enough for the example here

  # EPSG code for the coordinate system
  epsg_meuse = 28992

  # open river location data
  utils::data(meuse.riv)
  crs_meuse = sf::st_crs(epsg_meuse)[['wkt']]

  # reshape the river (edge) point data as a more densely segmented polygon
  colnames(meuse.riv) = c('x', 'y')
  meuse_river_points = sf::st_as_sf(as.data.frame(meuse.riv), coords=c('x', 'y'), crs=crs_meuse)
  meuse_river_seg = sf::st_cast(sf::st_combine(meuse_river_points), 'LINESTRING')
  meuse_river_poly = sf::st_cast(st_segmentize(meuse_river_seg, dfMaxLength), 'POLYGON')

  # skeletonization trick to get a single linestring at center of the river
  meuse_river_voronoi = sf::st_cast(sf::st_voronoi(meuse_river_poly, bOnlyEdges=TRUE), 'POINT')
  meuse_river_skele = sf::st_intersection(meuse_river_voronoi, meuse_river_poly)
  n_skele = length(meuse_river_skele)

  # compute distance matrix
  dmat_skele = units::drop_units(sf::st_distance(meuse_river_skele))

  # re-order to start from northernmost point
  idx_first = which.max(st_coordinates(meuse_river_skele)[,2])
  idx_reorder = c(idx_first, integer(n_skele-1L))
  for(idx_skele in seq(n_skele-1L))
  {
    # find least distance match
    idx_tocheck = seq(n_skele) != idx_first
    idx_first = which(idx_tocheck)[ which.min(dmat_skele[idx_tocheck, idx_first]) ]
    idx_reorder[1L+idx_skele] = idx_first

    # modify distance matrix so the matching point is not selected again
    dmat_skele[idx_first, ] = Inf
  }

  # connect the points to get the spine
  meuse_river = sf::st_cast(sf::st_combine(meuse_river_skele[idx_reorder]), 'LINESTRING')

  # load soil points data
  utils::data(meuse)
  meuse_soils = sf::st_as_sf(meuse, coords=c('x', 'y'), crs=epsg_meuse)

  # add 'distance' (to river) and 'logzinc' columns
  meuse_soils[['distance']] = units::drop_units( sf::st_distance(meuse_soils, meuse_river))
  meuse_soils[['log_zinc']] = log(meuse_soils[['zinc']])

  # crop the river objects to buffered bounding box of soils data
  bbox_padded = st_buffer(sf::st_as_sfc(sf::st_bbox(meuse_soils)), units::set_units(500, m))
  meuse_river_poly = sf::st_crop(meuse_river_poly, bbox_padded)
  meuse_river = sf::st_crop(meuse_river, bbox_padded)

  # return three geometry objects in a list
  return( list(soils=meuse_soils, river_poly=meuse_river_poly, river_line=meuse_river) )
}

## ----meuse_load---------------------------------------------------------------
# load the Meuse data into a convenient format
meuse_sf = get_meuse()

# extract the logarithm of the zinc concentration as sf points
pts = meuse_sf[['soils']]['log_zinc']

## ----meuse_source_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# set up a common color palette (this is the default in snapKrig)
.pal = function(n) { hcl.colors(n, 'Spectral', rev=TRUE) }

# plot source data using sf package
plot(pts, pch=16, reset=FALSE, pal=.pal, key.pos=1, main='Meuse log[zinc]')
plot(meuse_sf[['river_poly']], col='lightblue', border=NA, add=TRUE)
plot(st_geometry(pts), pch=1, add=TRUE)

## ----meuse_snap_default-------------------------------------------------------
# snap points with default settings
g = sk_snap(pts)
print(g)

## ----meuse_snap---------------------------------------------------------------
# snap again to 50m x 50m grid
g = sk_snap(pts, list(gres=c(50, 50)))
print(g)
summary(g)

## ----meuse_snapped_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# plot gridded version using the snapKrig package
plot(g, zlab='log(ppb)', main='snapped Meuse log[zinc] data')
plot(meuse_sf[['river_poly']], col='lightblue', border=NA, add=TRUE)


## ----meuse_make_river_dist----------------------------------------------------
# measure distances for every point in the grid
river_dist = sf::st_distance(sk_coords(g, out='sf'), meuse_sf[['river_line']])

## ----meuse_make_x-------------------------------------------------------------
# make a copy of g and insert the scaled distances as grid point values
X = g
X[] = scale( as.vector( units::drop_units(river_dist) ) )
summary(X)

## ----meuse_make_x_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
# plot the result
plot(X, zlab='distance\n(scaled)', main='distance to river covariate')
plot(meuse_sf[['river_line']], add=TRUE)

## ----meuse_fit_uk-------------------------------------------------------------
#fit the covariance model and trend with X
fit_result_uk = sk_fit(g, X=X, quiet=TRUE)

## ----meuse_fit_uk_likelihood--------------------------------------------------
# compute model likelihood
sk_LL(fit_result_uk, g, X)

## ----meuse_pred_uk------------------------------------------------------------
# compute conditional mean and variance
g_uk = sk_cmean(g, fit_result_uk, X)

## ----meuse_pred_uk_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
plot(g_uk, zlab='log[zinc]', main='universal kriging predictions')
plot(meuse_sf[['river_line']], add=TRUE)

## ----meuse_var_uk-------------------------------------------------------------
# compute conditional mean and variance
g_uk_var = sk_cmean(g, fit_result_uk, X, what='v', quiet=TRUE)

## ----meuse_var_uk_plot, out.width='50%', fig.dim=c(5,5), fig.align='center'----
plot(sqrt(g_uk_var), zlab='log[zinc]', main='universal kriging standard error')
plot(meuse_sf[['river_line']], add=TRUE)
plot(st_geometry(pts), pch=1, add=TRUE)

## ----meuse_uk_orig_plot, out.width='100%', fig.dim=c(10,10), fig.align='center'----
# prediction bias adjustment from log scale
g_uk_orig = exp(g_uk + g_uk_var/2)

# points on original scale
pts_orig = meuse_sf[['soils']]['zinc']

# prediction plot
zlim = range(exp(g), na.rm=TRUE)
plot(g_uk_orig, zlab='zinc (ppm)', main='[zinc] predictions and observations', cex=1.5, zlim=zlim)
plot(meuse_sf[['river_line']], add=TRUE)

# full plot
plot(g_uk_orig, zlab='zinc (ppm)', main='[zinc] predictions and observations', cex=1.5, zlim=zlim, reset=FALSE)
plot(meuse_sf[['river_line']], add=TRUE)

# overlay observation points
plot(pts_orig, add=TRUE, pch=16, pal=.pal)
plot(sf::st_geometry(pts_orig), add=TRUE)

