# sk.R
# Dean Koch, 2022
# S3 class for sk objects (grid lists)
#
#'
#' Make a snapKrig grid list object
#'
#' Constructs snapKrig ("sk") class list, representing a 2-dimensional regular spatial grid
#'
#' This function accepts 'RasterLayer' and 'RasterStack' inputs from the `raster` package,
#' 'SpatRaster' objects from `terra`, as well as any non-complex matrix, or a set of arguments
#' defining the vectorization of one. It returns a sk class list containing at least the
#' following three elements:
#'
#' * `gdim`: vector of two positive integers, the number of grid lines (n = their product)
#' * `gres`: vector of two positive scalars, the resolution (in distance between grid lines)
#' * `gyx`: list of two numeric vectors (lengths matching gdim), the grid line intercepts
#'
#' and optionally,
#'
#' * `crs`: character, the WKT representation of the CRS for the grid (optional)
#' * `gval`: numeric vector or matrix, the grid data
#' * `is_obs`: logical vector indicating non-NA entries in the grid data
#' * `idx_grid`: length-n numeric vector mapping rows of `gval` to grid points
#'
#' Supply some/all of these elements (including at least one of `gdim` or `gyx`) as named
#' arguments to `...`. The function will fill in missing entries wherever possible.
#'
#' If `gres` is missing, it is computed from the first two grid lines in `gyx`; If `gyx` is
#' missing, it is assigned the sequence `1:n` (scaled by `gres`, if available) for each `n`
#' in `gdim`; and if `gdim` is missing, it is set to the number of grid lines specified in
#' `gyx`. `gyx` should be sorted (ascending order), regularly spaced (with spacing `gres`),
#' and have lengths matching `gdim`.
#'
#' Scalar inputs to `gdim`, `gres` are duplicated for both dimensions. For example the call
#' `sk(gdim=c(5,5))` can be simplified to `sk(gdim=5)`, or `sk(5)`.
#'
#' For convenience, arguments can also be supplied together in a named list passed to `...`.
#' If a single unnamed argument is supplied (and it is not a list) the function expects it to
#' be either a numeric vector (`gdim`), a matrix, or a raster object.
#'
#' Alternatively, you can supply an `sk` object as (unnamed) first argument, followed by
#' individual named arguments. This replaces the named elements in the `sk` object then does
#' a validity check.
#'
#' Empty grids - with all data `NA` - can be initialized by setting `vals=FALSE`, in which case
#' `gval` will be absent from the returned list). Otherwise `gval` is the
#' column-vectorized grid data, either as a numeric vector (single layer case only) or as a
#' matrix with grid data in columns. `gval` is always accompanied by `is_obs`, which supplies
#' an index of `NA` entries (or rows)
#'
#' A sparse representation is used when `gval` is a matrix, where only the non-`NA` entries (or
#' rows) are stored. `idx_grid` in this case contains `NA`'s were `is_obs` is `FALSE`, and
#' otherwise contains the integer index of the corresponding row in `gval`. In the matrix case
#' it is assumed that each layer (ie column) has the same `NA` structure. `idx_grid` is only
#' computed for the first layer. If a point is missing from one layer, it should be missing
#' from all layers.
#'
#' @param ... raster, matrix, numeric vector, or list of named arguments (see details)
#' @param vals logical indicating to include the data vector `gval` in return list
#'
#' @return a "sk" class list object
#' @export
#' @family sk constructors
#'
#' @examples
#'
#' # simple grid construction from dimensions
#' gdim = c(12, 10)
#' g = sk(gdim)
#' summary(g)
#'
#' # pass result to sk and get the same thing back
#' identical(g, sk(g))
#'
#' # supply grid lines as named argument instead and get the same result
#' all.equal(g, sk(gyx=lapply(gdim, function(x) seq(x)-1L)))
#'
#' # display coordinates and grid line indices
#' plot(g)
#' plot(g, ij=TRUE)
#'
#' # same dimensions, different resolution, affecting aspect ratio in plot
#' gres_new = c(3, 4)
#' plot(sk(gdim=gdim, gres=gres_new))
#'
#' # single argument (unnamed) can be grid dimensions, with shorthand for square grids
#' all.equal(sk(gdim=c(2,2)), sk(c(2,2)))
#' all.equal(sk(2), sk(gdim=c(2,2)))
#'
#' # example with matrix data
#' gdim = c(25, 25)
#' gyx = as.list(expand.grid(lapply(gdim, seq)))
#' eg_vec = as.numeric( gyx[[2]] %% gyx[[1]] )
#' eg_mat = matrix(eg_vec, gdim)
#' g = sk(eg_mat)
#' plot(g, ij=TRUE, zlab='j mod i')
#'
#' # y is in descending order
#' plot(g, xlab='x = j', ylab='y = 26 - i', zlab='j mod i')
#'
#' # this is R's default matrix vectorization order
#' all.equal(eg_vec, as.vector(eg_mat))
#' all.equal(g, sk(gdim=gdim, gval=as.vector(eg_mat)))
#'
#' # multi-layer example from matrix
#' n_pt = prod(gdim)
#' n_layer = 3
#' mat_multi = matrix(stats::rnorm(n_pt*n_layer), n_pt, n_layer)
#' g_multi = sk(gdim=gdim, gval=mat_multi)
#' summary(g_multi)
#'
#' # repeat with missing data (note all columns must have consistent NA structure)
#' mat_multi[sample.int(n_pt, 0.5*n_pt),] = NA
#' g_multi_miss = sk(gdim=gdim, gval=mat_multi)
#' summary(g_multi_miss)
#'
#' # only observed data points are stored, idx_grid maps them to the full grid vector
#' max(abs( g_multi[['gval']] - g_multi_miss[['gval']][g_multi_miss[['idx_grid']],] ), na.rm=TRUE)
#'
#' # single bracket indexing magic does the mapping automatically
#' max(abs( g_multi[] - g_multi_miss[] ), na.rm=TRUE)
#'
#' # vals=FALSE drops multi-layer information
#' sk(gdim=gdim, gval=mat_multi, vals=FALSE)
#'
#' # raster import examples skipped to keep CMD check time < 5s on slower machines
#' \donttest{
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('external/rlogo.grd', package='raster')
#' r = raster::raster(r_path)
#'
#' # convert to sk (notice only first layer was loaded by raster)
#' g = sk(r)
#' summary(g)
#' plot(g)
#'
#' # open a RasterStack - gval becomes a matrix with layers in columns
#' r_multi = raster::stack(r_path)
#' g_multi = sk(r_multi)
#' summary(g_multi)
#' plot(g_multi, layer=1)
#' plot(g_multi, layer=2)
#' plot(g_multi, layer=3)
#'
#' # repeat with terra
#' if( requireNamespace('terra') ) {
#'
#' # open example file as SpatRaster (all layers loaded by default)
#' r_multi = terra::rast(r_path)
#' g_multi = sk(r_multi)
#' summary(g_multi)
#'
#' # open first layer only
#' g = sk(r[[1]])
#' summary(g)
#'
#' }
#' }
#' }
sk = function(..., vals=TRUE)
{
  # collapse the list when only one argument supplied to dots
  if( (...length() == 1 ) & all(is.null(...names()) ) ) { g = ..1 } else {

    # collect loose arguments into a list
    g = list(...)

    # case when an sk object is passed first, then some named arguments
    if( inherits(g[[1]], 'sk') & (length(g) > 1) ) {

      # drop attributes from sk object then replace any named arguments
      g = g[[1L]] |> c() |> utils::modifyList(g[-1L])
    }
  }

  # handle terra and raster objects
  is_terra = inherits(g, 'SpatRaster')
  is_raster_or_terra = inherits(g, c('RasterLayer', 'RasterStack')) | is_terra
  if( is_raster_or_terra )
  {
    # make sure the required package is loaded
    is_valid = if(is_terra) { requireNamespace('terra', quietly=TRUE) } else {
      requireNamespace('raster', quietly=TRUE)
    }
    if(!is_valid) stop(paste('package', ifelse(is_terra, 'terra', 'raster'), 'not loaded'))

    # order in dimensions is y, x, like in sk
    gdim = dim(g)[1:2]
    gval = NULL

    # terra and raster use another ordering for vectorization
    if(vals) idx_v = matrix(seq(prod(gdim)), gdim, byrow=TRUE)

    # look up the required package-specific functions
    res_ = if(is_terra) terra::res else raster::res
    yFromRow_ = if(is_terra) terra::yFromRow else raster::yFromRow
    xFromCol_ = if(is_terra) terra::xFromCol else raster::xFromCol
    crs_ = if(is_terra) terra::crs else raster::wkt
    nlyr_ = if(is_terra) terra::nlyr else raster::nlayers
    values_ = if(is_terra) terra::values else raster::getValues

    # copy CRS, number of layers, resolution, grid line locations
    if(vals) n_layer = nlyr_(g)
    gcrs = crs_(g)
    gres = res_(g)[2:1] # need this in order dy, dx
    gyx = lapply( list( y = yFromRow_(g, seq(gdim[1])), x = xFromCol_(g, seq(gdim[2])) ), sort)

    # copy vectorized data in correct order
    gval = NULL
    if(vals) gval = if(n_layer == 1) { values_(g)[idx_v] } else { values_(g)[idx_v,] }

    # initialize list
    g = list(gdim=gdim, gres=gres, gyx=gyx, crs=crs_(g), gval=gval)
  }

  # handle matrix objects
  if( is.matrix(g) ) g = list(gdim = dim(g), gval = if(vals) as.vector(g) else NULL )

  # handle numeric vectors
  if( is.numeric(g) )
  {
    # these are interpreted as grid dimensions
    if( length(g) > 2 ) stop('numeric vector g must be of length 1 or 2')
    g = list(gdim=g)
  }

  # halt on unrecognized objects
  if( !is.list(g) ) { stop('input g was not recognized') } else {

    # omit values if requested
    if(!vals) g[['gval']] = NULL

    # construct then validate the sk object
    return(sk_validate(sk_make(g)))
  }
}

#'
#'
#' Make a sk grid object
#'
#' This constructs a "sk" object from a named list containing at least the element `gdim`
#' or `gyx`. Users can optionally provide other list elements `gres`, `gval`, `crs`, `is_obs`,
#' and `idx_grid`.
#'
#' Input classes and lengths are checked before returning. `gdim` and `gres` are length-2
#' vectors (with y and x elements) but they can each be specified by a single number, as
#' shorthand to use the same value for y and x. `gval` should be a matrix or vector of grid
#' data, and `crs` should be a character string (the WKT representation).
#'
#' @param g list with any of the seven named arguments mentioned above (`gdim`, etc)
#'
#' @return a "sk" object
#' @export
#' @keywords internal
#'
#' @examples
#'
#' # auto-print reminds users to validate
#' sk_make(list(gdim=10, gres=0.5))
#'
sk_make = function(g)
{
  # compute gdim as needed
  if( is.null(g[['gdim']]) )
  {
    # gdim or gyx is required
    if( is.null(g[['gyx']]) ) stop('argument gdim not found')
    g[['gdim']] = sapply(g[['gyx']], length)
  }

  # check class and length of gdim and duplicate if it has length 1
  if( !is.numeric(g[['gdim']]) ) stop('gdim was not numeric')
  g[['gdim']] = sapply(g[['gdim']], as.integer)
  if( length(g[['gdim']]) == 1 ) g[['gdim']] = rep(g[['gdim']], 2)
  if( length(g[['gdim']]) !=2 ) stop('gdim must be length 2')

  # check class and length of grid resolution, duplicate if it has length 1
  if( !is.null(g[['gres']]) )
  {
    if( !is.numeric(g[['gres']]) ) stop('gres was not numeric')
    if( length(g[['gres']]) == 1 ) g[['gres']] = rep(g[['gres']], 2)
    if( length(g[['gres']]) !=2 ) stop('gres must be length 2')
  }

  # check class and length of grid line positions
  if( !is.null(g[['gyx']]) )
  {
    if( !is.list(g[['gyx']]) ) stop('g$gyx was not a list')
    if( !all( sapply(g[['gyx']], is.numeric) ) ) stop('non-numeric entries in gyx')
    if( length(g[['gyx']]) !=2 ) stop('gyx must be length 2')
  }

  # check class and length of CRS string
  if( !is.null(g[['crs']]) )
  {
    if( !is.character(g[['crs']]) ) stop('non-character crs')
    if( length(g[['crs']]) !=1 ) stop('crs must be a single character string')
  }

  # check for values and sparse representation index
  any_gval = !is.null(g[['gval']])
  is_sparse = !is.null(g[['idx_grid']])
  na_mapped = !is.null(g[['is_obs']])

  # check gval input when it is supplied
  if(!any_gval)
  {
    # drop unneeded indexing vectors
    g[['is_obs']] = NULL
    g[['idx_grid']] = NULL

  } else {

    # convert factor gval to character
    if( is.factor(g[['gval']]) ) g[['gval']] = as.character(g[['gval']])

    # check class of values in gval
    is_multi = is.matrix(g[['gval']])
    is_g_valid = is.vector(g[['gval']]) | is_multi
    if(!is_g_valid) stop('gval was not a vector or matrix')

    # validity checks and/or build idx_grid
    if(is_multi)
    {
      # number of grid points and layers
      nrow_multi = nrow(g[['gval']])
      n_layer = ncol(g[['gval']])

      # create sparse representation when it is expected (matrix gval) but not found
      if(is_sparse)
      {
        # validity checks for idx_grid
        is_obs = !is.na(g[['idx_grid']])
        if( length(g[['idx_grid']]) != prod(g[['gdim']]) ) stop('idx_grid has the wrong length')
        if( nrow(g[['gval']]) != sum(is_obs) ) stop('wrong number of non-NA indices in idx_grid')
        if( !all(g[['idx_grid']][is_obs] == seq(nrow_multi)) ) stop('mismatch in idx_grid and gval')

        # if valid idx_grid is supplied, the trimmed copy of gval should now have no NAs
        if( anyNA(g[['gval']]) ) stop('gval should have no NAs when supplied with idx_grid')

      } else {

        # identify observed data in first layer and build an indexing vector from it
        is_obs_first = !is.na(g[['gval']][,1L])
        g[['idx_grid']] = match(seq(nrow(g[['gval']])), which(is_obs_first))

        # check if there is any data
        if(length(g[['idx_grid']]) == 0)
        {
          # if not then we don't bother with sparse representation
          g[['idx_grid']] = NULL
          g[['gval']] = NULL
          is_sparse = FALSE

        } else {

          # omit NA rows and set indexing flag
          g[['gval']] = matrix(g[['gval']][is_obs_first,], ncol=n_layer)
          is_sparse = TRUE

          # if valid idx_grid is supplied, the trimmed copy of gval should now have no NAs
          if( anyNA(g[['gval']]) ) stop('inconsistent pattern of NAs among layers')
        }
      }
    }

    # check class of sparse matrix indexing vector
    if(is_sparse) g[['idx_grid']] = as.integer(g[['idx_grid']])

    # check class of na indicator vector
    if(na_mapped) g[['is_obs']] = as.logical(g[['is_obs']])

  }

  # set attributes for the class before returning
  structure(g, class='sk')
}

#'
#' Check compatibility of entries in a sk grid object, and fill in any missing ones
#'
#' This constructs the object and fills missing entries. It then does some sanity checks
#' and computes the index of NA points (in list entry `is_obs`).
#'
#' The function removes/introduces `idx_grid` depending on whether `gval` is a vector
#' (single-layer case) or a matrix (usually a multi-layer case). If `idx_grid` is missing
#' and `gval` is a matrix, it is assumed to contain all grid-points (including NAs)
#'
#' The function also assigns dimension names in the order 'y', 'x' (unless otherwise
#' specified) for `gdim`, `gres`, and `gyx`.
#'
#' `res_tol` is used to check if the resolution `gres` is consistent with the spacing
#' of the grid lines. Only the spacing of the first two lines is computed - if the
#' relative error along either dimensions is greater `res_tol`, the function throws
#' an error.
#'
#' @param g a "sk" object or a list accepted by `sk_make`
#' @param res_tol positive numeric, tolerance validating resolution (see details)
#'
#' @return a validated "sk" object
#' @export
#' @keywords internal

#' @examples
#'
#' sk_validate(list(gdim=10, gres=0.5))
#' sk_validate(list(gval=stats::rnorm(10^2), gdim=10, gres=0.5))
sk_validate = function(g, res_tol=1e-6)
{
  # order of list entries in the output
  nm_order = c('gdim', 'gres', 'gyx', 'crs', 'gval', 'idx_grid', 'is_obs')

  # check class of g
  if(!is.list(g)) stop('g must be a list')

  # construct sk object and copy its names
  g = sk_make(g)
  g_names = names(g)

  # names for dimensional components
  nm_dim = c('y', 'x')

  # check for grid point values and sparse representation index
  any_gval = !is.null(g[['gval']])
  na_mapped = !is.null(g[['is_obs']])
  is_indexed = !is.null(g[['idx_grid']])
  is_multi = is.matrix(g[['gval']])
  if( is_indexed & !is_multi ) stop('gval must be a matrix when idx_grid is supplied')
  if( is_multi & !is_indexed ) stop('gval is a matrix but idx_grid not found')

  # check names and put in the expected order
  if( is.null(names(g[['gdim']]) ) ) names(g[['gdim']]) = nm_dim
  if( !all( nm_dim %in% names(g[['gdim']]) ) ) stop('unknown names in gdim (expected "y", "x")')
  g[['gdim']] = g[['gdim']][nm_dim]

  # check for either grid dimension == 1
  if( !all(g[['gdim']] > 1) ) stop('grid must have two or more grid lines along each dimension')

  # when grid resolution is missing calculate it from gyx or else set up default
  if( is.null(g[['gres']]) )
  {
    # compute from gyx where available (or set unit default)
    g[['gres']] = c(1, 1)
    if( !is.null(g[['gyx']]) ) g[['gres']] = as.numeric(sapply(g[['gyx']], function(r) diff(r)[1]))
  }

  # check names and put in the expected order
  if( is.null(names(g[['gres']]) ) ) names(g[['gres']]) = nm_dim
  if( !all( nm_dim %in% names(g[['gres']]) ) ) stop('unknown names in gres (expected "y", "x")')
  g[['gres']] = g[['gres']][nm_dim]

  # set up grid line positions if they're missing
  if( is.null(g[['gyx']]) )
  {
    g[['gyx']] = Map(function(d, r) as.numeric(r*(seq(d)-1L)), d=g[['gdim']], r=g[['gres']])
  }

  # check names and put in the expected order
  if( is.null(names(g[['gyx']]) ) ) names(g[['gyx']]) = nm_dim
  if( !all( nm_dim %in% names(g[['gyx']]) ) ) stop('unknown names in gyx (expected "y", "x")')
  g[['gyx']] = g[['gyx']][nm_dim]

  # consistency checks for grid line locations
  gyx_mismatch = !all(sapply(g[['gyx']], length) == g[['gdim']])
  if( gyx_mismatch ) stop('number of grid lines (gyx) did not match gdim')

  # and for resolution
  gyx_error = abs(as.numeric(sapply(g[['gyx']], function(r) diff(r)[1])) - g[['gres']])
  if( any( (gyx_error / g[['gres']]) > res_tol) ) stop('resolution (gres) not consistent with gyx')

  # compute is_obs and check idx_grid as needed
  n = as.integer(prod(g[['gdim']]))
  if( !any_gval ) { g[['is_obs']] = logical(n) } else {

    # recycle scalar gval to get vector (or matrix) of expected size
    n_grid = length(g[['gval']])
    if(n_grid == 1) {

      g[['gval']] = rep(g[['gval']], n)
      if(is_indexed) g[['gval']] = matrix(g[['gval']], ncol=1)
      n_grid = n
    }

    # vector case without indexing vector
    if(!is_indexed) {

      g[['is_obs']] = !is.na(g[['gval']])

    } else {

      # sparse representation
      n_grid = length(g[['idx_grid']])
      n_obs = nrow(g[['gval']])

      # omit NAs specified in idx_grid if present
      if(n_obs == n_grid) {

        # trim to specified observed subset and recompute its size
        g[['gval']][ g[['idx_grid']][ !is.na(g[['idx_grid']]) ] ]
        n_obs = nrow(g[['gval']])
      }

      # check for wrong number of NAs in indexing grid
      g[['is_obs']] = !is.na(g[['idx_grid']])
      if( (n_grid - n_obs) != sum(!g[['is_obs']]) ) stop('gval and idx_grid are incompatible')
    }

    # check for wrong length in gval (or idx_grid) given gdim
    if(n_grid != n) stop('gdim and gval are incompatible')
  }

  # make the ordering consistent before returning ([] drops the class attribute)
  return(structure(g[ nm_order[nm_order %in% names(g)] ], class='sk'))
}
