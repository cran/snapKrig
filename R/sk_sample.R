# sk_sample.R
# Dean Koch, 2022
# sample variograms

#' Theoretical variogram function
#'
#' Computes the value of the variogram function w, defined by covariance model `pars`
#' at the component y and x lags supplied in `d`.
#'
#' By definition w is Var( Z(s1) - Z(s2) ), where s1 and s2 are a pair of spatial
#' locations, and Z is the spatial process value. If Z is second-order stationary
#' then w only depends on the relative displacement, s1 - s2 = (dx, dy). `snapKrig`
#' models the variogram as W = 2 ( `eps` + `psill` ( 1 - cy(dy) cx(dx) ) ).
#'
#' `sk_vario_fun` evaluates this function using the correlogram functions (cy and cx),
#' partial sill (`psill`) and nugget (`eps`) defined in `pars`, over the displacement
#' (dy and dx) supplied in `d`.
#'
#' NOTE: w is twice the semi-variogram, usually denoted by greek letter gamma. Variogram
#' w is therefore often written 2*gamma. This can (and does) lead to confusion in the
#' literature about whether to include a factor 2 in downstream calculations.
#' This function multiplies the semi-variogram function by 2, returning the variogram w
#' (ie 2*gamma), NOT the semi-variogram.
#'
#' If `d` is a list, its 'y' and 'x' components should supply the y and x component distances.
#' These must be equal-length non-negative numeric vectors. The function returns the corresponding
#' variogram values in a vector of the same length.
#'
#' If `d` is a numeric vector, it is interpreted as a set of distances at which to
#' evaluate the range of the variogram function. Anisotropic variograms will exhibit a range
#' of values for a given distance (depending on the relative sizes of the x and y components).
#' The function returns this range in a data frame with columns 'min' and 'max'.
#'
#' @param pars list of the form returned by `sk_pars` with entries 'y', 'x', 'eps', 'psill'
#' @param d numeric vector or list with vector entries 'y' and 'x', the distances to evaluate
#'
#' @return data frame (for list `d`) or numeric vector (for vector `d`) of variogram values
#'
#' @export
#' @keywords internal
#' @family variogram functions
#' @seealso sk_pars
#'
#' @examples
#' # set up example grid and parameters
#' gdim = c(10, 15)
#' d_max = sqrt(sum(gdim^2))
#' pars = sk_pars(gdim, 'mat')
#'
#' # set up test distances
#' d_test = seq(0, d_max, length.out=1e2)
#'
#' # evaluate and plot the variogram values for equal displacements along x and y
#' d_equal = stats::setNames(rep(list(sqrt(1/2)*d_test), 2), c('y', 'x'))
#' vario = sk_vario_fun(pars, d=d_equal)
#' plot(d_test, vario, pch=NA)
#' lines(d_test, vario, col='blue')
#'
#' # evaluate and plot the range of variogram values (for all possible x and y displacements)
#' vario_lims = sk_vario_fun(pars, d=d_test)
#' lines(d_test, vario_lims[,1])
#' lines(d_test, vario_lims[,2])
#'
sk_vario_fun = function(pars, d=NULL)
{
  # 1d vector case
  if( !is.list(d) )
  {
    # component distances to test along each axis
    d0 = rep(0, length(d))

    # compute three sets of component distances: varying x, y, and both equally
    cov_y = sk_vario_fun( pars, list(y=d, x=d0) )
    cov_x = sk_vario_fun( pars, list(y=d0, x=d) )
    cov_yx = sk_vario_fun( pars, list(y=d/sqrt(2), x=d/sqrt(2)))

    # compute range of theoretical semi-variogram at each test distance
    cov_min = pmin(cov_y, cov_x, cov_yx)
    cov_max = pmax(cov_y, cov_x, cov_yx)
    return( data.frame(min=cov_min, max=cov_max) )
  }

  # take product of correlation functions
  nm_yx = c('y', 'x')
  if( !all(nm_yx %in% names(d)) ) stop('list elements y and x (in d) must be named')
  corrvals = do.call('*', Map(function(p, dyx) sk_corr(p, dyx), p=pars[nm_yx], dyx=d[nm_yx]))

  # return the variogram (twice the semi-variogram)
  return( 2 * pars[['eps']] + pars[['psill']] * ( 1 - corrvals ) )
}


#' Sub-grid point sampler for grid data
#'
#' Sample `n` locations from the non-NA points in the input grid `g`, optionally using
#' them as centers to place `n` sub-grids of the specified size and resolution.
#'
#' When `sk_out=TRUE` (the default), the function returns an sk grid containing the sampled
#' points. If multiple samples are requested, a multi-layer grid is returned. When
#' `sk_out=FALSE`, the function returns the vector index of the sampled grid points,
#' or if multiple samples are requested, a list of vectors.
#'
#' By default the function simply draws a sample of `n` locations (uniformly at random)
#' from the non-NA points in the input grid `g`.
#'
#' When `lag_max > 1`, the function instead returns the the Moore neighbourhood of
#' radius `lag_max` around each of the sample points (including the center point). These
#' sub-grids are returned as distinct layers (or list entries, if `sk_out=FALSE`). Their
#' resolution can be coarsened (up-scaled) by increasing `up` from its default 0. `up`
#' must either be 0 or else a positive integer that evenly divides `lag_max`
#' (see `sk_rescale`).
#'
#'
#' For a given `up`, the grid `g` can be partitioned into `(up+1)^2` distinct
#' non-overlapping sub-grids. When `over=FALSE` (the default), the function apportions
#' its `n` point samples as evenly as possible among these disjoint subsets. This ensures
#' that if `n` is less than or equal to `(up+1)^2`, and there are no `NA`s, there can be
#' no repetition (overlap) of points in the returned sub-grids.
#'
#' Note that with the default `sk_out=TRUE`, `lag_max > 1` is only supported for complete
#' grids `g`. This is because with missing data it is hard (and sometimes impossible) to
#' ensure that the Moore neighborhoods have identical `NA` structure (and this is a
#' requirement for multi-layer sk grids).
#'
#' Note also that multi-layer sk grids are not fully supported yet. If you pass a
#' multi-layer grid to g, the function returns results for the first layer only.
#'
#' @param g an sk grid object or any other object accepted by `sk`
#' @param n integer > 0, the maximum number of center points to sample
#' @param lag_max integer, Moore neighborhood radius (ie the maximum queen's distance)
#' @param up integer, the up-scaling factor for sampling sub-grids of `g`
#' @param over logical, indicates to allow overlapping sub-grids (when they can be avoided)
#' @param sk_out logical, if TRUE (the default) the function returns an sk grid
#' @param seed integer seed, passed to `base::set.seed`
#'
#' @return If `lag_max == 0` (the default), the function returns a single-layer sk grid when
#' `sk_out=TRUE`, or else the sample indices in `g` as a length-`n` integer vector. If
#' `lag_max > 0`, the function returns a multi-layer sk grid `sk_out=TRUE`, or else a list of
#' `n` vectors indexing the sampled points in each sub-grid of `g`.
#'
#' @export
#' @keywords internal
#' @seealso sk sk_sample_vg
#'
#' @examples
#' # define an empty grid
#' g_empty = sk(gdim = c(100, 100))
#'
#' # get an ordinary random sample with default settings
#' g_sample = sk_sample_pt(g_empty)
#' plot(g_sample)
#'
#' # same call with index return mode
#' idx_sample = sk_sample_pt(g_empty, sk_out=FALSE)
#' str(idx_sample)
#'
#' # reduce or increase number of center points from default 100
#' g_sample = sk_sample_pt(g_empty, n=10)
#' plot(g_sample)
#'
#' # add some data to g and repeat
#' pars = sk_pars(g_empty)
#' pars$eps = 1e-6
#' g = sk_sim(g_empty, pars)
#' plot(g)
#' g_sample = sk_sample_pt(g)
#' plot(g_sample)
#'
#' # sample 3 subgrids from Moore neighbourhoods of radius 6 (index output mode)
#' n = 3
#' idx_sample = sk_sample_pt(g, n=n, lag_max=6L, sk_out=FALSE, seed=42)
#'
#' # plot each list element a different color
#' group_sample = rep(0L, length(g))
#' for(i in seq(n)) group_sample[ idx_sample[[i]] ] = i
#' sk_plot(group_sample, dim(g), breaks=c('not sampled', seq(n)), zlab='sub-grid')
#'
#' # plot all the sub-grid data
#' g_plot = g_empty
#' g_plot[unlist(idx_sample)] = g[unlist(idx_sample)]
#' plot(g_plot)
#'
#' # default sk_out=TRUE returns them as multi-layer grid object
#' g_sample = sk_sample_pt(g, n=n, lag_max=6L, seed=42)
#' plot(g_sample, layer=1, zlim=range(g_plot, na.rm=TRUE))
#' plot(g_sample, layer=2, zlim=range(g_plot, na.rm=TRUE))
#' plot(g_sample, layer=3, zlim=range(g_plot, na.rm=TRUE))
#'
#'
#'
#' # When up > 0 the function will attempts to avoid overlap whenever possible
#' up = 1
#' n = (up+1)^2 # to get disjoint results n must be less than or equal to (up+1)^2
#' lag_max = 10 * (up+1) # vary to get larger/smaller subsets. max allowable: min(gdim)/2
#' idx_sample = sk_sample_pt(g, n=n, up=up, lag_max=lag_max, sk_out=FALSE)
#' idx_overlap = rowSums( sapply(idx_sample, function(i) seq_along(g) %in% i) )
#'
#' # plot each list element a different color
#' group_sample = rep(0L, length(g))
#' for(i in seq(n)) group_sample[ idx_sample[[i]] ] = i
#' sk_plot(group_sample, dim(g), breaks=c('not sampled', seq(n)), zlab='sub-grid')
#'
#' # no overlap
#' sk_plot(as.integer(idx_overlap), dim(g), zlab='times sampled')
#'
#' # compare with over=TRUE (usually results in overlap - try running a few times)
#' idx_sample_compare = sk_sample_pt(g, n=n, up=up, lag_max=lag_max, over=TRUE, sk_out=FALSE)
#' idx_overlap_compare = rowSums( sapply(idx_sample_compare, function(i) seq_along(g) %in% i) )
#' sk_plot(as.integer(idx_overlap_compare), dim(g), zlab='times sampled')
#'
#' # incomplete input data example
#' g_sample = sk_sample_pt(g, n=10)
#' sk_plot(g_sample)
#'
#' # draw a sample of center points and indicate sub-grids in color
#' idx_sample = sk_sample_pt(g_sample, n=10, lag_max=6, up=1, over=FALSE, sk_out=FALSE)
#' g_sample_grid = g_empty
#' g_sample_grid[] = rep('not sampled', length(g_empty))
#' g_sample_grid[unlist(idx_sample)] = 'sub-grid sample'
#' plot(g_sample_grid)
#'
sk_sample_pt = function(g, n=1e2, lag_max=0, up=0, over=FALSE, sk_out=TRUE, seed=NULL)
{
  # initialize RNG state
  set.seed(seed)

  # unpack the grid object
  g = sk(g)
  gdim = dim(g)
  is_obs = !is.na(g)

  # extract only the first layer from multi-layer objects
  if( is.matrix(g[['gval']]) ) g = sk(utils::modifyList(g, list(gval=g[,1L])))

  # all-missing case treated as all-observed
  all_missing = !any(is_obs)
  if( all_missing )
  {
    is_obs = rep(TRUE, length(is_obs))
    g[] = rep(TRUE, length(g))
  }

  # check for valid parameters in request
  gdim_inner = gdim - 2 * lag_max
  if( !( all(gdim_inner > 0 ) ) ) stop('lag_max cannot exceed min(dim(g))/2')

  # find subset of points for which Moore neighbourhood is completely inside grid
  ij_inner = lapply(stats::setNames(gdim, c('i', 'j')), function(d) seq(1+lag_max, d-lag_max) )
  is_inner = sk_sub_idx(gdim, ij=ij_inner)

  # compute their vector index and grid positions
  is_eligible = is_obs & is_inner
  idx_eligible = which(is_eligible)
  ij_obs = sk_vec2mat(idx_eligible, gdim)

  # sample center points only
  n = min(sum(is_eligible), n)
  idx_center_all = sort(sample(idx_eligible, n))

  # handle single layer output
  if(lag_max == 0)
  {
    # return from index mode
    if(!sk_out) return(idx_center_all)

    # sk grid output
    is_selected = seq_along(g) %in% idx_center_all
    g[!is_selected] = ifelse(all_missing, FALSE, NA)
    return(g)
  }

  # check for valid `up` and define sequence interval for sub-grid lines
  sep = 1L + as.integer(up)
  if(sep < 1) stop('up must be a non-negative integer')

  # pick a sample that avoids overlap
  if(!over)
  {
    # the number of disjoint sub-grids and the index of points
    n_disjoint = sep^2
    idx_inner = which(is_inner)

    # identify all disjoint sub-grid origins within inner sub-grid,
    ij_origin = apply(expand.grid(seq(sep)-1L, seq(sep)-1L), 1L, identity, simplify=FALSE)
    ij_all = lapply(ij_origin, function(ij) Map(function(d, r) seq(1L+r, d, by=sep), d=gdim_inner, r=ij))
    idx_disjoint = lapply(ij_all, function(ij) sk_sub_idx(gdim_inner, unname(ij), idx=TRUE))

    # omit NA points
    idx_disjoint = lapply(idx_disjoint, function(ij) ij[ is_obs[idx_inner][ij] ] )

    # the number of non-NA points in each disjoint sub-grid
    n_obs_disjoint = sapply(idx_disjoint, length)

    # apportion samples evenly to remaining sample sites, building n_per in a loop
    n_per = integer(n_disjoint)
    n_remain = n
    while( (n_remain > 0) & any(n_obs_disjoint > 0) )
    {
      # find minimum number of non-NA points remaining in nonempty sub-grids
      is_nonempty = n_obs_disjoint > 0
      n_nonempty = sum(is_nonempty)

      # draw the same number of points from each one
      n_each = min(min(n_obs_disjoint[is_nonempty]), floor(n_remain/n_nonempty))

      # exit case: fewer points needed than available sub-grids
      if(n_each == 0)
      {
        # draw a single point from a subset of the sub-grids
        n_each = 1
        n_rem = sample(which(is_nonempty), n_remain)
        is_nonempty = seq(n_disjoint) %in% n_rem
      }

      # draw at most this many points from each non-empty sub-grid
      n_per[is_nonempty] = n_per[is_nonempty] + n_each
      n_obs_disjoint = n_obs_disjoint - n_each
      n_remain = n_remain - n_each * sum(is_nonempty)
    }

    # sample n_per[i] sub-grid origins from ith disjoint set
    idx_center_inner = unlist( Map(function(n, idx) sample(idx, n), n=n_per, idx=idx_disjoint) )

    # remap to original (not inner) grid
    idx_center_all = idx_inner[idx_center_inner]
  }

  # check for invalid up-scaling factor then make center box template
  if( ( lag_max %% sep ) != 0) stop('up+1 must divide lag_max')
  lag_max = min(c(lag_max, gdim))
  box_offset = seq(-lag_max, lag_max, by=sep)

  # make a list of box grid lines
  ij_list = lapply(idx_center_all, function(idx) {

    # find grid (i,j) indices for box template centered at center point
    ij_center = lapply(sk_vec2mat(idx, gdim, out='list'), function(idx) idx + box_offset)
  })

  # return from sk grid mode
  if(sk_out)
  {
    # create single-layer sk (sub)grids and extract their data into columns of matrix
    g_sub_mat = sapply(ij_list, function(ij) sk_sub(g, ij)[])

    # merge them into a multi-layer object
    g_sub = sk_sub(g, ij_list[[1]])
    return( sk(utils::modifyList(g_sub, list(gval=g_sub_mat))) )
  }

  # loop over boxes, computing index of all points in Moore neighbourhood
  idx_list = lapply(ij_list, function(ij) sk_sub_idx(gdim, ij, idx=TRUE))
  return(idx_list)
}


#' Sample point pair absolute differences for use in semi-variogram estimation
#'
#' Compute the absolute differences for point pairs in `g`, along with their separation
#' distances. If no sample point index is supplied (in `idx`), the function samples points
#' at random using `sk_sample_pt`.
#'
#' In a set of n points there are n_pp(n) = (n^2-n)/2 possible point pairs. This
#' expression is inverted to determine the maximum number of sample points in `g` to use
#' in order to satisfy the argument `n_pp`, the maximum number of point pairs to sample.
#' A random sub-sample of `idx` is taken as needed. By default `n_pp=1e4` which results
#' in `n=141`.
#'
#' The mean of the point pair absolute values ('dabs') for a given distance interval is the
#' classical estimator of the variogram. This and two other robust methods are implemented
#' in `sk_plot_semi`. These statistics are sensitive to the choice of distance bins. They
#' are added automatically by a call to `sk_add_bins` (with `n_bin`) but users can also set
#' up bins manually by adjusting the 'bin' column of the output.
#'
#' For multi-layer `g`, the function samples observed point locations once and re-uses this
#' selection in all layers. At most `n_layer_max` layers are sampled in this way (default is
#' the square root of the number of layers, rounded up)
#'
#' @param g any grid object accepted or returned by `sk`
#' @param n_pp integer maximum number of point pairs to sample
#' @param idx optional integer vector indexing the points to sample
#' @param n_bin integer number of distance bins to assign (passed to `sk_add_bins`)
#' @param n_layer_max integer, maximum number of layers to sample (for multi-layer `g`)
#' @param quiet logical, suppresses console output
#'
#'
#' @return A data frame with a row for each sampled point pair. Fields include 'dabs' and 'd',
#' the absolute difference in point values and the separation distance, along with the vector
#' index, row and column numbers, and component (x, y) distances for each point pair. 'bin'
#' indicates membership in one of `n_bin` categories.
#'
#' @export
#' @seealso sk sk_sample_pt sk_add_bins
#'
#' @examples
#'
#' # make example grid and reference covariance model
#' gdim = c(22, 15)
#' n = prod(gdim)
#' g_empty = sk(gdim)
#' pars = sk_pars(g_empty, 'mat')
#'
#' # generate sample data and sample semi-variogram
#' g_obs = sk_sim(g_empty, pars)
#' vg = sk_sample_vg(g_obs)
#' str(vg)
#'
#' # pass to plotter and overlay the model that generated the data
#' sk_plot_semi(vg, pars)
#'
#' # repeat with smaller sample sizes
#' sk_plot_semi(sk_sample_vg(g_obs, 1e2), pars)
#' sk_plot_semi(sk_sample_vg(g_obs, 1e3), pars)
#'
#' # use a set of specific points
#' n_sp = 10
#' ( n_sp^2 - n_sp ) / 2 # the number of point pairs
#' vg = sk_sample_vg(g_obs, idx=sample.int(n, n_sp))
#' sk_plot_semi(vg, pars)
#'
#' # non-essential examples skipped to stay below 5s exec time on slower machines
#' \donttest{
#'
#' # repeat with all point pairs sampled (not recommended for big data sets)
#' vg = sk_sample_vg(g_obs, n_pp=Inf)
#' sk_plot_semi(vg, pars)
#' ( n^2 - n ) / 2 # the number of point pairs
#'
#' ## example with multiple layers
#'
#' # generate five layers
#' g_obs_multi = sk_sim(g_empty, pars, n_layer=5)
#'
#' # by default, a sub-sample of sqrt(n_layers) is selected
#' vg = sk_sample_vg(g_obs_multi)
#' sk_plot_semi(vg, pars)
#'
#' # change this behaviour with n_layer_max
#' vg = sk_sample_vg(g_obs_multi, n_layer_max=5)
#' sk_plot_semi(vg, pars)
#'
#' }
#'
sk_sample_vg = function(g, n_pp=1e4, idx=NULL, n_bin=25, n_layer_max=NA, quiet=FALSE)
{
  # unpack expected inputs
  g = sk(g)
  is_obs = !is.na(g)

  # count observations and halt if there aren't any
  n = sum(is_obs)
  if(n == 0) stop('g had no non-NA values')

  # the number of sample points required to get n point pairs, the inverse of f(n)=(n^2-n)/2
  n_obs = min(floor( (1 + sqrt(1 + 8*n_pp)) / 2 ), n)

  # call point sampler if a subset of points was not specified
  if( is.null(idx) ) idx = sk_sample_pt(g, n_obs, sk_out=FALSE)

  # verify that input sample point locations are non-NA and sub-sample as needed
  idx = idx[ is_obs[idx] ]
  if( length(idx) > n_obs ) idx = sample(idx, n_obs)
  n_obs = length(idx)
  n_pp = (n_obs^2-n_obs)/2

  # console output
  msg_n = paste0('sampling ', n_obs, ' of ', n, ' non-NA points (', n_pp, ' point pairs)\n')
  if(!quiet) cat(msg_n)

  # compute grid indices (i, j) for the sample points
  ij_obs = sk_vec2mat(idx, dim(g))

  # anonymous function to compute lower triangular part of a 1d distance matrix
  vec_dist = function(x) c(stats::dist(x, method='manhattan', diag=TRUE))

  # compute dimension-wise grid line (i,j) distances between all point pairs
  dmat = apply(ij_obs, 2, vec_dist, simplify=FALSE)
  names(dmat) = c('di', 'dj')

  # indexes the distance matrix entries vectorized by the `c` in vec_dist
  idx_lower = lower.tri(matrix(NA, n_obs, n_obs))

  # i,j indices corresponding to vectorized entries of idx_lower
  idx_i_lower = matrix(rep(seq(n_obs), n_obs), n_obs)[idx_lower]
  idx_j_lower = matrix(rep(seq(n_obs), each=n_obs), n_obs)[idx_lower]

  # i,j indices for each point pair, and their vectorized indices
  ij_pp = list(p1=ij_obs[idx_i_lower,], p2=ij_obs[idx_j_lower,])
  idx_pp = lapply(ij_pp, function(ij) sk_mat2vec(ij, dim(g)))

  # compute absolute differences in data for each point pair, then separation distance
  is_multi = is.matrix(g[['gval']])
  if( is_multi )
  {
    # speed things up by directly indexing matrix of non-NAs
    idx_pp_sparse = lapply(idx_pp, function(i) g[['idx_grid']][i])

    # sub-sample among layers
    n_layer = ncol(g[['gval']])
    if( is.na(n_layer_max) ) n_layer_max = max(1L, ceiling(sqrt(n_layer)))
    n_layer_max = min(n_layer, n_layer_max)
    idx_sample_layer = sample.int(n_layer, n_layer_max)

    # console output
    if(!quiet) cat(paste('sampling', n_layer_max, 'of', n_layer, 'layers\n'))

    # loop over selected sample layers, drawing the same point pairs in each layer
    dabs_all = sapply(idx_sample_layer, function(idx_layer) {

      # data values for selected point pairs in a big matrix
      z_mat = sapply(idx_pp_sparse, function(i) as.vector(g[i, idx_layer]))

      # absolute differences for point pairs in this layer
      abs(apply(z_mat, 1, diff))
    })

    lyr = rep(idx_sample_layer, each=nrow(dabs_all))

  } else {

    # one call to do all of the above for a single layer only
    dabs_all = abs( apply( sapply(idx_pp, function(i) g[i]), 1, diff ) )
    lyr = 1L
  }

  # compile everything in a data frame then append distance and layer number
  vg = data.frame(dabs=c(dabs_all), idx_pp, ij_pp, dmat)
  vg[c('dy', 'dx')] = rep(g[['gres']], each=nrow(vg)) * vg[c('di', 'dj')]
  vg[['lyr']] = lyr
  vg[['d']] = sqrt( rowSums( vg[c('dy', 'dx')]^2 ) )

  # split into distance bins before returning
  return( sk_add_bins(vg, n_bin) )
}


#' Add bin labels to a variogram data frame
#'
#' Helper function for grouping the rows of input data frame `vg` into `n_bins` bins
#' according to the value of the (numeric) distance column `d`. This uses either `base::cut`
#' or, if `probs` is supplied, `stats::quantile`.
#'
#' By default, the function sets `probs` to a sequence of length `1+n_bin` evenly splitting
#' the interval \[0,1\] to ensure approximately equal sample sizes for each bin. Setting
#' `probs=NA` instead sets the bin endpoints such that the range of distances is split
#' evenly (note this may produce empty bins)
#'
#' The function is called by `sk_sample_vg` and `sk_plot_semi` (when column `bin` is
#' missing). It can also be used to recompute bins after an `rbind` of multiple variogram
#' data frames.
#'
#' @param vg data frame with numeric column 'd'
#' @param n_bin integer number of distance bins to assign
#' @param probs numeric vector of quantile probabilities to establish breakpoints (length `n_bin+1`)
#'
#' @return same as input `vg` but with integer column `bin` added/modified
#'
#' @export
#' @keywords internal
#' @family variogram functions
#'
#' @examples
#' distance_df = data.frame(d=runif(25))
#' sk_add_bins(distance_df)
#'
#' # specify fewer bins and set up quantiles explicitly
#' sk_add_bins(distance_df, n_bin = 5) # same as ...
#' sk_add_bins(distance_df, n_bin = 5, probs=seq(0, 1, length.out=6))
#'
#' # break range of distances into evenly spaced bins (of varying sample sizes)
#' sk_add_bins(distance_df, n_bin = 5, probs=NULL)
#'
sk_add_bins = function(vg, n_bin=25, probs=NULL)
{
  # check for distance column
  d = vg[['d']]
  if( is.null(d) ) stop('missing distance column "d" in variogram data frame')

  # check for valid input
  if(n_bin > nrow(vg))
  {
    n_bin = nrow(vg)
    warning(paste('n_bin reduced to', n_bin))
  }

  # assign default quantile breaks
  if( is.null(NULL) ) probs = seq(0, 1, length.out=1L+n_bin)

  # n_bin=1 puts everything in bin 1
  if(n_bin == 1) return( utils::modifyList(vg, list(bin=1L)) )

  # split range evenly into bins when `probs=NULL`
  if( anyNA(probs) ) return( utils::modifyList(vg, list( bin=as.integer(cut(d, breaks=n_bin))) ) )

  # otherwise find bins of roughly equal content
  return( utils::modifyList(vg, list(bin=findInterval(d, stats::quantile(d, probs=probs)))) )
}
