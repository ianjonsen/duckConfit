##' State-space filter
##'
##' Wrapper function that takes data prepared by \code{prefilter}, assigns appropriate
##' time steps for GPS, GLS and/or PTT devices & calls \code{ssmTMB::fit_ssm} to do the
##' filtering. \code{ssmTMB} can be installed via \code{devtools::install_github("ianjonsen/ssmTMB")}.
##'
##' @title sfilter
##' @param d input data from \code{prefilter} as a tibble of individual tracks grouped by id
##' @param ts specify a list of time steps (in h) for gps, gls, & ptt datasets
##' @param ... additional arguments passed to \code{ssmTMB::fit_ssm}
##' @return a list with 10 elements (see ??fit_ssm)
##'
##' @examples
##' \dontrun{
##' pfd <- prefilter(
##'   sp = "SOES",
##'   min_obs = 30,
##'   min_days = 5,
##'   vmax = 10,
##'   path2data = "~/Dropbox/r"
##'   )
##'
##' ssm_by_id <- pfd %>%
##'   dplyr::do(ssm = sfilter(., span = 0.4, nu = 5))
##'
##' ## try re-filtering tracks that failed to converge
##' ssm_by_id <- redo_sfilter(ssm_by_id, pfd, tries = 10)
##' }
##'
##' @importFrom dplyr mutate
##' @importFrom ssmTMB fit_ssm
##'

sfilter <- function(d,
                   ts = list(gps = 1, gls = 12, ptt = 2),
                   ...) {

  switch(unique(d$device_type),
         GPS = {
             ts = ts$gps
         },
         GLS = {
             ts = ts$gls
         },
         PTT = {
             ts = ts$ptt
         })

  ## when track straddles -180,+180 shift to 0,360
  f <- subset(d, keep)

  if(diff(range(f$lon)) > 300) {
    tmp <- f %>% mutate(lon = wrapLon(lon, lmin = 0))
    if(min(diff(tmp$lon)) < -310 || max(diff(d$lon)) > 310) {
      d[d$keep, ] <- f %>% mutate(lon = unwrapLon(lon, lmin = 0))
    }
    else{
      d[d$keep, ] <- tmp
    }
  }

  out <- try(fit_ssm(d, subset = d$keep, tstep = ts / 24, ...), silent = TRUE)

  ## if track on 0,360+ then shift back to -180,+180
  if(length(out) > 1) {
    out$predicted <- out$predicted %>% mutate(lon = wrapLon(lon, lmin = -180))
    out$fitted <- out$fitted %>% mutate(lon = wrapLon(lon, lmin = -180))
    out$data <- out$data %>% mutate(lon = wrapLon(lon, lmin = -180))
  }

  out
}

