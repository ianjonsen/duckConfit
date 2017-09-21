sfilter = function(d,
                   ts = list(gps = 1, gls = 12, ptt = 2),
                   ...) {
  
  ## dplyr-enabled wrapper function for calling fit_ssm from https://github.com/ianjonsen/ssmTMB
  ##
  ## d - input prefiltered data as a tibble of individual tracks grouped by id in the expected format
  ## ts - specify a list of time steps (in h) for gps, gls, & ptt datasets
  ## ... - additional arguments passed to fit_ssm
  
  
  ## add code back in here or in ssmTMB to fit model with common.tau in cases where 
  ##    gamma estimate > 0.9 & deployment duration is < 1 month (eg some KIPE deployments)
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

wrapLon <- function(lon, lmin = -180)
  (lon - lmin) %% 360 + lmin

unwrapLon <- function(lon, lmin = -180)
  cumsum(c(wrapLon(lon[1], lmin), wrapLon(diff(lon))))

redo_sfilter <-
  function(ssm_obj,
           data,
           s.inc = 0.1,
           n.inc = 1,
           tries = 5,
           common.tau = FALSE) {
    

  ## redo ssm filter for tracks that failed to converge
  ##  trying incrementally different span & nu values.
  ##  Up to tries re-filter attemps are made. Function then
  ##  searches for cases where hat(gamma) > 0.9, deployment < 30 d,
  ##  device_type == PTT & refits sfilter w common.tau = TRUE. This
  ##  may reduce over-smoothing of these short & typically sparsely
  ##  observed tracks
  ##
  ## ssm_obj - ssm filtered output from sfilter()
  ## pfd     - pre-filtered data (input to sfilter)
  
    cf <-
      which(sapply(ssm_obj$ssm, function(x)
        length(x) == 1 || x$opt$conv != 0))
    cat(paste("\n", length(cf), " individuals failed to converge"))
    i <- 0
    spn <- 0.5
    n <- 10
    while (length(cf) > 0) {
      i <- i + 1
      spn <- spn + s.inc
      if (spn > 0.7)
        spn <- 0.7
      n = n + n.inc
      if (i == tries + 1)
        break
      cat(paste("\nattempt", i, "\n"))
      redo <- slice(ssm_obj, cf)
      
      if(length(group_vars(data)) == 1) { 
        redo_ssm <- right_join(data, redo, by = "id") %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id) %>%
          do(ssm = sfilter(., span = spn, nu = n))
      }
      else if(length(group_vars(data)) == 2) {
        redo_ssm <- right_join(data, redo, by = c("id", "stage")) %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id, stage) %>%
          do(ssm = sfilter(., span = spn, nu = n))        
      }
      ssm_obj$ssm[cf] <- redo_ssm$ssm
      cf <-
        which(sapply(ssm_obj$ssm, function(x)
          length(x) == 1 || x$opt$conv != 0))
      cat(paste("\n", length(cf), " individuals still not converged"))
    }
    if (common.tau) {
    ## deal with cases where gamma > 0.9, deployment < 30 days, device_type = PTT
    ##  try fitting with common.tau = TRUE
    hg <-
      which(
        sapply(ssm_obj$ssm, function(x)
          x$par["gamma", 1] > 0.9 &&
            difftime(max(x$data$date), min(x$data$date), units = "days") < 30 &&
            x$data$device_type[1] == "PTT")
      )
    cat(paste(
      "\n try estimating a common tau for lon & lat - ",
      length(hg),
      " cases\n"
    ))
    redo <- slice(ssm_obj, hg)
    if(length(group_vars(data)) == 1) { 
      redo_ssm <- right_join(data, redo, by = "id") %>%
        select(-ssm) %>%
        ungroup() %>%
        group_by(id) %>%
        do(ssm = sfilter(., span = spn, nu = n))
    }
    else if( length(group_vars(data)) == 2) {
      redo_ssm <- right_join(data, redo, by = c("id", "stage")) %>%
        select(-ssm) %>%
        ungroup() %>%
        group_by(id, stage) %>%
        do(ssm = sfilter(., span = spn, nu = n))        
    }
    ssm_obj$ssm[hg] <- redo_ssm$ssm
  }
  
  ssm_obj
}
