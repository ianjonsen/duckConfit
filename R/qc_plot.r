##' Generate quality control plots
##'
##' Prepare standardised track data for state-space filtering by performing
##' further standardisation and culling unacceptably short deployments
##'
##' The track data are fetched automatically from the RAATD repository by specifying
##' a species' abbreviated name.
##'
##' @title qc_plot
##' @param ssm.tbl a compound tibble of SSM filter output grouped by \code{id} or
##'   \code{id} and \code{stage}
##' @param sp the species' 4-letter abbreviated name
##' @param filt_rnd the filtering round for which plots are generated (one of: 1, 1b, 2, 3)
##' @return a .pdf file of quality control plots
##'
##' @examples
##' \dontrun{
##' ssm_by_id %>% qc_plot(sp = "soes")
##' }
##' @importFrom ggplot2 aes map_data ggplot geom_map theme ylim xlim geom_path geom_point ggtitle geom_line geom_rug
##' @importFrom dplyr mutate do
##' @importFrom gridExtra grid.arrange

qc_plot <- function(ssm.tbl,
                    sp,
                    filt_rnd = 1) {
  ## generate a .pdf of SSM fits to data for qc

  require(mapdata, quietly = TRUE)

  map.world <- map_data(map = "world")


  plt.fn <- function(d) {

    ssm <- d$ssm[[1]]
    p <- ggplot() + geom_map(
      data = map.world,
      map = map.world,
      aes(map_id = "Antarctica"),
      fill = "snow"
    )

    if (length(ssm) > 1) {
      loc.dat <- subset(ssm$data, keep)
      loc.pred <- ssm$predicted
      wl <- FALSE

      x.rng <- extendrange(loc.pred$lon, f = 0.1)
      if (diff(x.rng) > 300)
        wl <- TRUE

      if (wl) {
        tlon.d <- loc.dat %>% mutate(lon = wrapLon(lon, lmin = 0))
        tlon.p <- loc.pred %>% mutate(lon = wrapLon(lon, lmin = 0))
        if (min(diff(tlon.d$lon)) < -180 ||
            max(diff(tlon.d$lon > 180))) {
          loc.dat <- loc.dat %>% mutate(lon = unwrapLon(lon, lmin = 0))
          loc.pred <-
            loc.pred %>% mutate(lon = unwrapLon(lon, lmin = 0))
        }
        else{
          loc.dat <- tlon.d
          loc.pred <- tlon.p
        }
        x.rng <- extendrange(loc.pred$lon, f = 0.1)
      }

      y.rng <- extendrange(loc.pred$lat, f = 0.05)

      p1 <- p +
        theme(legend.position = "none") + ylim(y.rng) + xlim(x.rng) +
        geom_point(
          aes(x = lon, y = lat),
          data = loc.dat,
          col = "firebrick",
          size = 1
        )

      p1 <-
        p1 + geom_path(data = loc.pred,
                       aes(x = lon, y = lat),
                       col = "dodgerblue",
                       lwd = 0.25) +
        geom_point(
          data = loc.pred,
          aes(x = lon, y = lat),
          col = "dodgerblue",
          size = 0.5
        ) +
        ggtitle(
          paste(
            "ID = ",
            d$id,
            "\n stage = ",
            d$stage,
            "\n device_type = ",
            unique(ssm$data$device_type),
            "; time step = ",
            ssm$tstep * 24,
            " h\n",
            "optimisation result: ",
            ssm$opt$convergence,
            "; ",
            ssm$opt$message,
            "\n gamma: ",
            round(ssm$par["gamma", "Estimate"], 3)
            ,
            sep = ""
          )
        )

      p2 <-
        ggplot(data = loc.dat, aes(y = lon, x = date)) +
        ylim(x.rng) +
        geom_point(col = "firebrick", size = 1)

      p2 <- p2 + geom_rug(col = "firebrick", sides = "b") +
        geom_line(
          data = loc.pred,
          aes(x = date, y = lon + lon.se * 1.96),
          lwd = 0.25,
          col = "dodgerblue"
        ) +
        geom_line(
          data = loc.pred,
          aes(x = date, y = lon - lon.se * 1.96),
          lwd = 0.25,
          col = "dodgerblue"
        ) +
        geom_point(
          data = loc.pred,
          aes(x = date, y = lon),
          col = "dodgerblue",
          size = 0.1
        ) +
        geom_rug(data = loc.pred,
                 col = "dodgerblue",
                 sides = "t")

      p3 <-
        ggplot(loc.dat, aes(y = lat, x = date)) + ylim(y.rng) +
        geom_point(col = "firebrick", size = 1)

      p3 <- p3 + geom_rug(col = "firebrick", sides = "b") +
        geom_line(
          data = loc.pred,
          aes(x = date, y = lat + lat.se * 1.96),
          lwd = 0.25,
          col = "dodgerblue"
        ) +
        geom_line(
          data = loc.pred,
          aes(x = date, y = lat - lat.se * 1.96),
          lwd = 0.25,
          col = "dodgerblue"
        ) +
        geom_point(
          data = loc.pred,
          aes(x = date, y = lat),
          col = "dodgerblue",
          size = 0.1
        ) +
        geom_rug(data = loc.pred,
                 col = "dodgerblue",
                 sides = "t")

      grid.arrange(p1, p2, p3, heights = c(2, 1, 1))
    }
    else if (length(ssm) == 1) {

      plot(
        0,
        0,
        type = 'n',
        main = paste(d$id, "\n", d$ssm),
        axes = F,
        xlab = "",
        ylab = ""
      )
    }
  }

  cat("\n generating pdf plots...\n")

  if(filt_rnd == 1) {
    pdf(
      file = paste(
        "../raatd_data/data_filtered/data_for_QC/qc_1/",
        sp,
        "_forQC_ssm.pdf",
        sep = ""
      ),
      width = 8,
      height = 10,
      pointsize = 16
    )
  }
  else if(filt_rnd == "1b"){
    pdf(
      file = paste(
        "../raatd_data/data_filtered/filtered_1/",
        sp,
        "_by_id_ssm.pdf",
        sep = ""
      ),
      width = 8,
      height = 10,
      pointsize = 16
    )
  }
  else if(filt_rnd == 2) {
    pdf(
      file = paste(
        "../raatd_data/data_filtered/data_for_QC/qc_2/",
        sp,
        "_stage_forQC_ssm.pdf",
        sep = ""
      ),
      width = 8,
      height = 10,
      pointsize = 16
    )
  }
  else if(filt_rnd == 3) {
    pdf(
      file = paste(
        "../raatd_data/data_filtered/filtered_by_stage/",
        sp,
        "_by_stage_ssm.pdf",
        sep = ""
      ),
      width = 8,
      height = 10,
      pointsize = 16
    )
  }

  ssm.tbl %>% do(p = plt.fn(.))
  dev.off()

}
