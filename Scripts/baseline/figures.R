# Baseline figures: test wealth paths and probit test probabilities.

# Test wealth at 10 bp for tuned and 0.50 cutoffs against the passive benchmarks.
draw_wealth_figure <- function(out, paths, cfg) {
  colors <- c(probit = "#2463A7", logit = "#BA5A31", all_stocks = "#293A45", all_bills = "#78916D")
  grDevices::png(file.path(out, "test_wealth.png"), width = 1800, height = 780, res = 150)
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (h in cfg$horizons) {
    d <- paths[paths$period == "test" & paths$h == h & paths$cost_rate == cfg$primary_cost, ]
    dates <- as.Date(paste0(d$month, "-01"))
    graphics::plot(range(dates), range(c(1, d$wealth)), type = "n", xlab = "Return month",
                   ylab = "Wealth (initial = 1)", main = paste(h, "month horizon"))
    graphics::grid(col = "#E7E7E7")
    for (s in names(colors)) {
      v <- d[d$strategy == s & d$cutoff_rule %in% c("validation_tuned", "not_applicable"), ]
      graphics::lines(as.Date(paste0(v$month, "-01")), v$wealth, col = colors[s], lwd = 2)
      if (s %in% cfg$links) {
        v <- d[d$strategy == s & d$cutoff_rule == "fixed_050", ]
        graphics::lines(as.Date(paste0(v$month, "-01")), v$wealth, col = colors[s], lwd = 1.5, lty = 2)
      }
    }
    graphics::legend("topleft", legend = c("Probit tuned", "Logit tuned", "All stocks", "All bills", "Probit 0.50", "Logit 0.50"),
                     col = c(colors, colors[c("probit", "logit")]), lty = c(1, 1, 1, 1, 2, 2), lwd = 2, bty = "n", cex = .68)
    # Note when the tuned models never leave equities.
    m <- d[d$strategy %in% cfg$links & d$cutoff_rule == "validation_tuned", ]
    if (all(abs(m$stock_weight - 1) < 1e-12))
      graphics::mtext("Probit = Logit = All stocks", side = 3, line = .15, cex = .65, col = "#555555")
  }
  graphics::mtext("Test portfolios | 10 bp costs | validation-tuned and fixed 0.50 cutoffs", outer = TRUE, cex = 1.05)
  grDevices::dev.off()
}

# Probit test probabilities by origin, shaded where bills beat stocks, with the tuned and 0.50 cutoffs.
draw_probability_figure <- function(out, forecasts, cfg) {
  grDevices::png(file.path(out, "test_probabilities.png"), width = 7, height = 3.3, units = "in", res = 300, pointsize = 11)
  graphics::par(mfrow = c(3, 1), cex = 1, mar = c(1.0, 3.0, 1.0, .3), oma = c(1.1, .7, .2, 0), mgp = c(1.7, .45, 0), tcl = -.2)
  limits <- as.Date(paste0(c(cfg$test_start, month_string(month_number(cfg$test_end) + 1L)), "-01"))
  for (h in cfg$horizons) {
    f <- forecasts[forecasts$model == "probit" & forecasts$period == "test" & forecasts$h == h, ]
    start <- as.Date(paste0(f$origin, "-01"))
    end <- as.Date(paste0(month_string(month_number(f$origin) + 1L), "-01"))
    midpoint <- start + (end - start) / 2
    graphics::plot(midpoint, f$probability, type = "n", xlim = limits, ylim = c(.35, .9),
                   xaxs = "i", yaxs = "i", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
    below <- which(!is.na(f$target) & f$target == 0 & f$classification_eligible)
    graphics::rect(start[below], .35, end[below], .9, col = "#E7E7E7", border = NA)
    graphics::abline(h = f$cutoff[1], col = "#A45A35", lty = 2, lwd = 1)
    graphics::abline(h = .50, col = "#4F5964", lty = 3, lwd = 1)
    graphics::lines(midpoint, f$probability, col = "#1D5C83", lwd = 1.3)
    graphics::axis(2, at = c(.4, .6, .8), las = 1, cex.axis = .8)
    ticks <- as.Date(paste0(seq(2013, 2025, by = 2), "-01-01"))
    graphics::axis(1, at = ticks, labels = if (h == max(cfg$horizons)) format(ticks, "%Y") else FALSE, cex.axis = .8)
    graphics::mtext(paste(h, "month horizon"), side = 3, line = .1, adj = 0, cex = .8)
    if (h == min(cfg$horizons))
      graphics::legend("topright", legend = c("Probit", sprintf("%.2f", f$cutoff[1]), "0.50", "Bills outperform"),
                       col = c("#1D5C83", "#A45A35", "#4F5964", "#E7E7E7"), lty = c(1, 2, 3, NA),
                       pch = c(NA, NA, NA, 15), lwd = c(1.3, 1, 1, NA), bty = "n", horiz = TRUE, cex = .75, x.intersp = .6)
  }
  graphics::mtext("Probability", side = 2, outer = TRUE, line = -.2, cex = .85)
  graphics::mtext("Forecast origin", side = 1, outer = TRUE, line = 0, cex = .85)
  grDevices::dev.off()
}
