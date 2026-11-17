# Checks of the portfolio engine and evaluation rules, run before every model.
shared_checks <- function(returns, cfg) {
  checks <- list(
    "Six overlapping cohorts reproduce a hand-computed 12-month path" = function() {
      signals <- c(1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,1)
      market <- c(.03,-.02,.04,-.06,.02,.05,-.08,.03,-.01,.06,-.04,.02)
      w <- vapply(1:12, function(i) mean(signals[i:(i+5)]), numeric(1))
      p <- simulate_weights(w, market, (1-w)*.003, 0)
      value <- 100; expected <- numeric(12)
      for (i in 1:12) { value <- w[i]*value*(1+market[i]) + (1-w[i])*value*1.003; expected[i] <- value }
      near(p$wealth*100, expected)
      near(p$signed_stock_trade[-1]*100, w[-1]*expected[-12] - w[-12]*c(100,expected[1:10])*(1+market[-12]))
    },
    "Trades hit the target weight exactly after costs" = function() {
      for (s in c(0, 20, 50, 100)) for (w in c(0, .2, .5, 1)) for (c in c(0, .001, .002)) {
        t <- rebalance_trade(100, s, w, c)
        near((s + t$trade) / t$capital, w); near(t$cost, abs(t$trade) * c); near(t$capital, 100 - t$cost)
      }
    },
    "A constant target still rebalances after market drift" = function() {
      p <- simulate_weights(rep(.5, 3), c(.1, -.05, .02), rep(.0015, 3), .001)
      ensure(p$signed_stock_trade[2] < 0 && p$cost[2] > 0, "Drift trade omitted")
    },
    "Cutoff ties go to 0.50 and undefined Sharpes are skipped" = function() {
      g <- data.frame(cutoff = c(.2, .49, .51, 1), sharpe_matched_ladder = c(.1, .2, .2, NA))
      ensure(choose_cutoff(g) == 2L, "Cutoff tie rule failed")
      g$sharpe_matched_ladder <- NA_real_
      ensure(inherits(try(choose_cutoff(g), silent = TRUE), "try-error"), "Undefined Sharpes must stop selection")
    },
    "Brier ties go to the first candidate in the tie order" = function() {
      ensure(best_candidate(c(.20, .19, .19 + 1e-13, .25), c(1, 3, 2, 4)) == 3L, "Brier tie rule failed")
      ensure(best_candidate(c(.20, .18, .19, .25), c(1, 3, 2, 4)) == 2L, "Lower Brier must win")
    },
    "Both Sharpe ratios use their own reference and the sample SD" = function() {
      p <- data.frame(month = paste0("2000-0", 1:4), net_return = c(.03, -.01, .02, .04),
        market_return = c(.04, -.02, .03, .05), common_rf_return = c(.001, .002, .001, .002),
        matched_ladder_return = c(.002, .004, .003, .005), traded_fraction = 0, cost_fraction = 0, stock_weight = .5)
      m <- portfolio_metrics(p)
      for (ref in c("common_rf", "matched_ladder")) {
        e <- p$net_return - p[[paste0(ref, "_return")]]
        near(m[[paste0("sharpe_", ref)]], sqrt(12) * mean(e) / stats::sd(e))
      }
    },
    "AUC handles ties and the cutoff is inclusive" = function() {
      m <- forecast_metrics(data.frame(classification_eligible = TRUE, target = c(0, 1, 0, 1), probability = .5), .5)
      near(m$AUC, .5); near(m$Brier, .25); near(m$recall, 1)
      near(auc(c(0, 0, 1, 1), c(.1, .2, .8, .9)), 1); near(auc(c(0, 0, 1, 1), c(.9, .8, .2, .1)), 0)
    },
    "The common reference is identical across horizons" = function() {
      refs <- lapply(cfg$horizons, function(h) {
        origins <- month_string(month_number(cfg$validation_start) + 0:95)
        inp <- portfolio_inputs(data.frame(origin = origins), returns, h, cfg$validation_end, cfg$common_rf_security)
        setNames(inp$common_rf_return, inp$ym)
      })
      months <- Reduce(intersect, lapply(refs, names))
      for (ref in refs[-1]) near(ref[months], refs[[1]][months])
    })
  for (h in cfg$horizons) checks[[paste("Bill accrual and passive strategies h", h)]] <- local({ h <- h; function() {
    d <- load_model_data(h, cfg)
    origins <- d$ym[period_rows(d, cfg, "validation")]
    inp <- portfolio_inputs(data.frame(origin = origins), returns, h, cfg$validation_end, cfg$common_rf_security)
    ensure(nrow(inp$active) == length(origins) - h && inp$ym[1] == month_string(month_number(cfg$validation_start) + h),
           "Wrong warm-up")
    ensure(all(inp$days >= 25 & inp$days <= 35), "Bills must accrue over one month, not their maturity")
    near(inp$days, as.numeric(inp$quote_end - inp$quote_start))
    if (h == 1) near(inp$bill_reference, inp$common_rf_return)
    # One vintage compounded over its h months matches the target's bill accrual.
    origin <- month_number("2010-01"); mi <- match(origin, returns$market$m)
    b <- returns$bills[returns$bills$KYTREASNOX == c(`1` = 2000061, `3` = 2000062, `6` = 2000063)[[as.character(h)]], ]
    days <- as.numeric(diff(returns$market$date[mi:(mi + h)]))
    near(prod(exp(b$TMYLD[match(origin, b$m)] * days)) - 1, d$metadata$bill_accrual_proxy[d$m == origin])
    stocks <- portfolio_path(inp, rep(1L, length(origins)), .002, cfg)
    bills <- portfolio_path(inp, rep(0L, length(origins)), .002, cfg)
    mix <- portfolio_path(inp, rep(0L, length(origins)), .002, cfg, fixed_weight = .5)
    near(stocks$net_return, inp$market); near(bills$net_return, inp$bill_reference)
    ensure(all(c(stocks$cost, bills$cost) == 0) && all(mix$stock_weight == .5), "Passive strategies must not trade")
    ensure(is.na(portfolio_metrics(bills)$sharpe_matched_ladder), "All-bills matched Sharpe must be undefined")
    near(portfolio_metrics(stocks)$beta, 1)
  }})
  checks
}

# A forecast must not change when anything unknown at its origin is altered.
# predict_at(d, i) returns the model's probability for row i using only mature rows.
perturbation_check <- function(d, h, predict_at, origins = c("2005-01", "2013-01"), tol = 1e-10) {
  for (o in origins) {
    origin <- month_number(o); i <- which(d$m == origin)
    rows <- mature_rows(d, origin)
    ensure(max(d$m[rows]) == origin - h && max(d$end[rows]) == origin, "Wrong maturity cutoff")
    changed <- d
    unknown <- which(is.na(d$end) | d$end > origin)
    changed$y[unknown] <- 1L - ifelse(is.na(d$y[unknown]), 0L, d$y[unknown])
    changed$X[d$m > origin, ] <- changed$X[d$m > origin, ] + 1000
    a <- predict_at(d, i); b <- predict_at(changed, i)
    if (tol == 0) ensure(identical(a, b), "Future data changed a forecast") else near(a, b, tol)
  }
}
