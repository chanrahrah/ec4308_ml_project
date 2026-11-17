# Portfolio engine shared by every model.

# Read monthly market returns and Treasury-bill yields.
load_return_data <- function() {
  market <- utils::read.csv(RETURN_FILES[1])
  market$date <- as.Date(market$mthcaldt)
  market$m <- month_number(market$date)
  market <- market[order(market$m), ]
  ensure(!anyDuplicated(market$m) && all(diff(market$m) == 1L), "Market calendar mismatch")
  bills <- utils::read.csv(RETURN_FILES[2])
  bills$m <- month_number(bills$MCALDT)
  bills <- bills[!is.na(bills$m), ]
  list(market = market, bills = bills)
}

# Align each return month with its h active signals and the bill yield each cohort locked in.
portfolio_inputs <- function(forecasts, returns, h, end, common_rf_security = 2000061L) {
  origin <- month_number(forecasts$origin)
  ensure(all(diff(origin) == 1L), "Forecasts must have a complete monthly calendar")
  b <- returns$bills[returns$bills$KYTREASNOX == c(`1` = 2000061, `3` = 2000062, `6` = 2000063)[[as.character(h)]], ]
  ensure(!anyDuplicated(b$m), "Duplicate matched bill quote")
  yield <- b$TMYLD[match(origin, b$m)]
  ensure(all(is.finite(yield)), "Missing matched bill yield")
  # The first h-1 forecasts only warm up the cohorts; no signals come from before the split.
  j <- which(seq_along(origin) >= h & origin + 1L <= month_number(end))
  ensure(length(j) > 1L, "Insufficient complete portfolio months")
  k <- match(origin[j] + 1L, returns$market$m)
  prev <- match(origin[j], returns$market$m)
  ensure(!anyNA(c(k, prev)), "Missing monthly return quote")
  days <- as.numeric(returns$market$date[k] - returns$market$date[prev])
  common <- returns$bills[returns$bills$KYTREASNOX == common_rf_security, ]
  ensure(!anyDuplicated(common$m), "Duplicate common bill quote")
  common_yield <- common$TMYLD[match(origin[j], common$m)]
  ensure(all(is.finite(common_yield)), "Missing common bill yield")
  # Accrue every bill yield over this month's quote interval, whatever its maturity.
  common_return <- expm1(common_yield * days)
  active <- t(matrix(unlist(lapply(j, function(i) seq.int(i - h + 1L, i))), nrow = h))
  bill_matrix <- expm1(sweep(matrix(yield[active], nrow = length(j), ncol = h), 1L, days, "*"))
  ensure(all(is.finite(bill_matrix)) && all(days > 0), "Invalid bill accrual")
  list(h = h, j = j, active = active, ym = month_string(origin[j] + 1L), days = days,
       market = returns$market$vwretd[k], bill_matrix = bill_matrix,
       bill_reference = rowMeans(bill_matrix), common_rf_return = common_return,
       quote_start = returns$market$date[prev], quote_end = returns$market$date[k],
       origin = forecasts$origin[j])
}

# Signed equity trade q that hits the target weight after paying cost c * |q|.
rebalance_trade <- function(wealth, stock_value, target_weight, cost_rate) {
  ensure(wealth > 0 && stock_value >= -1e-12 && stock_value <= wealth + 1e-12 &&
           target_weight >= 0 && target_weight <= 1 && cost_rate >= 0 && cost_rate < 1,
         "Invalid rebalancing inputs")
  gap <- target_weight * wealth - stock_value
  trade <- if (gap >= 0) gap / (1 + target_weight * cost_rate) else gap / (1 - target_weight * cost_rate)
  cost <- cost_rate * abs(trade)
  list(trade = trade, cost = cost, capital = wealth - cost)
}

# Rebalance to the target weights each month and compound net returns.
simulate_weights <- function(weights, market, bill_contribution, cost_rate, initial_cost = FALSE) {
  n <- length(weights)
  ensure(n == length(market) && n == length(bill_contribution) && n > 1L &&
           all(is.finite(c(weights, market, bill_contribution))) && all(market > -1),
         "Invalid portfolio vectors")
  wealth <- 1; stock <- 0
  out <- vector("list", n)
  for (i in seq_len(n)) {
    before <- wealth
    if (i == 1L && !initial_cost) {
      trade <- list(trade = 0, cost = 0, capital = wealth)
      stock <- weights[i] * wealth
    } else trade <- rebalance_trade(wealth, stock, weights[i], cost_rate)
    gross <- weights[i] * market[i] + bill_contribution[i]
    stock_start <- weights[i] * trade$capital
    stock <- stock_start * (1 + market[i])
    wealth <- trade$capital * (1 + gross)
    ensure(wealth > 0 && stock <= wealth + 1e-10, "Invalid portfolio value")
    out[[i]] <- data.frame(stock_weight = weights[i], pretrade_wealth = before,
      gross_return = gross, net_return = wealth / before - 1,
      signed_stock_trade = trade$trade, traded_fraction = abs(trade$trade) / before,
      cost = trade$cost, cost_fraction = trade$cost / before, wealth = wealth,
      stock_start = stock_start, stock_end = stock)
  }
  do.call(rbind, out)
}

# Each of the h active cohorts gets equal capital and keeps its signal for h months.
# fixed_weight replaces the signals with a constant stock weight (e.g. the 50/50 benchmark).
portfolio_path <- function(inputs, signals, cost_rate, cfg, fixed_weight = NULL) {
  ensure(length(signals) >= max(inputs$active) && all(signals %in% 0:1), "Invalid trading signals")
  active_signals <- matrix(signals[inputs$active], nrow = nrow(inputs$active), ncol = inputs$h)
  weights <- rowMeans(active_signals)
  bill <- rowMeans((1 - active_signals) * inputs$bill_matrix)
  if (!is.null(fixed_weight)) {
    weights <- rep(fixed_weight, length(weights))
    bill <- (1 - fixed_weight) * inputs$bill_reference
  }
  p <- simulate_weights(weights, inputs$market, bill, cost_rate, cfg$initial_allocation_cost)
  cbind(data.frame(month = inputs$ym, last_signal_origin = inputs$origin,
    quote_start = inputs$quote_start, quote_end = inputs$quote_end,
    quote_days = inputs$days, market_return = inputs$market,
    common_rf_return = inputs$common_rf_return,
    matched_ladder_return = inputs$bill_reference), p)
}

# Annualised Sharpe ratio of monthly excess returns (NA when the excess return is constant).
annual_sharpe <- function(excess) {
  ensure(length(excess) > 1L && all(is.finite(excess)), "Invalid excess-return sample")
  sigma <- stats::sd(excess)
  if (sigma > 1e-12) sqrt(12) * mean(excess) / sigma else NA_real_
}

# Returns, risk, equity exposure and trading costs over the months in p.
portfolio_metrics <- function(p, gamma = 3) {
  ensure(all(is.finite(as.matrix(p[, c("net_return", "common_rf_return", "matched_ladder_return")]))),
         "Missing or nonfinite returns")
  r <- p$net_return; excess <- r - p$matched_ladder_return
  sharpe <- annual_sharpe(excess)
  wealth <- cumprod(1 + r)
  dd <- wealth / cummax(c(1, wealth))[-1L] - 1
  maxdd <- max(0, -min(dd))
  cagr <- tail(wealth, 1)^(12 / length(r)) - 1
  downside <- sqrt(mean(pmin(excess, 0)^2))
  market_excess <- p$market_return - p$matched_ladder_return
  data.frame(n_months = length(r), start = p$month[1], end = tail(p$month, 1),
    total_return = tail(wealth, 1) - 1, CAGR = cagr,
    annual_volatility = stats::sd(r) * sqrt(12),
    sharpe_common_rf = annual_sharpe(r - p$common_rf_return), sharpe_matched_ladder = sharpe,
    Sortino = if (downside > 1e-12) sqrt(12) * mean(excess) / downside else NA_real_,
    max_drawdown = maxdd, Calmar = if (maxdd > 1e-12) cagr / maxdd else NA_real_,
    Ulcer_index = sqrt(mean(dd^2)),
    beta = stats::cov(excess, market_excess) / stats::var(market_excess),
    annual_CE_excess = 12 * (mean(excess) - gamma / 2 * stats::var(excess)),
    average_stock_weight = mean(p$stock_weight),
    annual_turnover = 12 * mean(p$traded_fraction), total_cost_fraction = sum(p$cost_fraction),
    rebalance_trade_months = sum(p$traded_fraction > 1e-12),
    allocation_change_months = sum(abs(diff(p$stock_weight)) > 1e-12))
}

# Highest validation Sharpe; ties go to the cutoff nearest 0.50, then the lower one.
choose_cutoff <- function(grid, objective = "sharpe_matched_ladder") {
  scores <- grid[[objective]]
  eligible <- which(is.finite(scores))
  ensure(length(eligible) > 0, "No cutoff has a defined validation Sharpe")
  best <- max(scores[eligible])
  tied <- eligible[abs(scores[eligible] - best) < 1e-12]
  tied[order(abs(grid$cutoff[tied] - 0.5), grid$cutoff[tied])][1L]
}
