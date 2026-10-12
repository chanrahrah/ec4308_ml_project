"""Build monthly CRSP predictors and month-end lagged features."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "Data"

STOCK_FILE = DATA_DIR / "crsp_monthly_stock_data" / "crsp_monthly_stock_data.csv"
MONTHLY_INDEX_FILE = DATA_DIR / "crsp_monthly_index_data" / "crsp_monthly_index_data.csv"
DAILY_INDEX_FILE = DATA_DIR / "crsp_daily_index_data" / "crsp_daily_index_data.csv"
TREASURY_DAILY_FILE = (
    DATA_DIR / "crsp_treasury_risk_free_series_data"
    / "crsp_treasury_risk_free_series_data.csv"
)
TREASURY_FIXED_FILE = (
    DATA_DIR / "crsp_treasury_fixed_term_index_data"
    / "crsp_treasury_fixed_term_index_data.csv"
)
CPI_FILE = DATA_DIR / "cpi_u_data" / "CPIAUCNS.csv"

OUTPUT_DIR = SCRIPT_DIR / "output"
CHECKS_DIR = OUTPUT_DIR / "checks"
STOCK_CACHE = OUTPUT_DIR / "cache" / "crsp_monthly_stock_selected.parquet"
OUTPUT_START = pd.Period("1970-01", freq="M")
FEATURE_START = pd.Period("1973-01", freq="M")

BASE_PREDICTORS = [
    "mkt_excess",
    "mom3",
    "mom6",
    "mom12",
    "log_rv",
    "downside_share",
    "log_dp",
    "bill_yld_h1",
    "bill_yld_h3",
    "bill_yld_h6",
    "bill_yld_chg",
    "term_spread",
    "short_slope",
    "cpi_infl",
    "breadth",
    "ma_part",
    "hi_lo",
    "dispersion",
    "ew_vw",
    "tail_imb",
    "turnover_dt",
    "illiq_dt",
]

STOCK_SOURCE_COLUMNS = [
    "PERMNO",
    "MthCalDt",
    "PrimaryExch",
    "ConditionalType",
    "TradingStatusFlg",
    "USIncFlg",
    "IssuerType",
    "SecurityType",
    "SecuritySubType",
    "ShareType",
    "MthRet",
    "MthRetx",
    "MthPrc",
    "MthVol",
    "ShrOut",
]
STOCK_COLUMN_NAMES = [column.lower() for column in STOCK_SOURCE_COLUMNS]
STOCK_KEY = ["permno", "ym"]
STOCK_NUMERIC_FOR_CONFLICT_CHECK = [
    "mthret",
    "mthretx",
    "mthprc",
    "mthvol",
    "shrout",
]


def log(message: str) -> None:
    """Print a concise progress message."""
    print(f"[build_predictors] {message}", flush=True)


def require_columns(frame: pd.DataFrame, columns: list[str], source: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def add_month_period(frame: pd.DataFrame, date_column: str) -> pd.DataFrame:
    frame = frame.copy()
    frame[date_column] = pd.to_datetime(frame[date_column], errors="raise")
    frame["ym"] = frame[date_column].dt.to_period("M")
    return frame


def load_stock_data() -> tuple[pd.DataFrame, dict[str, int]]:
    """Load selected CIZ fields, cache the source subset, and resolve duplicate keys."""
    STOCK_CACHE.parent.mkdir(parents=True, exist_ok=True)
    if STOCK_CACHE.exists():
        log(f"Reading stock cache: {STOCK_CACHE}")
        stock = pd.read_parquet(STOCK_CACHE)
        require_columns(stock, STOCK_COLUMN_NAMES, "Stock parquet cache")
        stock["mthcaldt"] = pd.to_datetime(stock["mthcaldt"], errors="raise")
    else:
        log(f"Reading selected stock columns from CSV: {STOCK_FILE}")
        if not STOCK_FILE.exists():
            raise FileNotFoundError(f"Stock CSV not found: {STOCK_FILE}")
        stock_dtypes: dict[str, Any] = {
            "PERMNO": np.int32,
            "MthCalDt": "string",
            "PrimaryExch": "category",
            "ConditionalType": "category",
            "TradingStatusFlg": "category",
            "USIncFlg": "category",
            "IssuerType": "category",
            "SecurityType": "category",
            "SecuritySubType": "category",
            "ShareType": "category",
            "MthRet": np.float64,
            "MthRetx": np.float64,
            "MthPrc": np.float64,
            "MthVol": np.float64,
            "ShrOut": np.float64,
        }
        stock = pd.read_csv(
            STOCK_FILE,
            usecols=STOCK_SOURCE_COLUMNS,
            dtype=stock_dtypes,
            low_memory=False,
        )
        stock.columns = stock.columns.str.lower()
        stock["mthcaldt"] = pd.to_datetime(stock["mthcaldt"], errors="raise")
        # The cache preserves source duplicates so their resolution is auditable on every run.
        stock.to_parquet(STOCK_CACHE, index=False)
        log(f"Wrote selected stock parquet cache: {STOCK_CACHE}")

    stock["ym"] = stock["mthcaldt"].dt.to_period("M")
    require_columns(stock, STOCK_COLUMN_NAMES + ["ym"], "Stock data")

    raw_duplicate_groups = stock.groupby(STOCK_KEY, observed=True, sort=False).size()
    raw_duplicate_keys = int((raw_duplicate_groups > 1).sum())

    before_exact = len(stock)
    stock = stock.drop_duplicates(subset=STOCK_COLUMN_NAMES).copy()
    exact_duplicate_rows_removed = before_exact - len(stock)

    duplicate_mask = stock.duplicated(STOCK_KEY, keep=False)
    duplicate_rows = stock.loc[duplicate_mask]
    duplicate_groups = duplicate_rows.groupby(STOCK_KEY, observed=True, sort=False)
    duplicate_keys_after_exact = int(duplicate_groups.ngroups)
    conflicting = duplicate_groups[STOCK_NUMERIC_FOR_CONFLICT_CHECK].nunique(
        dropna=False
    )
    conflict_keys = conflicting.index[(conflicting > 1).any(axis=1)]
    if len(conflict_keys):
        conflict_table = conflict_keys.to_frame(index=False)
        print("Conflicting stock keys (numeric fields disagree):")
        print(conflict_table.to_string(index=False))
        raise ValueError(
            f"{len(conflict_keys)} PERMNO/month keys disagree on stock numeric fields."
        )

    stock["_in_universe"] = (
        stock["sharetype"].eq("NS")
        & stock["securitytype"].eq("EQTY")
        & stock["securitysubtype"].eq("COM")
        & stock["usincflg"].eq("Y")
        & stock["issuertype"].isin(["ACOR", "CORP"])
        & stock["primaryexch"].isin(["N", "A", "Q"])
        & stock["conditionaltype"].eq("RW")
        & stock["tradingstatusflg"].eq("A")
    )

    # If one copy qualifies for the universe, retain that copy's flags.
    stock = stock.sort_values(
        STOCK_KEY + ["_in_universe"],
        ascending=[True, True, False],
        kind="mergesort",
    )
    before_key_collapse = len(stock)
    stock = stock.drop_duplicates(STOCK_KEY, keep="first").reset_index(drop=True)
    collapsed_rows = before_key_collapse - len(stock)

    duplicate_info = {
        "raw_duplicate_keys": raw_duplicate_keys,
        "exact_duplicate_rows_removed": exact_duplicate_rows_removed,
        "duplicate_keys_collapsed_after_exact_dedup": duplicate_keys_after_exact,
        "additional_rows_collapsed_by_key": collapsed_rows,
        "numeric_conflict_keys": int(len(conflict_keys)),
    }
    log(
        "Stock duplicate handling: "
        f"{raw_duplicate_keys:,} raw duplicate keys; "
        f"{exact_duplicate_rows_removed:,} exact duplicate rows removed; "
        f"{duplicate_keys_after_exact:,} keys collapsed after exact deduplication; "
        f"{collapsed_rows:,} additional rows removed."
    )
    log(f"Loaded {len(stock):,} unique PERMNO/month stock rows.")
    return stock, duplicate_info


def load_monthly_index() -> pd.DataFrame:
    log("Reading monthly CRSP index")
    frame = pd.read_csv(
        MONTHLY_INDEX_FILE,
        usecols=["mthcaldt", "vwretd", "vwretx", "ewretd", "vwtotval"],
    )
    frame = add_month_period(frame, "mthcaldt")
    if frame["ym"].duplicated().any():
        raise ValueError("Monthly index has more than one row for a month.")
    return frame


def load_daily_index() -> pd.DataFrame:
    log("Reading daily CRSP value-weighted returns")
    frame = pd.read_csv(DAILY_INDEX_FILE, usecols=["dlycaldt", "vwretd"])
    frame = add_month_period(frame, "dlycaldt")
    return frame


def load_treasury_and_cpi() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    log("Reading monthly Treasury bill series, 10-year series, and CPI")
    bills = pd.read_csv(
        TREASURY_DAILY_FILE,
        usecols=["KYTREASNOX", "TTERMLBL", "MCALDT", "TMYLD"],
    )
    bills.columns = bills.columns.str.lower()
    bills = bills.loc[bills["kytreasnox"].isin([2000061, 2000062, 2000063])].copy()
    bills = bills.dropna(subset=["mcaldt", "tmyld"])
    bills = add_month_period(bills, "mcaldt")

    fixed = pd.read_csv(
        TREASURY_FIXED_FILE,
        usecols=["KYTREASNOX", "TTERMLBL", "MCALDT", "TMYEARSTM", "TMYTM"],
    )
    fixed.columns = fixed.columns.str.lower()
    fixed = fixed.loc[fixed["kytreasnox"].eq(2000007)].copy()
    fixed = add_month_period(fixed, "mcaldt")

    cpi = pd.read_csv(CPI_FILE, usecols=["observation_date", "CPIAUCNS"])
    cpi.columns = cpi.columns.str.lower()
    cpi = add_month_period(cpi, "observation_date")
    return bills, fixed, cpi


def make_rate_inflation_predictors(
    bills: pd.DataFrame,
    fixed: pd.DataFrame,
    cpi: pd.DataFrame,
    calendar: pd.PeriodIndex,
) -> tuple[pd.DataFrame, pd.Series, dict[str, str]]:
    """Build bill, term, slope, and publication-lagged inflation predictors."""
    expected_bill_codes = {2000061, 2000062, 2000063}
    found_bill_codes = set(bills["kytreasnox"].dropna().astype(int).unique())
    if not expected_bill_codes.issubset(found_bill_codes):
        raise ValueError(
            f"Treasury bill series missing; found codes {sorted(found_bill_codes)}."
        )
    if fixed["kytreasnox"].nunique() != 1 or not fixed["kytreasnox"].eq(2000007).all():
        raise ValueError("Expected only the CRSP 10-year fixed-term Treasury series.")

    bill_raw_monthly = (
        bills.sort_values("mcaldt")
        .groupby(["kytreasnox", "ym"], observed=True, sort=False)
        .tail(1)
        .pivot(index="ym", columns="kytreasnox", values="tmyld")
        .reindex(calendar)
    )
    # TMYLD is a daily promised yield; this compounds that daily rate for 365 days.
    bill_monthly = (1.0 + bill_raw_monthly.astype(float)) ** 365.0 - 1.0
    bill_monthly = bill_monthly.rename(
        columns={
            2000061: "bill_yld_h1",
            2000062: "bill_yld_h3",
            2000063: "bill_yld_h6",
        }
    )
    daily_bill_h1 = bill_raw_monthly[2000061].rename("bill_daily_h1")

    # TMYTM is a percent-scale bond-equivalent yield; convert it to effective annual.
    ytm_decimal = fixed["tmytm"].astype(float) / 100.0
    fixed = fixed.assign(
        effective_annual_yield=(1.0 + ytm_decimal / 2.0) ** 2.0 - 1.0
    )
    ten_year = (
        fixed.sort_values("mcaldt")
        .groupby("ym", observed=True, sort=False)
        .tail(1)
        .set_index("ym")["effective_annual_yield"]
        .reindex(calendar)
        .rename("yield_10y")
    )

    result = bill_monthly.copy()
    result["bill_yld_chg"] = result["bill_yld_h3"].diff()
    result["term_spread"] = ten_year - result["bill_yld_h3"]
    result["short_slope"] = result["bill_yld_h6"] - result["bill_yld_h1"]

    cpi_values = cpi.set_index("ym")["cpiaucns"].reindex(calendar).astype(float)
    cpi_infl = np.log(cpi_values.shift(1) / cpi_values.shift(13))
    flags: dict[str, str] = {}
    missing_cpi_months = cpi.loc[cpi["cpiaucns"].isna(), "ym"].astype(str).tolist()
    missing_period = pd.Period("2025-10", freq="M")
    prior_infl_period = pd.Period("2025-10", freq="M")
    override_period = pd.Period("2025-11", freq="M")
    if missing_period in set(cpi.loc[cpi["cpiaucns"].isna(), "ym"]):
        # October CPI was unpublished; at November's decision point use the latest
        # inflation reading actually available (the October feature).
        if prior_infl_period in calendar and override_period in calendar:
            cpi_infl.loc[override_period] = cpi_infl.loc[prior_infl_period]
            flags["cpi_publication_gap"] = (
                "CPI 2025-10 is missing. Per instruction, cpi_infl 2025-11 "
                "was set equal to cpi_infl 2025-10; no interpolation was used."
            )
    if missing_cpi_months:
        flags["cpi_missing_months"] = ", ".join(missing_cpi_months)
    result["cpi_infl"] = cpi_infl
    return result, daily_bill_h1, flags


def make_market_predictors(
    monthly: pd.DataFrame,
    daily: pd.DataFrame,
    daily_bill_h1: pd.Series,
    calendar: pd.PeriodIndex,
) -> pd.DataFrame:
    """Build market trend, realized risk, dividend-price, and excess-return signals."""
    monthly_aligned = monthly.set_index("ym").reindex(calendar)
    result = pd.DataFrame(index=calendar)
    market_return = monthly_aligned["vwretd"].astype(float)

    for months, name in ((3, "mom3"), (6, "mom6"), (12, "mom12")):
        growth = (1.0 + market_return).rolling(months, min_periods=months).apply(
            np.prod, raw=True
        )
        result[name] = growth - 1.0

    daily = daily.copy()
    daily["squared_return"] = daily["vwretd"].astype(float).pow(2)
    daily["downside_squared_return"] = daily["squared_return"].where(
        daily["vwretd"] < 0.0, 0.0
    )
    realized = daily.groupby("ym", observed=True).agg(
        n_days=("vwretd", "count"),
        sum_squared=("squared_return", "sum"),
        downside_squared=("downside_squared_return", "sum"),
    ).reindex(calendar)
    annualized_rv = np.sqrt(
        (252.0 / realized["n_days"]) * realized["sum_squared"]
    )
    result["log_rv"] = np.log(annualized_rv.where(annualized_rv > 0.0))
    result["downside_share"] = (
        realized["downside_squared"] / realized["sum_squared"]
    ).where(realized["sum_squared"] > 0.0)

    # Use the prior month's total market value as the dividend-dollar base, as confirmed.
    dividend_dollars = (
        (monthly_aligned["vwretd"] - monthly_aligned["vwretx"])
        * monthly_aligned["vwtotval"].shift(1)
    )
    dividend_sum_12 = dividend_dollars.rolling(12, min_periods=12).sum()
    dividend_yield = dividend_sum_12 / monthly_aligned["vwtotval"]
    result["log_dp"] = np.log(dividend_yield.where(dividend_yield > 0.0))

    days_in_month = pd.Series(calendar.days_in_month, index=calendar)
    # The month-t bill return compounds the raw daily TMYLD observed at t-1
    # over the actual number of calendar days in month t.
    bill_return = (1.0 + daily_bill_h1.shift(1)).pow(days_in_month) - 1.0
    result["mkt_excess"] = market_return - bill_return
    return result


def build_price_window_signals(stock: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate valid MA/high/low price-index signals on complete PERMNO calendars."""
    ma_signal = np.full(len(stock), np.nan, dtype=np.float32)
    high_signal = np.full(len(stock), np.nan, dtype=np.float32)
    low_signal = np.full(len(stock), np.nan, dtype=np.float32)

    for _, group in stock.groupby("permno", sort=False, observed=True):
        # Rows are sorted by PERMNO/month by load_stock_data. Expand only this
        # PERMNO's observed span; absent months become NaN and invalidate windows.
        month_ordinals = group["ym"].array.asi8
        first_month = month_ordinals[0]
        positions = month_ordinals - first_month
        returns = np.full(month_ordinals[-1] - first_month + 1, np.nan)
        returns[positions] = group["mthretx"].to_numpy(dtype=np.float64)
        row_indices = group.index.to_numpy(dtype=np.int64)

        if len(returns) >= 9:
            windows = np.lib.stride_tricks.sliding_window_view(returns, 9)
            factors = 1.0 + windows
            paths = np.concatenate(
                [np.ones((len(windows), 1)), np.cumprod(factors, axis=1)],
                axis=1,
            )
            valid = (
                np.isfinite(windows).all(axis=1)
                & (factors >= 0.0).all(axis=1)
                & np.isfinite(paths).all(axis=1)
            )
            ending_signal = np.full(len(windows), np.nan, dtype=np.float32)
            ending_signal[valid] = (
                paths[valid, -1] > paths[valid].mean(axis=1)
            ).astype(np.float32)
            full_signal = np.full(len(returns), np.nan, dtype=np.float32)
            full_signal[8:] = ending_signal
            ma_signal[row_indices] = full_signal[positions]

        if len(returns) >= 11:
            windows = np.lib.stride_tricks.sliding_window_view(returns, 11)
            factors = 1.0 + windows
            paths = np.concatenate(
                [np.ones((len(windows), 1)), np.cumprod(factors, axis=1)],
                axis=1,
            )
            valid = (
                np.isfinite(windows).all(axis=1)
                & (factors >= 0.0).all(axis=1)
                & np.isfinite(paths).all(axis=1)
            )
            end_is_high = np.full(len(windows), np.nan, dtype=np.float32)
            end_is_low = np.full(len(windows), np.nan, dtype=np.float32)
            end_is_high[valid] = (
                paths[valid, -1] == paths[valid].max(axis=1)
            ).astype(np.float32)
            end_is_low[valid] = (
                paths[valid, -1] == paths[valid].min(axis=1)
            ).astype(np.float32)
            full_high = np.full(len(returns), np.nan, dtype=np.float32)
            full_low = np.full(len(returns), np.nan, dtype=np.float32)
            full_high[10:] = end_is_high
            full_low[10:] = end_is_low
            high_signal[row_indices] = full_high[positions]
            low_signal[row_indices] = full_low[positions]

    return ma_signal, high_signal, low_signal


def make_participation_predictors(
    stock: pd.DataFrame,
    calendar: pd.PeriodIndex,
) -> pd.DataFrame:
    """Build breadth and price-index participation shares over the month-t universe."""
    log("Building per-PERMNO complete-calendar price-index windows")
    ma_signal, high_signal, low_signal = build_price_window_signals(stock)
    work = stock.copy()
    work["ma_above"] = ma_signal
    work["at_12m_high"] = high_signal
    work["at_12m_low"] = low_signal

    eligible = work.loc[work["_in_universe"] & work["mthret"].notna()].copy()
    if eligible.empty:
        raise ValueError("No nonmissing stock returns pass the requested universe.")
    eligible["positive_return"] = eligible["mthret"].gt(0.0).astype(float)
    grouped = eligible.groupby("ym", observed=True)
    result = pd.DataFrame(index=calendar)
    result["breadth"] = grouped["positive_return"].mean().reindex(calendar)
    result["ma_part"] = grouped["ma_above"].mean().reindex(calendar)
    result["hi_lo"] = (
        grouped["at_12m_high"].mean() - grouped["at_12m_low"].mean()
    ).reindex(calendar)
    return result


def make_cross_section_predictors(
    stock: pd.DataFrame,
    monthly: pd.DataFrame,
    calendar: pd.PeriodIndex,
) -> pd.DataFrame:
    """Build winsorized return dispersion, equal-minus-value return, and tail imbalance."""
    eligible = stock.loc[
        stock["_in_universe"] & stock["mthret"].notna(), ["ym", "mthret"]
    ].copy()
    result = pd.DataFrame(index=calendar)
    dispersion: dict[pd.Period, float] = {}
    tail_imbalance: dict[pd.Period, float] = {}
    for ym, group in eligible.groupby("ym", observed=True, sort=False):
        returns = group["mthret"].astype(float)
        lower, upper = returns.quantile([0.01, 0.99]).to_numpy()
        dispersion[ym] = float(returns.clip(lower, upper).std(ddof=1))
        tail_imbalance[ym] = float(
            returns.lt(-0.20).mean() - returns.gt(0.20).mean()
        )
    result["dispersion"] = pd.Series(dispersion).reindex(calendar)
    result["tail_imb"] = pd.Series(tail_imbalance).reindex(calendar)

    market = monthly.set_index("ym").reindex(calendar)
    result["ew_vw"] = market["ewretd"] - market["vwretd"]
    return result


def make_trading_activity_predictors(
    stock: pd.DataFrame,
    calendar: pd.PeriodIndex,
) -> pd.DataFrame:
    """Build NYSE/AMEX turnover and Amihud illiquidity changes."""
    nyse_amex = stock.loc[
        stock["_in_universe"] & stock["primaryexch"].isin(["N", "A"])
    ].copy()
    result = pd.DataFrame(index=calendar)

    valid_turnover = nyse_amex.loc[
        nyse_amex["mthvol"].notna()
        & nyse_amex["shrout"].notna()
        & nyse_amex["shrout"].gt(0.0)
    ]
    totals = valid_turnover.groupby("ym", observed=True).agg(
        total_volume=("mthvol", "sum"),
        total_shares=("shrout", "sum"),
    ).reindex(calendar)
    turnover = (totals["total_volume"] / totals["total_shares"]).where(
        totals["total_shares"] > 0.0
    )
    prior_turnover_mean = turnover.rolling(12, min_periods=12).mean().shift(1)
    result["turnover_dt"] = np.log(turnover.where(turnover > 0.0)) - np.log(
        prior_turnover_mean.where(prior_turnover_mean > 0.0)
    )

    valid_illiq = nyse_amex.loc[
        nyse_amex["mthret"].notna()
        & nyse_amex["mthvol"].gt(0.0)
        & nyse_amex["mthprc"].abs().ge(5.0)
    ].copy()
    valid_illiq["illiq"] = (
        valid_illiq["mthret"].abs()
        / (valid_illiq["mthprc"].abs() * valid_illiq["mthvol"])
    )
    monthly_illiq: dict[pd.Period, float] = {}
    for ym, group in valid_illiq.groupby("ym", observed=True, sort=False):
        values = group["illiq"].replace([np.inf, -np.inf], np.nan).dropna()
        if not values.empty:
            lower, upper = values.quantile([0.01, 0.99]).to_numpy()
            monthly_illiq[ym] = float(values.clip(lower, upper).mean())
    illiq = pd.Series(monthly_illiq).reindex(calendar)
    prior_illiq_mean = illiq.rolling(12, min_periods=12).mean().shift(1)
    result["illiq_dt"] = np.log(illiq.where(illiq > 0.0)) - np.log(
        prior_illiq_mean.where(prior_illiq_mean > 0.0)
    )
    return result


def make_predictor_frame(
    groups: list[pd.DataFrame],
    calendar: pd.PeriodIndex,
) -> pd.DataFrame:
    base = pd.DataFrame(index=calendar)
    for group in groups:
        for column in group.columns:
            if column in base.columns:
                raise ValueError(f"Predictor {column!r} was created more than once.")
        aligned = group.reindex(calendar)
        for column in aligned.columns:
            base[column] = aligned[column]
    missing = sorted(set(BASE_PREDICTORS).difference(base.columns))
    extra = sorted(set(base.columns).difference(BASE_PREDICTORS))
    if missing or extra:
        raise ValueError(f"Base predictor mismatch; missing={missing}, extra={extra}")
    return base[BASE_PREDICTORS]


def apply_lags(base: pd.DataFrame) -> pd.DataFrame:
    """Apply month-index lags on the complete monthly calendar."""
    lagged = pd.DataFrame(index=base.index)
    for predictor in BASE_PREDICTORS:
        for lag in range(6):
            lagged[f"{predictor}_lag{lag}"] = base[predictor].shift(lag)
    lagged = lagged.loc[lagged.index >= FEATURE_START]
    lagged.index.name = "ym"
    return lagged.reset_index()


def write_outputs(base: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CHECKS_DIR.mkdir(parents=True, exist_ok=True)

    base_file = OUTPUT_DIR / "base_predictors_monthly.csv"
    feature_file = OUTPUT_DIR / "features_lagged.csv"
    base_out = base.reset_index(names="ym").copy()
    base_out["ym"] = base_out["ym"].astype(str)
    base_out.to_csv(base_file, index=False)

    lagged = apply_lags(base)
    lagged["ym"] = lagged["ym"].astype(str)
    lagged.to_csv(feature_file, index=False)
    log(
        f"Wrote {base_file.name} ({len(base_out):,} rows, "
        f"{len(BASE_PREDICTORS)} predictors) and {feature_file.name} "
        f"({len(lagged):,} rows, {len(lagged.columns) - 1} features)."
    )
    return base_out, lagged


def write_checks(
    base: pd.DataFrame,
    duplicate_info: dict[str, int],
    cpi_flags: dict[str, str],
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    messages: list[str] = []
    summary_records = []
    for predictor in BASE_PREDICTORS:
        values = base[predictor]
        valid = values.dropna()
        summary_records.append(
            {
                "predictor": predictor,
                "mean": valid.mean() if len(valid) else np.nan,
                "sd": valid.std(ddof=1) if len(valid) > 1 else np.nan,
                "min": valid.min() if len(valid) else np.nan,
                "max": valid.max() if len(valid) else np.nan,
                "first_nonmissing_ym": str(valid.index[0]) if len(valid) else "",
                "nan_count": int(values.isna().sum()),
            }
        )
        fig, ax = plt.subplots(figsize=(11, 3.2))
        ax.plot(base.index.to_timestamp(), values.to_numpy(), linewidth=0.85)
        ax.set_title(predictor)
        ax.set_xlabel("Month")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(CHECKS_DIR / f"{predictor}.png", dpi=130)
        plt.close(fig)

    summary = pd.DataFrame(summary_records)
    summary.to_csv(CHECKS_DIR / "predictor_summary.csv", index=False)

    correlation = base.corr(min_periods=2)
    fig, ax = plt.subplots(figsize=(14, 11))
    image = ax.imshow(correlation.to_numpy(), cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(BASE_PREDICTORS)), BASE_PREDICTORS, rotation=90)
    ax.set_yticks(np.arange(len(BASE_PREDICTORS)), BASE_PREDICTORS)
    ax.set_title("Lag-0 predictor correlation")
    fig.colorbar(image, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(CHECKS_DIR / "correlation_lag0.png", dpi=150)
    plt.close(fig)

    def record(label: str, passed: bool, detail: str) -> None:
        status = "PASS" if passed else "FLAG"
        message = f"{status}: {label}: {detail}"
        messages.append(message)
        print(message)

    bounded_01 = ["downside_share", "breadth", "ma_part"]
    for column in bounded_01:
        valid = base[column].dropna()
        ok = bool(valid.between(0.0, 1.0).all())
        record(
            f"{column} bounds [0, 1]",
            ok,
            f"observed [{valid.min():.6g}, {valid.max():.6g}]" if len(valid) else "all NaN",
        )
    for column in ["tail_imb", "hi_lo"]:
        valid = base[column].dropna()
        ok = bool(valid.between(-1.0, 1.0).all())
        record(
            f"{column} bounds [-1, 1]",
            ok,
            f"observed [{valid.min():.6g}, {valid.max():.6g}]" if len(valid) else "all NaN",
        )

    mean_excess = base["mkt_excess"].mean()
    record(
        "Mean monthly market excess return approximately 0.4%-0.8%",
        bool(pd.notna(mean_excess) and 0.004 <= mean_excess <= 0.008),
        f"mean={mean_excess:.6%}" if pd.notna(mean_excess) else "mean is NaN",
    )

    annualized_rv = np.exp(base["log_rv"])
    typical_rv = annualized_rv.median()
    record(
        "Typical annualized realized volatility approximately 10%-25%",
        bool(pd.notna(typical_rv) and 0.10 <= typical_rv <= 0.25),
        f"median={typical_rv:.2%}" if pd.notna(typical_rv) else "median is NaN",
    )
    crash_months = [pd.Period(x, freq="M") for x in ("1987-10", "2008-10", "2020-03")]
    for month in crash_months:
        rv = annualized_rv.get(month, np.nan)
        record(
            f"Realized volatility spike {month}",
            bool(pd.notna(rv) and pd.notna(typical_rv) and rv > typical_rv),
            f"annualized RV={rv:.2%}" if pd.notna(rv) else "no value",
        )

    inversion_windows = {
        "1979-1981": (pd.Period("1979-01", "M"), pd.Period("1981-12", "M")),
        "2000": (pd.Period("2000-01", "M"), pd.Period("2000-12", "M")),
        "2006-2007": (pd.Period("2006-01", "M"), pd.Period("2007-12", "M")),
        "2023": (pd.Period("2023-01", "M"), pd.Period("2023-12", "M")),
    }
    for label, (start, end) in inversion_windows.items():
        window = base.loc[(base.index >= start) & (base.index <= end), "term_spread"]
        avg = window.mean()
        record(
            f"Negative term spread during {label}",
            bool(pd.notna(avg) and avg < 0.0),
            f"window mean={avg:.4%}" if pd.notna(avg) else "no values",
        )

    for month in crash_months:
        prior = base.loc[
            (base.index >= month - 12) & (base.index < month),
            ["breadth", "tail_imb"],
        ]
        current = base.reindex([month])
        breadth_now = current["breadth"].iloc[0]
        tail_now = current["tail_imb"].iloc[0]
        breadth_prior = prior["breadth"].mean()
        tail_prior = prior["tail_imb"].mean()
        record(
            f"Crash-month participation pattern {month}",
            bool(
                pd.notna(breadth_now)
                and pd.notna(tail_now)
                and pd.notna(breadth_prior)
                and pd.notna(tail_prior)
                and breadth_now < breadth_prior
                and tail_now > tail_prior
            ),
            f"breadth {breadth_now:.3f} vs prior-12m {breadth_prior:.3f}; "
            f"tail_imb {tail_now:.3f} vs prior-12m {tail_prior:.3f}",
        )

    log_dp_yield = np.exp(base["log_dp"])
    valid_dp = log_dp_yield.dropna()
    early_80s = valid_dp.loc[
        (valid_dp.index >= pd.Period("1980-01", "M"))
        & (valid_dp.index <= pd.Period("1984-12", "M"))
    ]
    max_period = valid_dp.idxmax() if len(valid_dp) else None
    median_dp = valid_dp.median() if len(valid_dp) else np.nan
    record(
        "Dividend-price yield scale and early-1980s peak",
        bool(
            pd.notna(median_dp)
            and 0.01 <= median_dp <= 0.06
            and max_period is not None
            and pd.Period("1980-01", "M") <= max_period <= pd.Period("1989-12", "M")
        ),
        f"median={median_dp:.2%}; maximum={valid_dp.max():.2%} at {max_period}; "
        f"1980-84 max={early_80s.max():.2%}" if len(valid_dp) else "no valid values",
    )

    event_months = [pd.Period("1972-12", "M"), pd.Period("1973-01", "M")]
    measures = ["breadth", "ma_part", "hi_lo", "dispersion", "tail_imb"]
    pre_event = base.loc[
        (base.index >= pd.Period("1971-12", "M"))
        & (base.index < pd.Period("1972-12", "M")),
        measures,
    ]
    event = base.reindex(event_months)[measures]
    event_changes = event.iloc[1] - event.iloc[0]
    pre_sd = pre_event.std(ddof=1)
    notable = {
        col: {
            "dec72": event.iloc[0][col],
            "jan73": event.iloc[1][col],
            "change": event_changes[col],
            "pre12m_sd": pre_sd[col],
        }
        for col in measures
    }
    flagged_measures = [
        col
        for col in measures
        if pd.notna(event_changes[col])
        and pd.notna(pre_sd[col])
        and pre_sd[col] > 0.0
        and abs(event_changes[col]) > 2.0 * pre_sd[col]
    ]
    record(
        "NASDAQ-entry discontinuity check (1972-12 to 1973-01)",
        not flagged_measures,
        f"flagged (>2 prior-12m SD): {flagged_measures or 'none'}; "
        f"changes={notable}",
    )

    for key, value in duplicate_info.items():
        messages.append(f"Stock duplicate audit: {key}={value:,}")
    for value in cpi_flags.values():
        messages.append(f"Data-quality flag: {value}")

    checks_file = CHECKS_DIR / "sanity_checks.txt"
    with checks_file.open("w", encoding="utf-8") as output:
        output.write("Predictor construction sanity checks\n")
        output.write("====================================\n")
        output.write("\n".join(messages))
        output.write("\n\nStock duplicate audit\n")
        for key, value in duplicate_info.items():
            output.write(f"{key}: {value:,}\n")
        output.write("\nData-quality flags\n")
        if cpi_flags:
            for value in cpi_flags.values():
                output.write(f"- {value}\n")
        else:
            output.write("- None.\n")

    log(f"Wrote predictor summaries, plots, and checks to {CHECKS_DIR}")
    return messages


def main() -> None:
    log("Starting predictor construction")
    stock, duplicate_info = load_stock_data()
    monthly = load_monthly_index()
    daily = load_daily_index()
    bills, fixed, cpi = load_treasury_and_cpi()

    max_month = max(monthly["ym"].max(), stock["ym"].max())
    calendar = pd.period_range(monthly["ym"].min(), max_month, freq="M")
    if calendar.min() > OUTPUT_START:
        raise ValueError(f"Input data begins after requested output start {OUTPUT_START}.")

    rates, daily_bill_h1, cpi_flags = make_rate_inflation_predictors(
        bills, fixed, cpi, calendar
    )
    log("Building market trend, risk, and dividend-price predictors")
    market = make_market_predictors(monthly, daily, daily_bill_h1, calendar)
    log("Building rates and inflation predictors")
    log("Building market participation predictors")
    participation = make_participation_predictors(stock, calendar)
    log("Building cross-sectional condition predictors")
    cross_section = make_cross_section_predictors(stock, monthly, calendar)
    log("Building NYSE/AMEX trading-activity predictors")
    activity = make_trading_activity_predictors(stock, calendar)

    base_all = make_predictor_frame(
        [market, rates, participation, cross_section, activity],
        calendar,
    )
    base_all.index.name = "ym"
    base = base_all.loc[base_all.index >= OUTPUT_START].copy()
    base = base[BASE_PREDICTORS]

    if len(base) and not base.index.equals(
        pd.period_range(OUTPUT_START, max_month, freq="M")
    ):
        raise ValueError("Base predictor calendar is not complete and monthly.")

    base_out, lagged = write_outputs(base)
    expected_lagged_columns = 1 + 6 * len(BASE_PREDICTORS)
    if len(lagged.columns) != expected_lagged_columns:
        raise ValueError(
            f"Expected {expected_lagged_columns} lagged output columns; "
            f"found {len(lagged.columns)}."
        )

    print("\nPredictors at requested event months:")
    event_months = ["1987-10", "2008-10", "2020-03"]
    event_table = base_out.loc[base_out["ym"].isin(event_months)]
    if not event_table.empty:
        print(event_table.to_string(index=False, float_format=lambda value: f"{value:.6g}"))
    else:
        print("None of the requested event months is in the output range.")

    messages = write_checks(base, duplicate_info, cpi_flags)
    flagged = [message for message in messages if message.startswith("FLAG:")]
    print(f"\nFinished. Base rows: {len(base):,}; lagged rows: {len(lagged):,}.")
    if flagged:
        print(f"Sanity checks flagged {len(flagged)} item(s); see checks/sanity_checks.txt.")
    else:
        print("All numeric sanity checks passed.")


if __name__ == "__main__":
    main()
