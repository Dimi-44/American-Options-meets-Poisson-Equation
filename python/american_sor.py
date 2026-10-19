import numpy as np
import QuantLib as ql
def price_american_option_sor(
    S0=100.0, K=100.0, T=1.0, r=0.05, sigma=0.2,
    option_type="put", S_max=300.0, M=100, N=1000,
    omega=1.2, tol=1e-6, align_strike=False
):
    """Prices an American option using FDM Implicit Scheme + SOR."""
    if align_strike:
        k_node = int(round((K / S_max) * M))
        dS = K / k_node
        S_max = dS * M
    else:
        dS = S_max / M
    dt = T / N
    S = np.linspace(0, S_max, M + 1)
    if option_type.lower() == "put":
        payoff = np.maximum(K - S, 0.0)
    else:
        payoff = np.maximum(S - K, 0.0)
    V = payoff.copy()
    i = np.arange(1, M)
    if align_strike:
        omega = 2.0 / (1.0 + np.sin(np.pi / M))
    # PDE Coefficients
    a = -0.5 * dt * (sigma**2 * i**2 - r * i)
    b = 1.0 + dt * (sigma**2 * i**2 + r)
    c = -0.5 * dt * (sigma**2 * i**2 + r * i)
    for n in range(1, N + 1):
        V_old = V.copy()
        tau = n * dt
        if option_type.lower() == "put":
            V[0] = K * np.exp(-r * tau)
            V[M] = 0.0
        else:
            V[0] = 0.0
            V[M] = S_max - K * np.exp(-r * tau)
        for _ in range(10000):
            max_diff = 0.0
            for k in range(1, M):
                idx = k - 1
                gs_val = (V_old[k] - a[idx] * V[k - 1] - c[idx] * V[k + 1]) / b[idx]
                V_new = (1.0 - omega) * V[k] + omega * gs_val
                V_new = max(V_new, payoff[k])
                diff = abs(V_new - V[k])
                if diff > max_diff:
                    max_diff = diff
                V[k] = V_new
            if max_diff < tol:
                break
    return float(np.interp(S0, S, V))
def price_quantlib_american(S0=100.0, K=100.0, T=1.0, r=0.05, sigma=0.2, option_type="put"):
    """Reference Price via QuantLib Engine."""
    today = ql.Date.todaysDate()
    ql.Settings.instance().evaluationDate = today
    
    maturity = today + ql.Period(int(T * 365), ql.Days)
    day_count = ql.Actual365Fixed()
    
    ql_type = ql.Option.Put if option_type.lower() == "put" else ql.Option.Call
    payoff = ql.PlainVanillaPayoff(ql_type, K)
    exercise = ql.AmericanExercise(today, maturity)
    option = ql.VanillaOption(payoff, exercise)
    
    spot_handle = ql.QuoteHandle(ql.SimpleQuote(S0))
    flat_ts = ql.YieldTermStructureHandle(ql.FlatForward(today, r, day_count))
    flat_vol = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.NullCalendar(), sigma, day_count))
    
    process = ql.BlackScholesProcess(spot_handle, flat_ts, flat_vol)
    engine = ql.FdBlackScholesVanillaEngine(process, 200, 200)
    option.setPricingEngine(engine)
    
    return option.NPV()
if __name__ == "__main__":
    params = dict(S0=100.0, K=100.0, T=1.0, r=0.05, sigma=0.2)
    print("\n" + "="*80)
    print(" PYTHON AMERICAN OPTION BENCHMARK: BASE SOR vs. OPTIMIZED SOR vs. QUANTLIB")
    print("="*80)
    for opt_type in ["put", "call"]:
        base_price = price_american_option_sor(**params, option_type=opt_type, M=100, N=1000, omega=1.2)
        opt_price  = price_american_option_sor(**params, option_type=opt_type, M=300, N=3000, align_strike=True)
        ql_price   = price_quantlib_american(**params, option_type=opt_type)
        print(f"\n[{opt_type.upper()} OPTION]")
        print(f"Base SOR Price (M=100, N=1000):  {base_price:.4f}  | Diff vs QL: {abs(base_price - ql_price):.6f}")
        print(f"Opt SOR Price (M=300, N=3000):   {opt_price:.4f}  | Diff vs QL: {abs(opt_price - ql_price):.6f}")
        print(f"QuantLib Reference Price:       {ql_price:.4f}")
