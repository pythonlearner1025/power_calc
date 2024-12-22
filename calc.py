import math

def z_value_for_power(power: float) -> float:
    """
    Return the z-value for a given power (one - beta).
    Examples:
      - 70% power -> z ~ 0.524
      - 80% power -> z ~ 0.84
      - 90% power -> z ~ 1.28
    """
    # Common approximate z-values for typical powers:
    power_to_z = {
        0.70: 0.524,  # ~ 70%
        0.80: 0.84,   # ~ 80%
        0.90: 1.28    # ~ 90%
    }
    return power_to_z.get(power, None)

def required_sample_size(delta: float, sigma_p: float, power: float, alpha: float = 0.05) -> int:
    """
    Compute the sample size needed per group (two-sided test) for
    given delta, pooled SD (sigma_p), desired power, and alpha=0.05.
    
    Formula:
        n = 2 * ( (Z_(1-alpha/2) + Z_(1-beta)) / (delta / sigma_p) )^2
    """
    # Two-sided alpha => z for alpha/2
    z_alpha = 1.96  # for alpha = 0.05
    z_power = z_value_for_power(power)
    
    # ratio = (sigma_p / delta), using absolute difference in the denominator
    ratio = sigma_p / abs(delta)
    
    # (Z_(1-alpha/2) + Z_(1-beta))^2
    sum_z = z_alpha + z_power
    sum_z_sq = sum_z ** 2
    
    # n per group
    n = 2.0 * sum_z_sq * (ratio ** 2)
    
    # Round up to the next integer
    return math.ceil(n)

def compute_pooled_sd(delta: float, t_stat: float, n_each: int) -> float:
    """
    Back-calculate pooled SD from:
        t = delta / sqrt(2 * sigma_p^2 / n_each)
        => sigma_p = ( delta * sqrt(n_each) ) / (sqrt(2) * t_stat)
    """
    return (abs(delta) * math.sqrt(n_each)) / (math.sqrt(2) * t_stat)

# -----------------------------------------------------------------------------
# Define each metric’s parameters in a list of dicts
# -----------------------------------------------------------------------------

from scipy.stats import t as tdist

def t_exact(p,df):
    return tdist.ppf(1 - p/2, df)

metrics = [
    {
        "name": "Total mobile", 
        "delta": 90.0,     # (mutant - normal)
        "p_value": 0.0001, # reported
        "n_normal": 19, 
        "n_mutant": 19,
        "df": 36,          # degrees of freedom = n1 + n2 - 2
        "t_approx": t_exact(0.0001,36)    # approximate t from p=0.0001, df=36
    },
    {
        "name": "Total sleep duration",
        "delta": -60.0,   # negative difference
        "p_value": 0.05,
        "n_normal": 9,
        "n_mutant": 9,
        "df": 16,
        "t_approx": t_exact(0.05,16)  # approximate t from p=0.05, df=16
    },
    {
        "name": "REM duration",
        "delta": -10.0,
        "p_value": 0.05,
        "n_normal": 9,
        "n_mutant": 9,
        "df": 16,
        "t_approx": t_exact(0.05,16)
    },
    {
        "name": "NREM duration",
        "delta": -100.0,
        "p_value": 0.05,
        "n_normal": 9,
        "n_mutant": 9,
        "df": 16,
        "t_approx": t_exact(0.05, 16)
    },
    {
        "name": "Delta Power (NREM)",
        "delta": 1.0,     # ±1 difference
        "p_value": 0.001,
        "n_normal": 9,
        "n_mutant": 9,
        "df": 16,
        "t_approx": t_exact(0.001, 16)
    }
]

# -----------------------------------------------------------------------------
# Main calculation & printing
# -----------------------------------------------------------------------------

powers = [0.70, 0.80, 0.90]  # 70%, 80%, 90% power

# Table header
header = (
    f"{'Metric':<25} "
    f"{'Delta':>7} "
    f"{'p':>7} "
    f"{'n_normal':>9} "
    f"{'n_mutant':>9} "
    f"{'t_exact':>9} "
    f"{'PooledSD':>10} "
    + "  " + "  ".join([f"n@{int(p*100)}%" for p in powers])
)
print(header)
print("-" * len(header))

# For each metric, compute pooled SD & required sample sizes
for m in metrics:
    name      = m["name"]
    delta     = m["delta"]
    p_val     = m["p_value"]
    n_normal  = m["n_normal"]
    n_mutant  = m["n_mutant"]
    df        = m["df"]
    t_approx  = m["t_approx"]
    
    # Compute pooled SD
    # n_each = n_normal = n_mutant (given equal sizes in this scenario)
    n_each = n_normal
    sigma_p = compute_pooled_sd(delta=delta, t_stat=t_approx, n_each=n_each)
    
    # Compute sample sizes for each desired power
    # (same formula, just different Z_(1-beta))
    sample_sizes = [required_sample_size(delta, sigma_p, pwr) for pwr in powers]
    
    # Format row
    row = (
        f"{name:<25} "
        f"{delta:>7.1f} "
        f"{p_val:>7g} "
        f"{n_normal:>9} "
        f"{n_mutant:>9} "
        f"{t_approx:>9.2f} "
        f"{sigma_p:>10.2f}   "
        + "  ".join([f"{nreq:>3}" for nreq in sample_sizes])
    )
    print(row)

