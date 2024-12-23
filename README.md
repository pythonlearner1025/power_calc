# Sample Size Computation

This repository provides a Python script to compute required sample sizes for specific metrics (e.g., total mobility, sleep duration, etc.) at 70%, 80%, and 90% power.  

## Steps

1. **Floor Sample Sizes**: For each endpoint, we set the current \(n\) for WT (wild type) and mutant to the minimum of the two reported sample sizes.  
2. **Approximate \(t\)-Value**: Use the two-sided \(p\)-value and degrees of freedom (\(df = n_{\text{normal}} + n_{\text{mutant}} - 2\)) to approximate the observed \(t\)-statistic via `scipy.stats.t.ppf()`.  
3. **Compute Pooled SD**: Derive the pooled standard deviation (\(\sigma\)) using  
   \[
     \sigma = \frac{\delta}{t} \sqrt{\frac{n}{2}} 
   \]
   where \(\delta\) is the mean difference, and \(n\) is the (floored) current sample size in each group.  
4. **Compute Required \(n\)**: Use a \(z\)-based formula for power to find new sample sizes per group for 70%, 80%, and 90% power:  
   \[
     n_{\text{per group}} 
       = 2 \left(\frac{z_{\alpha/2} + z_{\text{power}}}{\delta / \sigma}\right)^2.
   \]
5. **Output**: Print a final table containing each metric’s: mean difference, \(p\)-value, current sample sizes, exact \(t\), pooled SD, and required sample sizes for 70%, 80%, and 90% power.

## Assumptions

- Sample sizes for each endpoint are floored to \(\min(n_{\text{WT}}, n_{\text{mutant}})\).  
- Pooled SD is used for the standard deviation.

## Output

| Metric               | Delta  | p      | n_normal | n_mutant | t_exact | PooledSD | n@70% | n@80% | n@90% |
|----------------------|--------|--------|----------|----------|---------|----------|-------|-------|-------|
| Total mobile         | 90.0   | 0.0001 | 19       | 19       | 4.37    | 63.43    | 7     | 8     | 11    |
| Total sleep duration | -60.0  | 0.05   | 9        | 9        | 2.12    | 60.04    | 13    | 16    | 22    |
| REM duration         | -10.0  | 0.05   | 9        | 9        | 2.12    | 10.01    | 13    | 16    | 22    |
| NREM duration        | -100.0 | 0.05   | 9        | 9        | 2.12    | 100.07   | 13    | 16    | 22    |
| Delta Power (NREM)   | 1.0    | 0.001  | 9        | 9        | 4.01    | 0.53     | 4     | 5     | 6     |

