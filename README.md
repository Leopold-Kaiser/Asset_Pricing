# Demand Aggregation in Asset Pricing

This project proposes a portfolio construction framework that estimates **asset demand by aggregating multiple portfolio strategies**. The idea is that the demand for each asset can be inferred from the allocation decisions of different portfolio optimization methods.

By combining the weights assigned by several portfolio strategies, the model estimates the implicit demand for each asset and derives optimal weights for the strategies themselves. The resulting portfolio, called the **Demand Matrix Portfolio (DMP)**, dynamically allocates capital based on estimated demand signals.

The approach is evaluated through **backtesting on historical financial data** and compared to benchmark strategies such as an equally weighted portfolio and the S&P 500 index.



# Repository Structure

```
Demand_Aggregation_Asset_Pricing.R
---
Main R script implementing the portfolio strategies, demand aggregation
methodology, and backtesting framework.

Demand_Aggregation_Asset_Pricing.pdf
---
Project report explaining the theoretical motivation, mathematical model,
estimation approach, and backtesting results.
```



# Model Idea

The core assumption of the model is that **asset demand can be approximated by the relative returns of assets within a portfolio**.

Let

- \(R_i\) = return of asset \(i\)
- \(D_i\) = demand for asset \(i\)
- \(w_{ij}\) = weight assigned to asset \(i\) by portfolio strategy \(j\)
- \(p_j\) = weight of portfolio strategy \(j\)

The model assumes:

```
R_i / Σ R_j = D_i = Σ w_ij p_j
```

This relationship implies that **asset demand results from the weighted combination of portfolio strategies**.

The goal is to estimate the vector \(p_j\), which represents the importance of each portfolio strategy in determining market demand.



# Estimation Approach

Three potential approaches are considered for estimating the portfolio strategy weights:

### Weight Combination Search

A brute-force method that evaluates all possible combinations of portfolio weights.  
This approach quickly becomes computationally infeasible due to the exponential growth of possible combinations.

### Analytical Solution

An analytical solution would require a square system of equations (equal number of assets and strategies).  
Since the number of assets exceeds the number of strategies, the system is **overidentified**, making an exact solution impossible.

### Regression Approach

The final implementation uses a **least squares regression framework** to estimate the strategy weights.

The regression model takes the form:

```
D_i = w_i^T β + u_i
```

where

- \(β\) represents the contribution of each portfolio strategy
- \(u_i\) represents unexplained demand components.

Negative coefficients are removed to avoid short selling, and the remaining coefficients are **normalized to obtain valid portfolio weights**.



# Portfolio Strategies

The model aggregates several well-known portfolio construction methods:

- Equal Weight Portfolio (EWP)
- Markowitz Mean–Variance Portfolio
- Global Minimum Variance Portfolio (GMVP)
- Maximum Sharpe Ratio Portfolio
- Inverse Volatility Portfolio
- Risk Parity Portfolio
- Most Diversified Portfolio
- Maximum Decorrelation Portfolio

Each strategy generates a set of asset weights.  
These weights form the **demand matrix**, which is used to estimate the overall demand structure of the market.



# Demand Matrix Portfolio (DMP)

The **Demand Matrix Portfolio** combines the weights from the different portfolio strategies using the estimated global strategy weights.

The resulting portfolio reflects the **aggregated demand implied by multiple portfolio optimization frameworks**.

Weights are recalibrated periodically using updated market data.



# Backtesting Framework

The strategy is evaluated using historical financial data with the following setup:

- 22 assets across multiple sectors and asset classes
- Random subsamples of 10 assets per dataset
- 100 datasets generated through resampling
- Two-year investment horizons
- Rebalancing every **63 trading days**

Performance is compared against:

- **Equal Global Weight Portfolio (EGWP)**
- **S&P 500 Index**



# Results

Backtesting results show that the **Demand Matrix Portfolio (DMP)** delivers strong performance in terms of return and risk-adjusted metrics.

Key performance indicators:

| Metric | DMP | EGWP | S&P 500 |
|------|------|------|------|
| Sharpe Ratio | 1.959 | 0.743 | 0.698 |
| Sortino Ratio | 6.345 | 1.051 | 0.965 |
| Annual Return | 6.919% | 0.103% | 0.116% |

While the strategy achieves strong returns and risk-adjusted performance, it also exhibits **higher volatility and drawdowns**, suggesting that the strategy may capture speculative or demand-driven dynamics.



# Interpretation

The Demand Matrix Portfolio can be interpreted in two ways:

1. **Demand-driven asset allocation**  
   The strategy identifies assets experiencing increasing demand from common portfolio strategies.

2. **Fundamental value aggregation**  
   Since many portfolio strategies implicitly incorporate economic signals (risk, diversification, expected return), the aggregated demand may reflect a combination of different definitions of fundamental value.



# Requirements

The implementation uses the following R packages:

```
portfolioBacktest
quantmod
PerformanceAnalytics
CVXR
quadprog
riskParityPortfolio
pacman
DT
```

These libraries provide functionality for financial data acquisition, portfolio optimization, convex optimization, and performance evaluation.



# Reproducibility

All scripts required to reproduce the analysis and backtesting results are included in this repository.

Running the R script will:

1. Download historical financial data
2. Generate resampled datasets
3. Construct multiple portfolio strategies
4. Estimate strategy weights via regression
5. Build the Demand Matrix Portfolio
6. Perform backtesting and compute performance metrics



# Academic Context

This project was completed as part of the **Asset Pricing & Portfolio Theory** course in the Postgraduate Program in Data Science for Finance.

The work explores how portfolio strategies can be combined to infer market demand and improve portfolio construction.



# License

This repository is provided for academic and educational purposes.
