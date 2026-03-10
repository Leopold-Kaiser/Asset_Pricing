#install.packages("pacman")
# install & load packages
library(pacman)
p_load(portfolioBacktest)     # portfolio Backtesting
p_load(quantmod)              # to download stock data
p_load(PerformanceAnalytics)  # to compute performance measures
p_load(CVXR)                  # Convex Optimization in R
p_load(DT)                    # R interface to the JavaScript library DataTables  # to compute performance measures
p_load(quadprog)
p_load(riskParityPortfolio)

set.seed(321)


startDate <- "2010-01-01"
endDate <- Sys.Date()

# Download S&P 500 Index data separately
getSymbols("^GSPC", from = startDate, to = endDate)
sp500_index_data <- GSPC

# Corrected version of the stockDataDownload call
BackT_daily_data <- stockDataDownload(c("HPQ",        # HP Inc.
                  "MSTR",       # MicroStrategy Inc.
                  "AMD",        # Advanced Micro Devices
                  "GOOGL",      # Alphabet Inc. (Google Class A)
                  "MSFT",       # Microsoft Corporation
                  "AAPL",       # Apple Inc.
                  "WMT",         
                  "KO",          
                  "PLD",        # Prologis, Inc.
                  "JNJ",        # Johnson & Johnson
                  "V",          # Visa Inc.
                  "MA",         # Mastercard Inc.
                  "NSRGY",      # Nestle Nigeria Plc (Nigerian Stock Exchange)
                  "TAN",        # Invesco Solar ETF
                  "GLD",        # Gold ETF (for Gold price in USD)
                  "VNQ",        # Vanguard Real Estate Index Fund ETF
                  "DBA",        #Invesco DB Agriculture Fund
                  "EEM",        # iShares MSCI Emerging Markets ETF
                  "BRK-A",
                  "JPM",
                  "STLD",
                  "HBAN"),
  index_symbol = "^GSPC",   # Include S&P 500 Index for benchmarking
  from = startDate,
  to = endDate,
  rm_stocks_with_na = TRUE,
  #local_file_path = getwd()
)



# Print the result to verify
BackT_daily_data

my_dataset_list <- financialDataResample(BackT_daily_data,
                                         N_sample = 10,     # Desired number of financial instruments in each resample
                                         T_sample = 252*2,  # Desired length of each resample (consecutive samples with a random initial time)
                                         num_datasets = 100) # Number of resampled datasets

head(my_dataset_list$`dataset 5`$adjusted,3)
tail(my_dataset_list$`dataset 5`$adjusted,3)

equal_weight_portfolio_fun <- function(dataset, ...) {
  num_assets <- ncol(dataset$adjusted)
  weights <- rep(1 / num_assets, num_assets)
  names(weights) <- colnames(dataset$adjusted)
  return(weights)
}

# Markowitz mean-variance portfolio

Markowitz_portfolio_fun <- function(dataset, lambda=0.25, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu    <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)       # compute the SCM
  # design mean-variance portfolio
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lambda*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}



#GMVP (with heuristic not to allow short-selling)

GMVP_portfolio_fun <- function(dataset, ...) {
  #that is a function of a data set, that has to be called for a specific one (from 20) for calculating the quintile weights
  
  # Calculate log-returns
  returns <- diff(log(dataset$adjusted))[-1, ] # Remove the first NA row
  
  Sigma <- cov(returns)
  # compute SCM
  # computes the covariance matrix --> N * N matrix (in our case 24 * 24)
  
  # design GMVP
  
  w <- solve(Sigma, rep(1, nrow(Sigma)))
  #is used to compute a portfolio's weights that minimize risk (variance)
  #subject to certain constraints, such as being fully invested (the sum of the weights equals 1).
  
  
  #rep(1, nrow(Sigma)): This creates a vector of 1's with the same length as the number of rows (or columns) in Sigma (since Sigma is a square matrix).
  
  #The purpose of this is to indicate that the sum of the portfolio weights should equal 1,
  #implying that the portfolio is fully invested.
  
  # solve(Sigma, ...): solve(A, b) solves a system of linear equations A * x = b, where A is a square matrix and b a vector
  # consequently b in our example w (the weights) has to be a vector
  
  w <- abs(w)/sum(abs(w))
  # normalizes the vector w to ensure that the sum of the absolute values of its elements equals 1,
  #while preserving the relative proportions between the elements.
  
  #abs(w): This function takes the absolute value of each element in the vector w. If any weights are negative (which could happen in certain portfolio construction scenarios, like allowing short selling), abs(w) ensures all values are positive.
  
  return(w)
}

p_load(quadprog)

max_sharpe_ratio_portfolio_fun <- function(dataset, risk_free_rate = 0,...) {
  
  # Calculate log-returns
  returns <- diff(log(dataset$adjusted))[-1]
  
  # Calculate expected returns and covariance matrix
  exp_returns <- colMeans(returns, na.rm = TRUE)
  cov_matrix <- cov(returns, use = "complete.obs")
  
  # Define the optimization problem to maximize Sharpe ratio
  Dmat <- cov_matrix                        # Quadratic component (risk)
  dvec <- rep(0, length(exp_returns))       # No linear component
  
  # Adjust constraint matrix to enforce fully invested and non-negative weights
  # Constraint matrix: fully invested, expected returns, and non-negative weights
  Amat <- cbind(1, exp_returns, diag(length(exp_returns)))  
  bvec <- c(1, risk_free_rate, rep(0, length(exp_returns)))  # Sum to 1 and non-negative constraints
  meq <- 1  # Only the sum of weights constraint is equality
  
  # Solve the quadratic program
  result <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  weights <- result$solution
  
  # Set weights that are close to zero to exactly zero
  weights[weights < 1e-10] <- 0
  names(weights) <- colnames(dataset$adjusted)
  
  return(weights)
}


# Inverse volatility portfolio

inverse_volatility_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)  # compute SCM
  w <- riskParityPortfolio(Sigma, formulation='diag')$w
  return(w)
}

# C - Risk parity portfolio

risk_parity_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)
  # design RPP
  w <- riskParityPortfolio(Sigma)$w
  return(w)
}

# Define the most diversified portfolio function
MSRP <- function(mu, Sigma) {
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w_)/sum(result$getValue(w_)))
  names(w) <- colnames(Sigma)
  return(w)
}

most_diversified_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)    
  mu = sqrt(diag(Sigma))
  w <- MSRP(mu = sqrt(diag(Sigma)), Sigma)
  return(w)
}



maximum_decorrelation_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  corr <- cor(X)  # compute SCM
  # design GMVP
  w <- solve(corr, rep(1, nrow(corr)))
  w <- abs(w)/sum(abs(w))
  return(w)
}


equal_weight_portfolio_fun(my_dataset_list$`dataset 1`)
max_sharpe_ratio_portfolio_fun(my_dataset_list$`dataset 1`)
GMVP_portfolio_fun(my_dataset_list$`dataset 1`)
Markowitz_portfolio_fun(my_dataset_list$`dataset 1`)
inverse_volatility_portfolio_fun(my_dataset_list$`dataset 1`)
risk_parity_portfolio_fun(my_dataset_list$`dataset 1`)
most_diversified_portfolio_fun(my_dataset_list$`dataset 1`)
maximum_decorrelation_portfolio_fun(my_dataset_list$`dataset 1`)



# -------------------------------------------------------------------------------------------
# Demand Matix Function
# -------------------------------------------------------------------------------------------

portfolios_our <- list(
  "EWP"= equal_weight_portfolio_fun,
  "MSRP"  = max_sharpe_ratio_portfolio_fun,
  "GMVP" = GMVP_portfolio_fun,
  "Markowitz" = Markowitz_portfolio_fun,
  "IVP"=inverse_volatility_portfolio_fun,
  "RPP"=risk_parity_portfolio_fun,
  "MDP"=most_diversified_portfolio_fun,
  "MDCP"=maximum_decorrelation_portfolio_fun
)

demand_matrix <- function(dataset, portfolio_funcs=portfolios_our, ...) {
  # Initialize an empty list to store adjusted weights for each portfolio
  dataset_weights <- list()
  
  # Loop through each portfolio function
  for (portfolio_name in names(portfolio_funcs)) {
    portfolio_func <- portfolio_funcs[[portfolio_name]]
    
    # Calculate weights using the portfolio function
    weights <- portfolio_func(dataset)
    
    # Multiply the weights by the corresponding global weight
    dataset_weights[[portfolio_name]] <- weights  # Store adjusted weights
  }
  
  return(dataset_weights)
}

demand_matrix_portfolio <- function(dataset, portfolio_funcs=portfolios_our, ...) {


#returns 
  returns <- diff(log(dataset$adjusted))[-1, ]
  returns_cumul<- Return.cumulative(returns)
  returns_rel<-returns_cumul/rowSums(returns_cumul)



# Define the function to calculate portfolio weights with global weights for a single dataset


weight_matrix <- demand_matrix(dataset, portfolios_our)



# weight_matrix is the 22 x 8 matrix of weighted portfolio weights
# and `target_returns` is a vector of length 22

# Convert the weight matrix to a data frame for use in glm()
weight_matrix_df <- as.data.frame(weight_matrix)


returns_rel_df <- as.data.frame(returns_rel)

# Make sure returns_rel_df is a vector (it's currently a data frame)
returns_rel_vector <- as.vector(t(returns_rel_df))  

# Combine the target returns with the weight matrix to prepare data for regression
data <- cbind(returns_rel = returns_rel_vector, weight_matrix_df)

# Run the glm regression (remove the intercept term with -1)
fit <- glm(returns_rel ~ . - 1, data = data)

# Get the coefficients (the global weights)
global_weights <- fit$coefficients

#let's assume that short selling is not allowed: 
# Set negative weights to zero
global_weights[global_weights < 0] <- 0

# Normalize the weights so they sum to 1
normalized_global_weights <- global_weights / sum(global_weights)


# Convert to list and assign appropriate names
global_portfolio_weights <- as.list(normalized_global_weights)
names(global_portfolio_weights) <- names(normalized_global_weights)

global_portfolios_weights <- list(
  "EWP" = global_portfolio_weights$EWP,
  "MSRP" = global_portfolio_weights$MSRP,
  "GMVP" = global_portfolio_weights$GMVP,
  "Markowitz" = global_portfolio_weights$Markowitz,
  "IVP" = global_portfolio_weights$IVP,
  "RPP" = global_portfolio_weights$RPP,
  "MDP" = global_portfolio_weights$MDP,
  "MDCP" = global_portfolio_weights$MDCP
)


 # Initialize an empty list to store adjusted weights for each portfolio
  dataset_weights <- list()
  
  # Loop through each portfolio function
  for (portfolio_name in names(portfolio_funcs)) {
    portfolio_func <- portfolio_funcs[[portfolio_name]]
    
    # Calculate weights using the portfolio function
    weights <- portfolio_func(dataset)
    
    # Multiply the weights by the corresponding global weight
    adjusted_weights <- weights * global_weights[[portfolio_name]]
    dataset_weights[[portfolio_name]] <- adjusted_weights  # Store adjusted weights
  }
  
  # Sum the adjusted weights across portfolios by asset to create demand_vector
  demand_vector <- Reduce(`+`, dataset_weights)
  
  # Return both the adjusted weights and demand vector for the dataset
  return(demand_vector)
}

demand_matrix_portfolio(BackT_daily_data, portfolios_our)




# -------------------------------------------------------------------------------------------
# Backtesting and Plotting
# -------------------------------------------------------------------------------------------



 
 
 portfolios <- list(
                   #"EWP"= equal_weight_portfolio_fun,
                   #"MSRP"  = max_sharpe_ratio_portfolio_fun,
                   #"GMVP" = GMVP_portfolio_fun,
                   #"Markowitz" = Markowitz_portfolio_fun,
                   #"IVP"=inverse_volatility_portfolio_fun,
                   #"RPP"=risk_parity_portfolio_fun,
                   #"MDP"=most_diversified_portfolio_fun,
                   #"MDCP"=maximum_decorrelation_portfolio_fun)
                   "DMP"=demand_matrix_portfolio
                   )



bt<- portfolioBacktest(portfolios, my_dataset_list, benchmark = c("1/N","index"),
                       rebalance_every = 63,
                       optimize_every = 63, #The more frequent u rebalance
                       lookback = 252*0.25,
                       shortselling=F,
                       cost =list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4))
# Test strategis accounting and excluding trading costs
bt$DMP$`dataset 1`$cpu_time*60
names(bt)

# Portfolio's performance in tabless

# Summary of performance measures:
res_sum <- backtestSummary(bt, summary_fun = median, show_benchmark = TRUE)
names(res_sum)

round(res_sum$performance_summary,3)


