# Command Function (ParetoElnet) - Pareto-Optimization via with Elastic net-like Regularization
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com
# Song, Q. C., Tang, C., Newman, D. A., & Wee, S. (2023). Adverse impact reduction and job performance optimization via Pareto-optimal weighting: A shrinkage formula and regularization technique using machine learning. Journal of Applied Psychology. doi:10.1037/apl0001085.

#' ParetoElnet
#'
#' Command function to run the regularized tradeoff curve algorithm.
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @param d Subgroup difference
#' @param R Correlation matrix
#' @param Spac Number of Pareto points; default is 10
#' @param lambda.param lambda parameter, adjusts the magnitude of regularization within the elastic net function; default value is 0
#' @param alpha.param alpha parameter, adjusts the extent of ridge and LASSO within the elastic net function; default value is 1
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @param display_solution If TRUE, Pareto-optimal solutions will be printed
#' @return A list containing the criterion values (Pareto_Fmat) and predictor weights (Pareto_Xmat) of Pareto-optimal solutions
#' @examples
#' # Specify inputs
#' # (1) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (2) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (3) Subgroup differences (d): standardized mean differences between minority
#' # and majority subgroups (i.e., majority - minority), on each predictor (in applicant pool)
#' d <- c(1.00, 0.23, 0.09, 0.33)
#' # (4) Correlation matrix (R) = criterion & predictor inter-correlation matrix (in applicant pool)
#' # Format: Predictor_1, ..., Predictor_n, Criterion
#' R <- matrix(c(1, .24, .00, .19, .30,
#'               .24, 1, .12, .16, .30,
#'               .00, .12, 1, .51, .18,
#'               .19, .16, .51, 1, .28,
#'               .30, .30, .18, .28, 1),
#'             (length(d)+1),(length(d)+1))
#' # Fit Regularized Pareto-optimal model
#' out = ParetoElnet(prop, sr, d, R)
#'
#' @export
ParetoElnet <- function(prop, sr, d, R, Spac = 10, lambda.param = 0, alpha.param = 1, graph = TRUE, display_solution = TRUE){

  # Hyperparameters
  a <<- alpha.param # alpha, adjusts the extend of ridge and LASSO within the elastic net function
  lambda <<- lambda.param # lambda, adjusts the magnitude of regularization within the elastic net function

  # Input variables
  prop_ParetoR <<- prop
  sr_ParetoR <<- sr
  d_ParetoR <<- d
  R_ParetoR <<- R

  # Tolerance Level for Algorithm
  TolCon 	= 1.0000e-6 # tolerance of constraint
  TolF 	= 1.0000e-6 # tolerance of objective function
  TolX 	= 1.0000e-7 # tolerance of predictor variable

  # Do not change
  Fnum 	= 2
  Xnum = length(d)
  VLB = c(-(Xnum+1),rep(0,Xnum))
  VUB = c((Xnum+1),rep(1,Xnum))
  X0 = c(1,rep(1/Xnum,Xnum))

  ###### Find Pareto-Optimal Solution ######

  out = quiet(NBI_Elnet(X0, (Spac-1), Fnum, VLB, VUB, TolX, TolF, TolCon, graph=graph, display_solution=display_solution))
  return(out)

}


#' cvParetoElnet
#'
#' Cross-validation procedure to choose alpha and lambda values for regularized tradeoff curve algorithm
#' @param data Calibration data: a dataframe with first p columns as standardized predictor scores, the (p+1)th column as standardized job performance ratings, and the (p+2)th column (rightmost column) as a race dummy; all data should be numeric
#' @param prop Proportion of minority applicants in total applicant pool (#minority applicants/#total applicants)
#' @param sr Selection ratio
#' @param k Number of folds, default is 5
#' @param rep Number of repeats if doing repeated k-fold CV
#' @param Spac Number of Pareto-optimal solutions, default is 21
#' @param alpha.grid Vector of alpha values to try
#' @param lambda.grid Vector of lambda values to try
#' @return A list containing the following:
#' @return \describe{
#'   \item{weights}{Predictor weights obtained based on best alpha and lambda values}
#'   \item{alpha}{Best alpha parameter value}
#'   \item{lambda}{Best lambda parameter value}
#'   }
#' @examples
#' # Calibration data
#' # (1) Calibration data
#' # Format: Predictor_1, ..., Predictor_n, Job Performance Validity, Race dummy variable
#' # (e.g., 0-minority; 1-majority)
#' # All data should be numeric
#' # sample_data
#' # (3) Alpha values to try
#' alpha.grid <- seq(0, 1, length = 2)
#' # (4) Lambda values to try
#' lambda.grid <- 10^seq(1, -2, length = 4)
#' # (5) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (6) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (7) Spac = number of Pareto points
#' Spac <- 21
#' # Fit Regularized Pareto-optimal model with parameter selection via cross-validation
#' ## May take a while
#' out <- cvParetoElnet(data = sample_data, prop = prop, sr = sr, Spac = Spac,
#'                      k = 5, rep = 1,
#'                      lambda.grid = lambda.grid, alpha.grid = alpha.grid)
#' @export
cvParetoElnet <- function(data, prop, sr, k = 5, rep = 1, Spac = 21,
                          alpha.grid = 0, lambda.grid = 0) {

  data <- as.matrix(data)

  n <- dim(data)[1]
  p <- dim(data)[2] - 2

  race_dummy <- data[, dim(data)[2]]

  # parameter tuning
  best_alpha <- 0
  best_lambda <- 0
  best_score <- -1

  for(alp in alpha.grid) {
    for(lam in lambda.grid) {
      scores <- vector()
      for(i in 1:rep) {
        folds <- split(x = sample(1:n, n),
                       f = sort(rep_len(x = 1:k, length.out = n)))
        score <- vector()
        for(j in 1:k) {
          idx_test <- folds[[j]]
          data_test <- data[idx_test, ]
          data_train <- data[-idx_test, ]
          train_race <- race_dummy[-idx_test]

          R_train <- cor(data_train[, 1:(p+1)])
          d_train <- colMeans(data_train[train_race == 1, 1:p]) -
            colMeans(data_train[train_race == 0, 1:p])
          train <- quiet(ParetoElnet(prop = prop, sr = sr,
                                     d = d_train, R = R_train,
                                     alpha.param = alp, lambda.param = lam,
                                     Spac = Spac,
                                     graph = F))

          r_perf_test <- 0
          for(w in 1:dim(train$Pareto_Xmat)[1]) {
            b <- as.vector(train$Pareto_Xmat[w, ])
            yhat_test <- data_test[, 1:p] %*% b
            r_perf_test[w] <- cor(yhat_test, data_test[, p+1])
          }
          score[j] <- max(r_perf_test)
        }
        scores <- c(scores, mean(score))
      }
      current_score <- mean(scores)
      if(best_score < current_score) {
        best_score <- current_score
        best_alpha <- alp
        best_lambda <- lam
      }
    }
  }

  # refit
  R <- cor(data[, 1:(p+1)])
  d_samp <- colMeans(data[race_dummy == 1, 1:p]) - colMeans(data[race_dummy == 0, 1:p])
  mod <- quiet(
    ParetoElnet(prop = prop, sr = sr, d = d_samp, R = R, Spac = Spac,
                alpha.param = best_alpha, lambda.param = best_lambda, graph = FALSE))

  return(list(weights = mod$Pareto_Xmat,
              alpha = best_alpha,
              lambda = best_lambda))
}

# this function is from the SimDesign package (Chalmers & Adkins, 2020) to suppress messages
quiet <- function (..., messages = FALSE, cat = FALSE)
{
  if (!cat) {
    tmpf <- tempfile()
    sink(tmpf)
    on.exit({
      sink()
      file.remove(tmpf)
    })
  }
  out <- if (messages)
    eval(...)
  else suppressMessages(eval(...))
  out
}
