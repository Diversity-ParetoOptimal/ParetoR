### Pareto-Optimal Shrinkage Formula ###
# Song, Q. C., Tang, C., Newman, D. A., & Wee, S. (2023). Adverse impact reduction and job performance optimization via Pareto-optimal weighting: A shrinkage formula and regularization technique using machine learning. Journal of Applied Psychology. doi:10.1037/apl0001085.

#' Pareto_Adj function
#'
#' Estimate shrunken Pareto-optimal solution based on Pareto-optimal shrinkage formulas
#' @param Ncal Calibration sample size
#' @param wpred Matrix of predictor weights, where each row displays a Pareto-optimal solution and each column displays a weight for each predictor in that solution
#' @param Rperf_cal Vector of calibration sample job performance validity
#' @param Rrace_cal Vector of calibration sample race bivariate correlation (i.e., correlation between race dummy variable (0-minority, 1-majority) and predictor composite score)
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @return A data frame, each row shows a Pareto-optimal solution prior to and after the adjustment using the Pareto-optimal shrinkage formula
#' @return \describe{
#'   \item{Rperf_cal, Rperf_adj}{Job performance validity prior to and after adjustment, respectively}
#'   \item{Rrace_cal, Rrace_adj}{Race bivariate correlation prior to and after adjustment, respectively}
#'   \item{AIratio_cal, AIraio_adj}{Adverse impact ratio prior to and after adjustment, respectively}
#'   \item{w_perf, w_race}{The extent to which, respectively, job performance validity and race bivariate correlation is optimized in a given Pareto-optimal solution}
#'   }
#' @examples
#' # Example corresponds to Song, Tang, Newman, & Wee (2023), Supplemental Material 4
#' # Song, Q. C., Tang, C., Newman, D. A., & Wee, S. (2023). Adverse impact reduction
#' # and job performance optimization via Pareto-optimal weighting: A shrinkage formula
#' # and regularization technique using machine learning. Journal of Applied
#' # Psychology. doi:10.1037/apl0001085.
#' # (1) Calibration sample size
#' Ncal <- 40
#' # (2) Predictor weights for each Pareto solution.
#' # Rows: Pareto solutions; Columns: Predictors
#' wpred <- matrix(c(0,    0,    0,    0,    1,
#'                   0,    0,    0, 0.07, 0.93,
#'                   0,    0,    0, 0.13, 0.87,
#'                   0,    0,    0, 0.18, 0.82,
#'                   0,    0,    0, 0.23, 0.77,
#'                   0,    0,    0, 0.28, 0.72,
#'                   0,    0,    0, 0.32, 0.68,
#'                   0,    0,    0, 0.37, 0.63,
#'                   0,    0,    0, 0.41, 0.59,
#'                   0, 0.01,    0, 0.45, 0.55,
#'                   0, 0.04,    0, 0.44, 0.52,
#'                   0, 0.07,    0, 0.43,  0.5,
#'                   0,  0.1,    0, 0.42, 0.48,
#'                   0, 0.14,    0, 0.41, 0.46,
#'                   0, 0.17,    0,  0.4, 0.43,
#'                   0, 0.21,    0, 0.38, 0.41,
#'                   0, 0.25,    0, 0.37, 0.38,
#'                   0,  0.3,    0, 0.35, 0.35,
#'                   0, 0.33, 0.03, 0.34,  0.3,
#'                   0, 0.37, 0.08, 0.31, 0.25,
#'                   0, 0.42, 0.15, 0.27, 0.15), ncol = 5, nrow = 21)
#' # (3) Vector of calibration sample job performance validity
#' Rperf_cal = c(.20, .24, .27, .30, .33, .36, .39, .42, .45, .48, .51, .54,
#' .57, .60, .63, .66, .69, .72, .74, .76, .78)
#' # (4) Vector of calibration sample race bivariate correlation
#' # (i.e., correlation between race dummy variable (0-minority, 1-majority)
#' # and predictor composite score)
#' Rrace_cal = c(-.12, -.11, -.10, -.10, -.09, -.08, -.07, -.06, -.05, -.03,
#' -.02, .00, .01, .03, .05, .07, .09, .12, .15, .19, .24)
#' # (5) proportion of minority
#' prop <- 1/6
#' # (6) selection ratio
#' sr <- .15
#' # Estimate shrunken Pareto-optimal solution
#' ParetoAdj(Ncal = Ncal, wpred = wpred, Rperf_cal = Rperf_cal,
#' Rrace_cal = Rrace_cal, prop = prop, sr = sr)
#' @export

ParetoAdj <- function(Ncal, wpred, Rperf_cal, Rrace_cal, prop, sr){

  # Total number of predictors
  numpred = dim(wpred)[2]

  # Total number of Pareto solutions
  p_pareto = dim(wpred)[1]

  # Create space for output
  cnames = c("Rperf_cal","Rperf_adj","Rrace_cal","Rrace_adj", "AIratio_cal", "AIratio_adj", "w_perf","w_race")
  out = matrix(NA, p_pareto, length(cnames))
  colnames(out) = cnames

  out[, "Rperf_cal"] <- Rperf_cal
  out[, "Rrace_cal"] <- Rrace_cal
  out[, "AIratio_cal"] <- r2AI(Rrace_cal, prop, sr)

  # R2_perf where performance was maximized
  Rperf_max = Rperf_cal[p_pareto]

  # R2_perf where diversity was minimized
  Rperf_min = Rperf_cal[1]

  # R2_race where diversity was maximized
  Rrace_max = Rrace_cal[1]

  # R2_race where performance was minimized
  Rrace_min = Rrace_cal[p_pareto]

  # Number of predictors with non-zero weights on the endpoints

  p_perf0 <- sum(wpred[p_pareto, ] > 0)

  p_race0 <- sum(wpred[1, ] > 0)

  p_perf <- ifelse(p_perf0 < 3,
                   p_perf0 + (70 / (Ncal + p_perf0)),
                   p_perf0)

  p_race <- ifelse(p_race0 < 3,
                   p_race0 + (70 / (Ncal + p_race0)),
                   p_race0)

  # Obtain shrunken R^2 values of the endpoints via Claudy and Browne shrinkage formulas

  # Shrink Rperf_max
  Rperf_adj_max <- sqrt(R2_CB(N = Ncal, p = p_perf, R2 = Rperf_max^2))

  # Shrink Rrace_max
  R2race_adj_max <- R2_CB(N = Ncal, p = p_race, R2 = Rrace_max^2)
  Rrace_adj_max <- get_adj(Rrace_max, R2race_adj_max)

  # Interpolate

  for(l in 1:p_pareto){

    Rperf_p = Rperf_cal[l]
    Rrace_p = Rrace_cal[l]

    ## obtain criterion weights "w_perf", "w_race" ##

    w_perf = (Rperf_p - Rperf_max)/(Rperf_min - Rperf_max)
    w_race = (Rrace_p - Rrace_max)/(Rrace_min - Rrace_max)

    # Performance
    Rperf_adj = (w_perf*(Rperf_min - Rperf_adj_max)) + Rperf_adj_max
    # Diversity
    Rrace_adj = (w_race*(Rrace_min - Rrace_adj_max)) + Rrace_adj_max

    out[l, "Rperf_adj"] <- Rperf_adj
    out[l, "Rrace_adj"] <- Rrace_adj
    out[l, "AIratio_adj"] <- r2AI(Rrace_adj, prop, sr)
    out[l, "w_perf"] <- w_perf
    out[l, "w_race"] <- w_race
  }

  return(out)

}

##### Claudy & Browne Shrinkage Formula ######

#' R2_CB function
#'
#' Estimate shrunken R2 based on Browne (1975) formula
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @return R2_B formula-adjusted R2 based on Browne (1958) shrinkage formula
#' @examples
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_CB(N = N, p = p, R2 = R2)
#' @import hypergeo
#' @export
R2_CB <- function(N, p, R2){
  R2_C <- 1 - ((N - 4) / (N - p - 1) * (1 - R2)) * (1 + 2 * (1 - R2)/(N - p + 1))
  R2_C[R2_C < 0] <- 0
  R2_B <- ((N-p-3)*get_rho4(N, p, R2_C) + R2_C)/((N - 2*p - 2) * R2_C + p)

  return(R2_B)
}

## r2AI ##
# Convert r to AIratio
r2AI <- function(rrace,prop,sr){
  r <- rrace #r_race; point-biserial correlation
  p <- prop #proportion of Black applicants
  z_cut <- qnorm(1-sr) #standard normal cut score across all applicants
  #(which is a direct transformation of the overall selection ratio, SR)

  # d: subgroup d (Black-White)
  d <- r/sqrt((1-r^2)*p*(1-p))
  x_cut <- z_cut*sqrt(1+d^2*(p*(1-p))) - d*p
  trans.AIratio <- (1.64*x_cut+sqrt(0.76*(x_cut)^2+4))/(sqrt((exp(d^2+2*x_cut*d)))*(1.64*(x_cut+d)+sqrt(0.76*(x_cut+d)^2+4)))
  return(trans.AIratio)
}

get_rho4 <- function(N, p, rho2) {
  rho4 <- rho2^2 - (2*p*(1 - rho2)^2 / ((N - 1) * (N - p + 1)))
  rho4[rho4 < 0] <- 0
  return(rho4)
}

get_adj <- function(R, R2_adj) {
  return(R - (sqrt(R2_adj) - abs(R)))
}


