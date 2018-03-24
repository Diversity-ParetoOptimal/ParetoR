# Command Function (ParetoElnet) - Pareto-Optimization via with Elasticnet-like Regularization
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com
# Last Update: 03/2018

#' ParetoElnet
#'
#' Command function to run Regularized Pareto-Optimal algorithm.
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @param d Subgroup difference
#' @param R Correlation matrix
#' @param Spac Number of Pareto points (e.g., 21)
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @param display_solution If TRUE, Pareto-optimal solutions will be printed
#' @return out Pareto-Optimal solution with criterion values (Criterion) and predictor weights (ParetoWeights)
#' @examples
#' # Specify inputs
#' # (1) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (2) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (3) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor (in applicant pool)
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
ParetoElnet = function(prop, sr, d, R, Spac = 21, graph = TRUE, display_solution = TRUE){

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

  out = NBI_Elnet(X0, (Spac-1), Fnum, VLB, VUB, TolX, TolF, TolCon, graph=graph, display_solution=display_solution)
  return(out)

}


#' cv.ParetoElnet
#'
#' Cross-validation procedure to choose alpha and lambda values for regularized Pareto-optimal algorithm
#' @param data_cal calibration data
#' @param n_cal Calibration data sample size
#' @param D Subgroup difference: standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor and criterion 
#' @param R Correlation matrix: Criterion predictor inter-correlation matrix (in applicant pool); format: Predictor_1, ..., Predictor_n, Criterion
#' @param n_val validation data sample size
#' @param alpha.grid grid of alpha values to try
#' @param lambda.grid grid of lambda values to try
#' @param prop proportion of minority applicants in total applicant pool (#minority applicants/#total applicants)
#' @param sr selection ratio
#' @param Spac number of Pareto-optimal solutions
#' @import lm.beta
#' @return alpha.best best alpha parameter value
#' @return lambda.best best lambda parameter value
#' @return out_cal.best calibration model output based on best alpha and lambda values
#' @note Users could choose to use calibration sample data as input or calibration sample characteristic [i.e., calibration sample size (n_cal), subgroup difference in predictor and criterion (D), criterion and predictor correlation matrix (R)]
#' @examples
#' #' ### Example 1 ###
#' # Calibration data as input
#' # (1) Calibration data
#' # Format: Predictor_1, ..., Predictor_n, Job Performance Validity, Race dummy variable (e.g., 0-minority; 1-majority)
#' # Example
#' load(data_cal)
#' # (2) Validation sample size
#' n_val = 10000
#' # (3) Grid of alpha values to try
#' # Example
#' alpha.grid <- seq(0,1,length=3)
#' # (4) Grid of lambda values to try
#' # Example
#' lambda.grid <- 10^seq(1,-2,length=11)
#' # (5) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (6) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (7) Spac = number of Pareto points
#' Spac <- 21
#' # Fit Regularized Pareto-optimal model with parameter selectoin via cross-validation
#' cv.out = cv.ParetoElnet(data_cal = data_cal, n_val = 10000, 
#'                         lambda.grid = lambda.grid, alpha.grid = alpha.grid,
#'                         prop = prop, sr = sr, Spac = Spac)
#' ### Example 2 ###
#' # n_cal, D, R as input (input do not include calibration sample data set)
#' # (1) Calibration sample size
#' n_cal = 100
#' # (2) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor and criterion (in applicant pool)
#' d <- c(1.00, 0.23, 0.09, 0.33, .30)
#' # (3) Correlation matrix (R) = Criterion predictor inter-correlation matrix 
#' # Format: Predictor_1, ..., Predictor_n, Criterion
#' R <- matrix(c(1, .24, .00, .19, .30,
#'               .24, 1, .12, .16, .30,
#'               .00, .12, 1, .51, .18,
#'               .19, .16, .51, 1, .28,
#'               .30, .30, .18, .28, 1),
#'             (length(d)+1),(length(d)+1))
#' # (4) Validation sample size
#' n_val = 10000
#' # (5) Grid of alpha values to try
#' # Example
#' alpha.grid <- seq(0,1,length=3)
#' # (6) Grid of lambda values to try
#' # Example
#' lambda.grid <- 10^seq(1,-2,length=11)
#' # (7) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (8) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (9) Spac = number of Pareto points
#' Spac <- 21
#' # Fit Regularized Pareto-optimal model with parameter selection via cross-validation
#' cv.out = cv.ParetoElnet(n_cal = n_cal, D = D, R = R,  
#'                         n_val = n_val, 
#'                         lambda.grid = lambda.grid, alpha.grid = alpha.grid,
#'                         prop = prop, sr = sr, Spac = Spac)
#'
#' @export
cv.ParetoElnet = function(data_cal = NULL, 
                          n_cal = NULL, D = NULL, R = NULL, 
                          n_val = 10000, 
                          lambda.grid = 10^seq(1, -2, length = 21),
                          alpha.grid = seq(0, 1, length = 6),
                          # lambda2.grid = seq(30, 50, length = 2),
                          lambda2.grid = 50,
                          prop, sr,  Spac = 20)
{ require(lm.beta)

  if(!(is.null(data_cal) | (is.null(R) & is.null(D) & is.null(n_cal)))){
    stop('Need data_cal value or R and D values.')
  }else{
    if(is.null(data_cal)==FALSE){
      d_cal = (colMeans(subset(data_cal, Race_dummy == 1, select = -c(Race_ind, Race_dummy)), na.rm=T)
               - colMeans(subset(data_cal, Race_dummy == 0, select = -c(Race_ind, Race_dummy)), na.rm=T))
      R_cal = cor(subset(data_cal, select = -c(Race_ind, Race_dummy)), use = "pairwise.complete.obs") 
    }else{
      d_cal = D
      R_cal = R
      # Generate calibration sample
      data_cal = gen_data(prop = prop, sr = sr, n = n_cal, 
                          D = d_cal, R = R_cal)
    }
  
  }
  
  for(alp in alpha.grid){
    
    for(lam in lambda.grid){
      
      for(lam2 in lambda2.grid){
        
        # cat(paste0("Trying alpha = ", round(alp,2), "; lambda = ", round(lam,2), "; lambda2 = ", round(lam2,2)))
        cat(paste0("Trying alpha = ", round(alp,2), "; lambda = ", round(lam,2)))
        
        a <<- alp
        lambda <<- lam
        lambda2 <<- lam2
        
        prop.minority = prop
        prop.majority = 1 - prop
        
        # Generate validation sample
        data_val = gen_data(prop = prop, sr = sr, n = n_val, 
                            D = d_cal, R = R_cal)
        
        # Obtain calibration sample ParetoElnet solution
        model_cal = ParetoElnet(prop = prop, sr = sr, 
                             d = d_cal[-length(d_cal)], R = R_cal,
                             Spac = Spac, graph = FALSE, display_solution = FALSE)
        
        # Obtain calibration sample predictor weights
        w_cal = model_cal$Pareto_Xmat
        
        # Calibration sample
        out_cal = calculate_perf_rrace(data_cal, w_cal, prop, sr)
        # Validation sample
        out_val = calculate_perf_rrace(data_val, w_cal, prop, sr)
        
        # Regression validation Perf and rrace results
        model_reg = lm(Perf~., data = subset(data_cal, select = -c(Race_ind, Race_dummy)),
                       na.action = na.omit)
        w_reg = coef(lm.beta(model_reg))[-1]
        
        out_val_reg = calculate_perf_rrace(data_val, t(w_reg), prop, sr)
        
        RMSE = sqrt(mean((out_cal$Perf-out_val$Perf)^2))+sqrt(mean((out_cal$rrace-out_val$rrace)^2))

        diff = (max(out_val$Perf)-out_val_reg$Perf)
        
        # Threshold for RMSE-diff comparison
        comp = 1000
        
        if(comp > (RMSE-diff)){
          alpha.best = alp
          lambda.best = lam
          # lambda2.best = lam2
          model_cal.best = model_cal
          comp = (RMSE-diff)
        }
        
      }
      
    }
    
  }
  
  if(comp == 1000){
    stop("Best alpha and lambda values not found, try other grid ranges.")
  }else{
    cv.out = list("alpha.best" = alpha.best, "lambda.best" = lambda.best,
                  # "lambda2.best" = lambda2.best, 
                  "model_cal.best" = model_cal.best)
  }
 
 return(cv.out)
  
}

# Additional function -- generate validation sample

#' gen_data
#'
#' Function to generate sample (used in cv.ParetoElnet)
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @param D Subgroup difference: standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor and criterion 
#' @param R Correlation matrix: Criterion predictor inter-correlation matrix (in applicant pool); format: Predictor_1, ..., Predictor_n, Criterion
#' @return data_out Generated sample data
#' @examples
#' # Specify inputs
#' # (1) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (2) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (3) Sample size
#' n <- 100
#' # (4) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor 
#' D <- c(1.00, 0.23, 0.09, 0.33, .30)
#' # (5) Correlation matrix (R) = criterion & predictor inter-correlation matrix (in applicant pool)
#' # Format: Predictor_1, ..., Predictor_n, Criterion
#' R <- matrix(c(1, .24, .00, .19, .30,
#'               .24, 1, .12, .16, .30,
#'               .00, .12, 1, .51, .18,
#'               .19, .16, .51, 1, .28,
#'               .30, .30, .18, .28, 1),
#'             (length(d)+1),(length(d)+1))
#' # Generate data
#' data_out = gen_data(prop, sr, n, D, R)
#'
#' @export
gen_data = function(prop, sr, n, D, R){
  
  require(MASS)
  
  numvars = length(D)
  prop.minority = prop
  prop.majority = 1 - prop
  n.minority = ceiling(prop.minority*n)
  n.majority = ceiling(prop.majority*n)
  n.all = n.minority+n.majority
  
  #1) Convert d to r_pointbiserial (Newman, Jacobs, Bartram, 2007, p.1397)
  r_pb_race = matrix(0,numvars,1)
  for (i in 1:numvars){
    di = as.numeric(D[i])
    r_pb_race[i] = di/(sqrt(di^2+1/((1-prop.majority)*prop.majority)))
  }
  #2)Convert r_pointbiserial to r_biserial
  r_bis_race = matrix(0,numvars,1)
  for (i in 1:numvars){
    r_pb = r_pb_race[i]
    r_bis_race[i] = r_pb*(sqrt(prop.majority*prop.minority)/abs(dnorm(qnorm(prop.majority,0,1),0,1)))
  }
  
  r_bis_race = c(r_bis_race,1)
  
  #combine race population input with RV
  RVrace <- cbind(rbind(R, r_bis_race[1:numvars]), r_bis_race)
  
  # obtain validation sample
  sample_1 = mvrnorm(n = n.all, mu = rep(0,numvars+1), Sigma = RVrace, tol=1e-6, empirical = F)
  colnames(sample_1) = colnames(RVrace)
  #create race dummy variable
  sample_2 = sample_1[order(-sample_1[,"r_bis_race"]),]
  sample_3 = cbind(sample_2,c(rep(1,n.majority),rep(0,n.minority)))
  data_out = sample_3
  colnames(data_out) = c(colnames(data_out)[-(numvars:(numvars+2))],"Perf","Race_ind","Race_dummy")
  
  return(data_out)
  
}
