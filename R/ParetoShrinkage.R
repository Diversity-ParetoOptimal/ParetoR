### Pareto-Optimal Shrinkage Formulae ###

#' ParetoShrinkage function
#'
#' Estimate shrunken Pareto-optimal solution based on Pareto-optimal shrinkage formulae
#' @param Ncal Calibration sample size
#' @param numpred number of predictors
#' @param p_pareto number of Pareto-optimal points (e.g., 21)
#' @param R_perf_cal Vector of calibration sample job performance validity
#' @param R_race_cal Vector of calibration sample race bivariate correlation [i.e., correlation between race dummy variable (0-minority, 1-majority) and predictor composite score]
#' @param use Shrinkage formula to be used at first step when generalizing to population. Options for "use": "Wherry" - Wherry (1931) shrinkage formula; "OlkinPratt" - simplified Olkin & Pratt (1958) formula (see Cattin, 1980, p. 409), "OlkinPratt Hypergeo" - Olkin & Pratt (1958) formula with hypergeometric function; "Claudy" - Claudy (1978) formula. Default is "Wherry".
#' @return out formula-adjusted Pareto-optimal solution
#' @examples 
#' # (1) Calibration sample size
#' Ncal <- 100
#' # (2) Number of predictors
#' numpred <- 4
#' # (3) Number of Pareto-optimal points (i.e., number of sets of predictor weights)
#' p_pareto <- 21
#' # (4) Vector of calibration sample job performance validity
#' load(R_perf_cal)
#' # (5) Vector of calibration sample race bivariate correlation [i.e., correlation between race dummy variable (0-minority, 1-majority) and predictor composite score]
#' load(R_race_cal)
#' # Estimate shrunken Pareto-optimal solution
#' ParetoShrinkage(Ncal = Ncal, numpred = numpred, p_pareto = p_pareto, R_perf_cal = R_perf_cal, R_race_cal = R_race_cal, use = "Wherry)
#' @export

ParetoShrinkage <- function(Ncal, numpred, p_pareto, R_perf_cal, R_race_cal, use= "Wherry"){
  
  # Create space for output
  cnames = c("R_perf_cal","R_perf_adj","R_race_cal","R_race_adj","w_perf","w_race")
  out = matrix(NA, p_pareto, length(cnames))
  colnames(out) = cnames
  
  out[, "R_perf_cal"] <- R_perf_cal
  out[, "R_race_cal"] <- R_race_cal
  
  # R2_perf where performance was maximized
  R_perf_max = R_perf_cal[p_pareto]
  # R2_perf where diversity was maximized
  R_perf_min = R_perf_cal[1] 
  # R2_race where diversity was maximized
  R_race_max = R_race_cal[1]
  # R2_race where performance was maximized
  R_race_min = R_race_cal[p_pareto]
  
  
  # Obtain shrunken R2 values of endpionts via Browne shrinkage formula
  R_perf_adj_max = sqrt(R2_Browne(N = Ncal, p = (numpred), R2 = (as.numeric(R_perf_max))^2, use = use))
  if(R_race_max>=0){
    R_race_adj_max = sqrt(1-(R2_Browne(N = Ncal, p = (numpred), R2 = (1-(as.numeric(R_race_max))^2), use=use)))
  }else{
    R_race_adj_max = -sqrt(R2_Browne(N = Ncal, p = (numpred), R2 = (as.numeric(R_race_max))^2, use=use))
  }
 
  for(l in 1:p_pareto){
    
    R_perf_p = R_perf_cal[l]
    R_race_p = R_race_cal[l]
    
    ## obtain criterion weights "w_perf", "w_race" ## 
    
    w_perf = (R_race_p - R_race_max)/(R_race_min - R_race_max) 
    w_race = (R_perf_p - R_perf_max)/(R_perf_min - R_perf_max)
    out[l, "w_perf"] <- w_perf
    out[l, "w_race"] <- w_race
    
    ## Obtain "R_perf_Bh", "R_race_Bh" ##
    
    # Performance
    R_perf_adj = (w_race*(R_perf_min - R_perf_adj_max)) + R_perf_adj_max
    # Diversity
    R_race_adj = (w_perf*(R_race_min - R_race_adj_max)) + R_race_adj_max
    
    out[l, "R_perf_adj"] <- R_perf_adj
    out[l, "R_race_adj"] <- R_race_adj 
    
  }
  
  return(out)
  
}

###### Wherry formula ######
# generalizing to population

#' R2_Wherry function
#'
#' Estimate shrunken R2 based on Wherry (1931) formula
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @return R2_W formula-adjusted R2 based on Wherry (1931) shrinkage formula
#' @examples 
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_Wherry(N = N, p = p, R2 = R2)
#' @export
R2_Wherry = function(N, p, R2){
  
  # N: sample size
  # p: number of predictors
  # R2: calibration sample R-squared
  
  R2_W = 1-((N-1)/(N-p))*(1-R2)
  R2_W[R2_W<0]=0 # set to zero if output value is less than zero
  
  return(R2_W)
}

###### Olkin & Pratt shrinkage formula ######
# generalizing to population

#' R2_OP function
#'
#' Estimate shrunken R2 based on simplified Olkin & Pratt (1958) formula (see Cattin, 1980, p. 409)
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @return R2_OP formula-adjusted R2 based on simplified Olkin & Pratt (1958) shrinkage formula (see Cattin, 1980, p. 409)
#' @examples 
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_OP(N = N, p = p, R2 = R2)
#' @export
R2_OP = function(N, p, R2){
  
  # N: sample size
  # p: number of predictors
  # R2: calibration sample R-squared
  
  R2_OP = 1-((N-3)/(N-p-1))*(1-R2)*(1+(2*(1-R2)/(N-p+1))+((8*(1-R2)^2)/((N-p-1)*(N-p+3))))
  R2_OP[R2_OP<0]=0 # set to zero if output value is less than zero
  
  return(R2_OP)
}

# Olkin & Pratt formula with hypergeometric function

#' R2_OPh function
#'
#' Estimate shrunken R2 based on Olkin & Pratt (1958) formula with hypergeometric function
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @return R2_OPh formula-adjusted R2 based on Olkin & Pratt (1958) shrinkage formula with hypergeometric function
#' @examples 
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_OPh(N = N, p = p, R2 = R2)
#' @import hypergeo
#' @export
R2_OPh = function(N, p, R2){
  
  require("hypergeo")
  
  # N: sample size
  # p: number of predictors
  # R2: calibration sample R-squared
  
  R2_OPh = 1 - (N-2)/(N-p)*(1-R2)*(to_real(hypergeo(A=1, B=1, C=(N-p+2)/2, z=1-R2))[1])
  R2_OPh[R2_OPh<0]=0 # set to zero if output value is less than zero
  
  return(R2_OPh)
}

###### Claudy formula ######
# generalizing to population

#' R2_Claudy function
#'
#' Estimate shrunken R2 based on Claudy (1978) formula 
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @return R2_C formula-adjusted R2 based on Claudy (1978) shrinkage formula
#' @examples 
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_Claudy(N = N, p = p, R2 = R2)
#' @import hypergeo
#' @export
R2_Claudy = function(N, p, R2){
  
  # N: sample size
  # p: number of predictors
  # R2: calibration sample R-squared
  
  R2_C = 1-((N-4)*(1-R2))/(N-p-1)*(1+2*(1-R2)/(N-p+1))
  R2_C[R2_C<0]=0 # set to zero if output value is less than zero
  
  return(R2_C)
}

##### Browne shrinkage formula ######
# generalizing to population

#' R2_Browne function
#'
#' Estimate shrunken R2 based on Browne (1975) formula 
#' @param N Sample size
#' @param p number of predictors
#' @param R2 R-squared
#' @param use Shrinkage formula to be used at first step when generalizing to population. Options for "use": "Wherry" - Wherry (1931) shrinkage formula; "OlkinPratt" - simplified Olkin & Pratt (1958) formula (see Cattin, 1980, p. 409), "OlkinPratt Hypergeo" - Olkin & Pratt (1958) formula with hypergeometric function; "Claudy" - Claudy (1978) formula. Default is "Wherry".
#' @return R2_B formula-adjusted R2 based on Browne (1958) shrinkage formula
#' @examples 
#' # (1) Sample size
#' N <- 100
#' # (2) Number of predictors
#' p <- 5
#' # (3) R2 R-squared
#' R2 <- 0.30
#' # Estimate shrunken R2
#' R2_Browne(N = N, p = p, R2 = R2, use = "Wherry")
#' @import hypergeo
#' @export
R2_Browne = function(N, p, R2, use = "Wherry"){
  
  # N: sample size
  # p: number of predictors
  # R2: calibration sample R-squared
  
  # use 1st formula to generalize to population first
  if(use=="Wherry"){R2_pop = R2_Wherry(N,p,R2)}
  if(use=="OlkinPratt"){R2_pop = R2_OP(N,p,R2)}
  if(use=="OlkinPratt Hypergeo"){require('hypergeo'); R2_pop = R2_OPh(N,p,R2)}
  if(use=="Claudy"){R2_pop = R2_Claudy(N,p,R2)}
  
  # then, use Browne formula to generalize to next sample
  R2_B = ((N-p-3)*R2_pop^2+R2_pop)/((N-2*p-2)*R2_pop+p)
  R2_B[R2_B<0]=0 # set to zero if output value is less than zero
  
  return(R2_B)
}

#### Other Functions ####

## which.median() ##

# Find the index of median number
which.median <- function(x) {
  if (length(x) %% 2 != 0) {
    which(x == median(x))
  } else if (length(x) %% 2 == 0) {
    a = sort(x)[c(length(x)/2, length(x)/2+1)]
    c(which(x == a[1]), which(x == a[2]))
  }
}

## rpb2rb ##

# Convert point-biserial correlation to biserial correlation
rpb2rb <- function(r_pb, prop.maj, prop.min){
  r_b = r_pb*(sqrt(prop.maj*prop.min)/abs(dnorm(qnorm(prop.maj,0,1),0,1)))
  return(r_b)
}

## r2AI ##

# Convert r to AIratio
r2AI <- function(rrace,prop,sr){
  r <- rrace #r_race; point-biserial correlation
  p <- prop #proportion of Black applicants
  z_cut <- qnorm(sr) #standard normal cut score across all applicants 
  #(which is a direct transformation of the overall selection ratio, SR)
  
  #formula 1
  ## d: subgroup d (Black-White)
  d = r/sqrt((1-r^2)*p*(1-p))
  
  #formula 2
  x_cut = z_cut*sqrt(1+d^2*(p*(1-p))) - d*p
  
  #formula 3
  trans.AIratio = (1.64*x_cut+sqrt(0.76*(x_cut)^2+4))/(sqrt((exp(d^2+2*x_cut*d)))*(1.64*(x_cut+d)+sqrt(0.76*(x_cut+d)^2+4)))
  
  trans.AIratio[trans.AIratio<0] = 0
  trans.AIratio[trans.AIratio>2] = 0
  
  return(trans.AIratio)
}

## AI2r ##

AI2r <- function(AI, prop, sr, tol = 10^(-5)){
  
  r_temp = seq(-1,1,tol)
  AI_temp = r2AI(r_temp, prop, sr)
  trans.r = r_temp[which.min(abs(AI_temp - AI))]
  return(trans.r)
  
}

## wq ##
# Relative criterion weight in Pareto weighting
# (e.g., relative criterion weight related to AI ratio out of combined criterion of AI ratio and validity/performance)

wq <- function(c_p, c_max, c_min, prop, sr){
  # c_p: Criterion value at current Pareto point
  # c_max: Criterion value at endpoint where this Criterion was maximized
  # c_min: Criterion value at endpoint where the other criterion was maximized
  wq = (c_p - c_max)/(c_min - c_max)
  return(wq)
}

## AI2rwq ##
# Obtain relative criterion weight related to r_race from relative criterion weight for AI ratio

AI2rwq <- function(a_p,a_max,a_min,prop,sr){
  
  r_p = AI2r(a_p,prop,sr) # a_p: AI ratio at current Pareto point; r_p: r_race at current Pareto point
  r_max = AI2r(a_max,prop,sr) # a_max: AI ratio at endpoint where diversity was maximized; r_max: r_race at endpoint where diversity was maximized
  r_min = AI2r(a_min,prop,sr) # a_min: AI ratio at endpoint where performance was maximized; r_min: r_race at endpoint where performance was maximized
  r_wq = (r_p - r_max)/(r_min - r_max)
  return(r_wq)
  
}

## AI2rwq2 ##
# Obtain squared relative criterion weight related to r_race from relative criterion weight for AI ratio

AI2rwq2 <- function(a_p,a_max,a_min,prop,sr){
  
  r_p = AI2r(a_p,prop,sr) # a_p: AI ratio at current Pareto point; r_p: r_race at current Pareto point
  r_max = AI2r(a_max,prop,sr) # a_max: AI ratio at endpoint where diversity was maximized; r_max: r_race at endpoint where diversity was maximized
  r_min = AI2r(a_min,prop,sr) # a_min: AI ratio at endpoint where performance was maximized; r_min: r_race at endpoint where performance was maximized
  r_wq2 = (r_p^2 - r_max^2)/(r_min^2 - r_max^2)
  return(r_wq2)
  
}
