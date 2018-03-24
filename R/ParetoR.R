# Command Function (ParetoR) - Pareto-Optimization via Normal Boundary Intersection
# Developer: Temporarily removed for blind review
# Contact: Temporarily removed for blind review
# Last Update: 02/12/2018

###### Instruction ######

# 1. Open R console or RStudio
# 2. Install R package "ParetoR" through Github by running the following command in R console or RStudio:
#    install.packages("devtools")
#    library("devtools")
#    install_github("Diversity-ParetoOptimal/ParetoR")
# 3. Specify four input variables (details in "Input Description" below):
#    a) prop (i.e., proportion of minority applicants)
#    b) sr (i.e., selection ratio)
#    c) d (i.e., subgroup difference)
#    d) R (i.e., correlation matrix of predictor and criterions)
# 4. Run "out = ParetoR(prop, sr, R, d)" in R console or RStudio

## Input Description ##

# 1) Proportion of minority applicants (prop):
# prop = (# of minority applicants)/(total # of applicants)
# 2) Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)
# 3) Subgroup difference (d): criterion & predictor mean difference
# between minority and majority subgroups
# 4) Correlation matrix (R): criterion & predictor inter-correlation

## Output Description ##

# 1) Pareto Optimal solution (i.e., AI ratio, Criterion Validity, Predictor Weights);
# 2) Plots (i.e., AI ratio - Criterion Validity trade-off & predictor weights trade-off).

###### Example Input ######

# # DeCorte, Lievens & Sackett (2007) application example is used as example input below.
#
# # Proportion of minority applicants in full applicant pool
# prop <<- 1/4
#
# # Selection ratio
# sr <<- 0.10
#
# # Subgroup difference
# d <<- c(1.00, 0.23, 0.09, 0.33)
#
# # Correlation matrix
# # Format: Predictor_1, ..., Predictor_n, Criterion
# R <<- matrix(c(1, .24, .00, .19, .30,
#                .24, 1, .12, .16, .30,
#                .00, .12, 1, .51, .18,
#                .19, .16, .51, 1, .28,
#                .30, .30, .18, .28, 1),
#              (length(d)+1),(length(d)+1),byrow=T)
#
# ###### Additional Settings ######

#' ParetoR
#'
#' Command function to run Pareto-Optimal algorithm.
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @param d Subgroup difference
#' @param R Correlation matrix
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
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
#' # Fit Pareto-optimal model
#' out = ParetoR(prop, sr, d, R)
#'
#' @export
ParetoR = function(prop, sr, d, R, Spac = 20, graph = TRUE, display_solution = TRUE){

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

  out = NBI(X0, Spac, Fnum, VLB, VUB, TolX, TolF, TolCon, graph=graph, display_solution=display_solution)
  return(out)

}

