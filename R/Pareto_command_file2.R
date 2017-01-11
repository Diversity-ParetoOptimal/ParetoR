# Command File - Pareto-Optimization via Normal Boundary Intersection
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com
# Last Update: 11/18/2016

###### User Input ######

#### Instruction ####

# 1. Unzip 'Pareto_R.zip' directly under C drive (i.e., 'C:/') through option "Extract Here";
# 2. Open 'Pareto_command_file.r' in R;
# 3. Change the following 4 input variables (i.e., 'prop', 'sr', 'd', 'R')
# based on specific selection scenario. See 'Input Description' below
# for more details. You may also specify 'Additional Settings' if needed;
# 4. Select all and run 'Pareto_command_file.r'.

## Input Description ##

# 1) Population correlation matrix (R): criterion & predictor inter-correlation
# 2) Population subgroup difference (d): criterion & predictor mean difference
# between minority and majority subgroups
# 3) Proportion of minority applicants (prop):
# prop = (# of minority applicants)/(total # of applicants)
# 4) Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)

## Output Description ##

# 1) Pareto Optimal solution (i.e., AI ratio, Criterion Validity, Predictor Weights);
# 2) Plots (i.e., AI ratio - Criterion Validity trade-off & predictor weights trade-off).

###### Input ######

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
#' @return out Pareto-Optimal solution
#' @export
ParetoR = function(prop,sr,d,R){

  prop <<- prop
  sr <<- sr
  d <<- d
  R <<- R

  # Number of Pareto-Points
  Spac 	= 20

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

  out = NBI(X0,Spac,Fnum,VLB,VUB,TolX,TolF,TolCon)
  return(out)

}


# Pareto_Fmat = out$Pareto_Fmat
# Pareto_Xmat = out$Pareto_Xmat

