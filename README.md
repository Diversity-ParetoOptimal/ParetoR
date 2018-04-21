# ParetoR

Pareto-Optimization via Normal Boundary Intersection Method in Diversity Hiring <br />
Developer: Q. Chelsea Song <br />
Contact: qianqisong@gmail.com <br />
Last Update: 03/23/2018 

## Objective ##

The current R package provides a set of Pareto-optimal solutions that simultaneously optimize both diversity and criterion validity in a personnel selection scenario [see Song, Wee, & Newman (2017). The current package allows for implementation of (1) Pareto-optimal method that was adapted from De Corte, Lievens & Sackett (2007); (2) Pareto-optimal shrinkage formulae to estimate formula-adjusted shrunken Pareto-optimal solutions (see Study 2 of Song (2018; dissertation); (3) regularized Pareto-optimal method (see Study 3 of Song (2018; dissertation)). 

## Instructions ##

### Install and Load Package ###

1. Open an R console or RStudio window. (R can be downloaded for free from https://cran.r-project.org; RStudio can be downloaded for free from https://www.rstudio.com/)
2. Install R package "ParetoR" through Github by pasting and running the following commands in R console or RStudio:
   install.packages("devtools") <br />
   library("devtools") <br />
   install_github("Diversity-ParetoOptimal/ParetoR") <br />
   library("ParetoR") <br />

### Main Functions ###

**ParetoR function**  <br />
*Pareto-optimal method introduced by De Corte, Lievens & Sackett (2007)*  <br /> 

#### Example Implementation ####

1. Specify four inputs (example from DeCorte, Lievens & Sackett (2007) is given below): <br />
   &nbsp; # (1) Proportion of minority applicants (**prop**) = (# of minority applicants)/(total # of applicants) <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; prop <- 1/4 <br />
   &nbsp; # (2) Selection ratio (**sr**) = (# of selected applicants)/(total # of applicants) <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; sr <- 0.10 <br />
   &nbsp; # (3) Subgroup differences (**d**): standardized mean differences between minority and majority subgroups, on each predictor (in applicant pool) <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp;  d <- c(1.00, 0.23, 0.09, 0.33) <br />
   &nbsp; # (4) Correlation matrix (**R**) = criterion & predictor inter-correlation matrix (in applicant pool) <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; # Format: Predictor_1, ..., Predictor_n, Criterion <br />
&nbsp; &nbsp; &nbsp; &nbsp; R <- matrix(c(1, .24, .00, .19, .30, <br /> 
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .24, 1, .12, .16, .30, <br /> 
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .00, .12, 1, .51, .18, <br /> 
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .19, .16, .51, 1, .28, <br /> 
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .30, .30, .18, .28, 1), <br /> 
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; (length(d)+1),(length(d)+1)) <br /><br />
2. Paste and run the following command in R console or RStudio: <br />
&nbsp; &nbsp; &nbsp; &nbsp; out = ParetoR(prop, sr, d, R)

#### Output Description ####

1. Pareto Optimal solutions (i.e., 21 equally-spaced solutions that characterize the Criterion validity – AI ratio tradeoff curve, and Predictor Weights at each point along tradeoff curve).
2. Plots (i.e., Criterion validity – AI ratio tradeoff curve, Predictor weights across trade-off points).

**ParetoShrinkage** function <br />
*Estimate shrunken Pareto-optimal solution based on Pareto-optimal shrinkage formulae introduced in Study 2 of Song (2018; dissertation)* <br />

#### Example Implementation ####

1. Specify inputs <br />
 &nbsp; # (1) Calibration sample size  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; Ncal <- 100 <br />
 &nbsp; # (2) Number of predictors  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; numpred <- 4 <br />
 &nbsp; # (3) Number of Pareto-optimal points (i.e., number of sets of predictor weights)  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; p_pareto <- 21 <br />
 &nbsp; # (4) Vector of calibration sample job performance validity  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; load(R_perf_cal) <br />
 &nbsp; # (5) Vector of calibration sample race bivariate correlation [i.e., correlation between race dummy variable (0-minority, 1-majority) and predictor composite score]  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; load(R_race_cal) <br />

2. Paste and run the following command in R console or RStudio: <br />
 &nbsp; # Estimate shrunken Pareto-optimal solution <br />
 &nbsp; &nbsp; &nbsp; &nbsp; ParetoShrinkage(Ncal = Ncal, numpred = numpred, p_pareto = p_pareto,  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; R_perf_cal = R_perf_cal, R_race_cal = R_race_cal,  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; use = "Wherry) <br />

#### Output Description ####

1. Formula-adjusted Pareto-optimal solution

**ParetoElnet** function <br />
*Regularized Pareto-optimal method introduced in Study 3 of Song (2018; dissertation)* <br />

#### Example Implementation ####

1. Specify inputs (example from De Corte, Lievens & Sackett (2007) is given below): <br />
 &nbsp; # (1) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor and criterion (in applicant pool) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; d <- c(1.00, 0.23, 0.09, 0.33) <br />
 &nbsp; # (2) Correlation matrix (R) = Criterion predictor inter-correlation matrix <br />
 &nbsp; # Format: Predictor_1, ..., Predictor_n, Criterion <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; R <- matrix(c(1, .24, .00, .19, .30, <br />
 &nbsp; &nbsp; &nbsp; &nbsp;           .24, 1, .12, .16, .30, <br />
 &nbsp; &nbsp; &nbsp; &nbsp;           .00, .12, 1, .51, .18, <br />
 &nbsp; &nbsp; &nbsp; &nbsp;           .19, .16, .51, 1, .28, <br />
 &nbsp; &nbsp; &nbsp; &nbsp;           .30, .30, .18, .28, 1), <br />
 &nbsp; &nbsp; &nbsp; &nbsp;          (length(d)+1),(length(d)+1)) <br />
 &nbsp; # (3) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; prop <- 1/4 <br />
 &nbsp; # (4) Selection ratio (sr) = (# of selected applicants)/(total # of applicants) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; sr <- 0.10	<br />
 &nbsp; # (5) Spac = number of Pareto points <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; Spac <- 21 <br />

2. Paste and run the following command in R console or RStudio:  <br />
 &nbsp; # Fit Regularized Pareto-optimal model with parameter selection via cross-validation <br />
 &nbsp; &nbsp; &nbsp; &nbsp; out = ParetoElnet(d = d, R = R, prop = prop, sr = sr, Spac = Spac) <br />

#### Output Description ####

1. Regularized Pareto-optimal solutions (i.e., 21 equally-spaced solutions that characterize the Criterion validity – AI ratio trade-off curve, and Predictor Weights at each point along trade-off curve).
2. Plots (i.e., Criterion validity – AI ratio trade-off curve, and Predictor weights across trade-off points).


**cv.ParetoElnet** function <br />
*Regularized Pareto-optimal method introduced in Study 3 of of Song (2018; dissertation) with parameter selection based on cross-validation method* <br />

#### Example Implementation ####

Example: Statistics of calibration sample data as input (i.e., calibration sample size, standardized subgroup mean difference, predictor and criterion correlation matrix). Input do not include calibration sample raw data set. <br /> 

1. Specify inputs (example from De Corte, Lievens & Sackett (2007) is given below): <br />
 &nbsp; # (1) Calibration sample size <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; n_cal = 100 <br />
 &nbsp; # (2) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor and criterion (in applicant pool) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; d <- c(1.00, 0.23, 0.09, 0.33, .30) <br />
 &nbsp; # (3) Correlation matrix (R) = Criterion predictor inter-correlation matrix  <br />
 &nbsp; # Format: Predictor_1, ..., Predictor_n, Criterion <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; R <- matrix(c(1, .24, .00, .19, .30, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;       .24, 1, .12, .16, .30, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;       .00, .12, 1, .51, .18, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;       .19, .16, .51, 1, .28, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;   	   .30, .30, .18, .28, 1), <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 	   (length(d)+1),(length(d)+1)) <br />
 &nbsp; # (4) Validation sample size <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; n_val = 10000 <br />
 &nbsp; # (5) Grid of alpha values to try <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; alpha.grid <- seq(0,1,length=3) <br />
 &nbsp; # (6) Grid of lambda values to try <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; lambda.grid <- 10^seq(1,-2,length=11) <br />
 &nbsp; # (7) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp;  prop <- 1/4 <br />
 &nbsp; # (8) Selection ratio (sr) = (# of selected applicants)/(total # of applicants) <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp;  sr <- 0.10 <br />
 &nbsp; # (9) Spac = number of Pareto points <br />
 &nbsp; ## Example <br />
 &nbsp; &nbsp; &nbsp; &nbsp; Spac <- 21 <br />

2. Paste and run the following command in R console or RStudio:
 &nbsp; # Fit Regularized Pareto-optimal model with parameter selection via cross-validation <br />
 &nbsp; &nbsp; &nbsp; &nbsp;  cv.out = cv.ParetoElnet(n_cal = n_cal, D = D, R = R,  <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;     		n_val = n_val, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;     		lambda.grid = lambda.grid, alpha.grid = alpha.grid, <br />
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;      		prop = prop, sr = sr, Spac = Spac) <br />

#### Output Description ####

1. Regularized Paret-optimal solutions (i.e., 21 equally-spaced solutions that characterize the Criterion validity – AI ratio trade-off curve, and Predictor Weights at each point along trade-off curve) with the model using best parameters selected using cross-validation.
2. Selected model parameters (i.e., alpha parameter and lambda parameter)
3. Plots (i.e., Criterion validity – AI ratio trade-off curve, and Predictor weights across trade-off points).

#### Note ####

The program is updated based on the ParetoR package that was introduced in Song et al. (2017). It is partially modeled after De Corte's (2006) TROFSS Fortran program and Zhou's (2006) NBI Matlab program (version 0.1.3). The current version only supports scenarios where AI ratio and one other criterion are being optimized.

#### References ####

Song, Q. C., Wee, S., & Newman, D. (provisionally accepted). Diversity Shrinkage: Cross-Validating Pareto-Optimal Weights to  Enhance Diversity via Hiring Practices. *Journal of Applied Psychology*. <br />
Das, I., & Dennis, J. E. (1998). Normal-boundary intersection: A new method for generating the Pareto surface in nonlinear multicriteria optimization problems. *SIAM Journal on Optimization*, 8, 631-657. <br />
De Corte, W. (2006). *TROFSS User's Guide*. <br />
De Corte, W., Lievens, F., & Sackett, P. (2007). Combining predictors to achieve optimal trade-offs between selection quality and adverse impact. *Journal of Applied Psychology*, 92, 1380-1393. <br />
Wee, S., Newman, D. A., & Joseph, D. L. (2014). More than g: Selection quality and adverse impact implications of considering second-stratum cognitive abilities. *Journal of Applied Psychology*, 99, 547-563. <br />

#### Acknowledgements ####

Great appreciation to Dr. Serena Wee, Dr. Dan Newman, Dr. Wilfred De Corte, and Dr. Victoria Stodden for guidance and feedback on development of the program.

#### Web Application ####

We also developed a user-friendly web application to implement the Pareto-Optimal technique described in the current package (https://qchelseasong.shinyapps.io/ParetoR/). The web application (like the ParetoR package) uses only a correlation matrix, selection ratio, proportion of applicants from the minority group, and subgroup d values as input. It then provides a full set of Pareto solutions and their corresponding predictor weights.
