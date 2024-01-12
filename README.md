# ParetoR

Pareto-Optimization via Normal Boundary Intersection Method in Diversity Hiring <br />
Developer: Q. Chelsea Song <br />
Contact: qianqisong@gmail.com <br />
Last Update: 03/23/2023

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

**ParetoAdj function**  <br />
*Adjusts Pareto-optimal solutions using the shrinkage formula introduced by Song, Tang, Newman, & Wee (2023)*  <br /> 

#### Example Implementation

1.  Specify six inputs (example from Song, Tang, Newman, & Wee (2023) is given below): <br />
   &nbsp; # (1)  Calibration sample size (**Ncal**) <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; Ncal <- 40 <br />
   &nbsp; # (2) Predictor weights for each Pareto solution <br />
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; # Rows: Pareto solutions; Columns: Predictors <br />
&nbsp; &nbsp; &nbsp; &nbsp; wpred <- matrix(c(0,    0,    0,    0,    1,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.07, 0.93,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.13, 0.87,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.18, 0.82,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.23, 0.77,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.28, 0.72,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.32, 0.68,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.37, 0.63,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,    0,    0, 0.41, 0.59,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.01,    0, 0.45, 0.55,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.04,    0, 0.44, 0.52,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.07,    0, 0.43,  0.5,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,  0.1,    0, 0.42, 0.48,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.14,    0, 0.41, 0.46,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.17,    0,  0.4, 0.43,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.21,    0, 0.38, 0.41,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.25,    0, 0.37, 0.38,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0,  0.3,    0, 0.35, 0.35,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.33, 0.03, 0.34,  0.3,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.37, 0.08, 0.31, 0.25,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 0, 0.42, 0.15, 0.27, 0.15), ncol = 5, nrow = 21)
   &nbsp; # (3) Vector of calibration sample job performance validity <br />
      &nbsp; ## *Example*: <br />
&nbsp; &nbsp; &nbsp; &nbsp; Rperf_cal = c(.20, .24, .27, .30, .33, .36, .39, .42, .45, .48, .51, .54,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .57, .60, .63, .66, .69, .72, .74, .76, .78)
   &nbsp; # (4) Vector of calibration sample race bivariate correlation <br />
      &nbsp; &nbsp; &nbsp; &nbsp; # (i.e., correlation between race dummy variable (0-minority, 1-majority) <br />
      &nbsp; &nbsp; &nbsp; &nbsp; # and predictor composite score)
      &nbsp; ## *Example*: <br />
&nbsp; &nbsp; &nbsp; &nbsp; Rrace_cal = c(-.12, -.11, -.10, -.10, -.09, -.08, -.07, -.06, -.05, -.03,
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; -.02, .00, .01, .03, .05, .07, .09, .12, .15, .19, .24)
   &nbsp; # (5) proportion of minority
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; prop <- 1/6
   &nbsp; # (6) selection ratio
      &nbsp; ## *Example*: <br />
      &nbsp; &nbsp; &nbsp; &nbsp; sr <- .15

2. Paste and run the following command in R console or RStudio: <br />
   &nbsp; # Estimate shrunken Pareto-optimal solution
&nbsp; &nbsp; &nbsp; &nbsp; out <- ParetoAdj(Ncal = Ncal, wpred = wpred, Rperf_cal = Rperf_cal, Rrace_cal = Rrace_cal, prop = prop, sr = sr)

#### Output Description

1. The criterion validity and AI ratio results of the Pareto-optimal trade-off curve that is adjusted using the shrinkage formula


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
