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

1. Pareto Optimal solutions (i.e., 21 equally-spaced solutions that characterize the Criterion validity – AI ratio tradeoff curve, and Predictor Weights at each point along tradeoff curve)
2. Plots (i.e., Criterion validity – AI ratio tradeoff curve, Predictor weights across trade-off points)

**ParetoAdj function**  <br />
*Adjusts Pareto-optimal solutions using the shrinkage formula introduced by Song, Tang, Newman, & Wee (2023)*  <br /> 

#### Output Description

1. The criterion validity and AI ratio results of the Pareto-optimal trade-off curve that is adjusted using the shrinkage formula

**ParetoElnet function**  <br />
*Estimates Pareto-optimal solutions using the regularized tradeoff curve algorithm introduced by Song, Tang, Newman, & Wee (2023)*  <br /> 

#### Output Description

1. Predictor weights at each point along the tradeoff curve

**cvParetoElnet function**  <br />
*Estimates Pareto-optimal solutions using the regularized tradeoff curve algorithm introduced by Song, Tang, Newman, & Wee (2023); implements hyperparameter tuning*  <br /> 

#### Output Description

1. Predictor weights at each point along the tradeoff curve
2. lambda and alpha parameter values chosen through hyperparameter tuning

#### Note ####

The program is updated based on the ParetoR package that was introduced by Song et al. (2017). It is partially modeled after De Corte's (2006) TROFSS Fortran program and Zhou's (2006) NBI Matlab program (version 0.1.3). The current version only supports scenarios where AI ratio and one other criterion are being optimized.

#### References ####

Song, Q. C., Wee, S., & Newman, D. (2017). Diversity shrinkage: Cross-validating Pareto-optimal weights to enhance diversity via hiring practices. *Journal of Applied Psychology*. <br />
Song, Q. C., Tang, C., Newman, D. A., & Wee, S. (2023). Adverse impact reduction and job performance optimization via Pareto-optimal weighting: A shrinkage formula and regularization technique using machine learning. *Journal of Applied Psychology*, 108(9), 1461–1485.  <br />
Das, I., & Dennis, J. E. (1998). Normal-boundary intersection: A new method for generating the Pareto surface in nonlinear multicriteria optimization problems. *SIAM Journal on Optimization*, 8, 631-657. <br />
De Corte, W. (2006). *TROFSS User's Guide*. <br />
De Corte, W., Lievens, F., & Sackett, P. (2007). Combining predictors to achieve optimal trade-offs between selection quality and adverse impact. *Journal of Applied Psychology*, 92, 1380-1393. <br />
Wee, S., Newman, D. A., & Joseph, D. L. (2014). More than g: Selection quality and adverse impact implications of considering second-stratum cognitive abilities. *Journal of Applied Psychology*, 99, 547-563. <br />

#### Acknowledgements ####

Great appreciation to Dr. Serena Wee, Dr. Dan Newman, Dr. Chen Tang, Dr. Wilfred De Corte, and Dr. Victoria Stodden for guidance and feedback on the development of the program.

#### Web Application ####

We also developed a user-friendly web application to implement the Pareto-Optimal technique described in the current package (https://qchelseasong.shinyapps.io/ParetoR/). The web application (like the ParetoR package) uses only a correlation matrix, selection ratio, proportion of applicants from the minority group, and subgroup d values as input. It then provides a full set of Pareto solutions and their corresponding predictor weights.
