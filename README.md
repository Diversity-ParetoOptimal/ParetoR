# ParetoR

Pareto-Optimization via Normal Boundary Intersection Method
Developer: Temporarily removed for blind review
Contact: Temporarily removed for blind review
Last Update: 01/11/2017

## Instructions ##

1. Open R console or RStudio
2. Install R package "ParetoR" through Github by running the following commands in R console or RStudio:
   install.packages("devtools") <br />
   library("devtools") <br />
   install_github("Diversity-ParetoOptimal/ParetoR") <br />
   library("ParetoR") <br />
3. Specify four inputs (details in "Input Description" and example in "Example Input" below): <br />
   a) prop # proportion of minority applicants <br />
   b) sr # selection ratio <br />
   c) d # subgroup difference <br />
   d) R # correlation matrix of predictor and criterions <br /> 
4. Run "out = ParetoR(prop, sr, R, d)" in R console or RStudio

## Example Input ##

DeCorte, Lievens & Sackett (2007) application example is used as example input below.

#### Proportion of minority applicants in full applicant pool
prop <- 1/4

#### Selection ratio
sr <- 0.10

#### Predictor subgroup mean difference
d <- c(1.00, 0.23, 0.09, 0.33)

#### Correlation matrix
####### Format: Predictor_1, ..., Predictor_n, Criterion
R <- matrix(c(1, .24, .00, .19, .30, <br /> 
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .24, 1, .12, .16, .30, <br /> 
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .00, .12, 1, .51, .18, <br /> 
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .19, .16, .51, 1, .28, <br /> 
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; .30, .30, .18, .28, 1), <br /> 
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; (length(d)+1),(length(d)+1)) 
             
## Input Description ##

1. Proportion of minority applicants (prop):
prop = (# of minority applicants)/(total # of applicants)
2. Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)
3. Subgroup difference (d): criterion & predictor mean difference
between minority and majority subgroups
4. Correlation matrix (R): criterion & predictor inter-correlation

## Output Description ##

1. Pareto Optimal solution (i.e., AI ratio, Criterion Validity, Predictor Weights);
2. Plots (i.e., AI ratio - Criterion Validity trade-off & predictor weights trade-off).

#### Note ####

The program is modeled after DeCorte's (2006) TROFSS Fortran program Zhou (2006)'s NBI Matlab program (version 0.1.3).
Current version only supports scenario where AI ratio and one other criterion is being maximized.

#### References ####

1. Das, I., & Dennis, J. E. (1998). Normal-boundary intersection: A new method for generating the Pareto surface in nonlinear multicriteria optimization problems. SIAM Journal on Optimization, 8, 631-657.
2. De Corte, W. (2006). TROFSS User's Guide.
3. De Corte, W., Lievens, F., & Sackett, P. (2007). Combining predictors to achieve optimal trade-offs between selection quality and adverse impact. Journal of Applied Psychology, 92, 1380-1393. 

