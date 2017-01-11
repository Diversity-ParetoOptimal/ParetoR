# ParetoR

Pareto-Optimization via Normal Boundary Intersection Method
Developer: Q. Chelsea Song
Contact: qianqisong@gmail.com
Last Update: 11/18/2016

## Objective ##

The current R program provides a set of Pareto-optimal solution that simultaneously optimize both 
diversity and criterion validity in a personnel selection scenario (see Song, Wee & Newman (under review);
Wee, Newman & Joseph (2014); De Corte, Lievens & Sackett (2007) for more details) 
Pareto-optimal solution is estimated using Normal-Boundary Intersection method (Das & Dennis, 1996).

## Instruction ##

1. Unzip 'Pareto_R.zip' directly under C drive (i.e., 'C:/') using option "Extract Here";
2. Open 'Pareto_command_file.r' in R;
3. Change the following 4 input variables (i.e., 'prop', 'sr', 'd', 'R')
based on specific selection scenario. See 'Input Description' below
for more details. You may also specify 'Additional Settings' if needed;
4. Select all and run 'Pareto_command_file.r'.

## Input Description ##

1) Population correlation matrix (R): criterion & predictor inter-correlation
2) Population subgroup difference (d): criterion & predictor mean difference 
between minority and majority subgroups
3) Proportion of minority applicants (prop): 
prop = (# of minority applicants)/(total # of applicants)
4) Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)

## Output Description ##

1) Pareto Optimal solution (i.e., AI ratio, Criterion Validity, Predictor Weights);
2) Plots (i.e., AI ratio - Criterion Validity trade-off & predictor weights trade-off).

#### Note ####

The program is modeled after DeCorte's (2006) TROFSS Fortran program and XXX's NBI MATLAB program.
Current version only supports scenario where AI ratio and one other criterion is being maximized.

#### References ####

Das, I., & Dennis, J. E. (1998). Normal-boundary intersection: A new method for generating the Pareto surface in nonlinear multicriteria optimization problems. SIAM Journal on Optimization, 8, 631-657.
De Corte, W. (2006). TROFSS User's Guide.
De Corte, W., Lievens, F., & Sackett, P. (2007). Combining predictors to achieve optimal trade-offs between selection quality and adverse impact. Journal of Applied Psychology, 92, 1380-1393. 

