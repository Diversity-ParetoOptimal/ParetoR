# Pareto-Optimization via Normal Boundary Intersection
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com
# Last Update: 03/28/2017

####################################### NBI Main Function ####################################

#' NBI Main Function
#'
#' Main function for obtaining pareto-optimal solution via normal boundary intersection.
#' @param X0 Initial input for preditor weight vector
#' @param Spac Number of Pareto spaces (i.e., number of Pareto points minus one)
#' @param Fnum Number of criterions
#' @param VLB Lower boundary for weight vector estimation
#' @param VUB Upper boundary for weight vector estimation
#' @param TolX Tolerance index for estimating weight vector, default is 1e-4
#' @param TolF Tolerance index for estimating criterion, default is 1e-4
#' @param TolCon Tolerance index for constraint conditions, default is 1e-7
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @param display_solution If TRUE, Pareto-optimal solution will be displayed
#' @import nloptr
#' @return Pareto-Optimal solutions
#' @export
NBI = function(X0,Spac,Fnum,VLB=vector(),VUB=vector(),TolX=1e-4,TolF=1e-4,TolCon=1e-7,graph=TRUE,display_solution=TRUE){

cat('\n Estimating Pareto-Optimal Solution ... \n')

#------------------------------Initialize Options-------------------------------#

  X0 = assert_col_vec(X0)
  VLB = assert_col_vec(VLB)
  VUB = assert_col_vec(VUB)

  # Number of variables
  nvars = length(X0)

  # Set options
  # algorithm: sequential (least-squares) quadratic programming algorithm
  # (SQP is algorithm for nonlinearly constrained, gradient-based optimization,
  # supporting both equality and inequality constraints.)
  # maxeval: Max Iterations
  # xtol_abs: Tolerance for X
  # ftol_abs: Tolerance for F
  # tol_constraints_ineq/eq: Tolerance for inequal/equal constraints
  # (for reference) MATLAB constraints:
  # options = optimoptions('fmincon','Algorithm','sqp','MaxIter',(nvars+1)*1000,'TolFun',TolF,'TolX',TolX,'TolCon',TolCon,'Display','off')

  nloptr::nl.opts(optlist = list(
             maxeval = (nvars+1)*1000
             ,xtol_rel = TolX
             ,ftol_rel = TolF
             ))

  #Initialize PHI

  PHI = matrix(0,Fnum,Fnum)

  #------------------------------Shadow Point-------------------------------#

  # cat('\n ----Step 1: find shadow minimum---- \n')

  ShadowF = matrix(0,Fnum)
  ShadowX = matrix(0,nvars,Fnum)
  xstart  = X0
  out = WeightsFun(Fnum,Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight/Spac

  for(i in 1:Fnum){

    temp = c(1,dim(Weight)[2])
    j = temp[i]
    g_Weight <<- Weight[,j]
    fmin = 9999
    ntr = nvars-1
    fminv = matrix(0,ntr,1)
    fminx = matrix(0,nvars,ntr)

    for(k in 1:ntr){

      xstart = runif(length(X0))

      out = nloptr::slsqp(x0 = X0, fn = myLinCom
                 ,lower = VLB, upper = VUB
                 ,hin = myCon_ineq
                 ,heq = myCon_eq
                  )
      x = out$par
      f = out$value
      rm(out)

      fminv[k] = f
      fminx[,k] = x

      if(f <= fmin){

        fmin = f
        reps = k

      }

    }

    x = fminx[,reps]
    som = 0

    for(k in 2:nvars){
      som = som + x[k]
    }

    for(k in 2:nvars){
      x[k] = x[k]/som
    }
    # to make sum of x = 1

    ShadowX[,i] = x
    ShadowX = round(ShadowX,4)

    tempF = -myFM(x)
    ShadowF[i] = round(tempF[i],4)

  }

  # cat( '\n Shadow Minimum-F: \n')
  # print(round(ShadowF,3))
  # cat('\n Shadow Minimum--X(column) \n')
  # print(round(ShadowX,3))

  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for(i in 1:Fnum){

    PHI[,i] = myFM(ShadowX[,i]) + ShadowF
    PHI[i,i] = 0

  }

  # print(round(PHI,3))

  #Check to make sure that QPP is n-1 dimensional
  if(rcond(PHI) < 1e-8){stop(' Phi matrix singular, aborting.')}

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  g_Normal <<- -PHI%*%matrix(1,Fnum,1)

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun(Fnum,Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight/Spac
  num_w = dimFun(Weight)[2]

  # cat('\n Weights in row: \n')
  # print(round(Weight,3))

  #------------------------------NBI Subproblems-------------------------------#

  # cat('\n ----Step 5: solve NBI sub-problems---- \n')

  # Starting point for first NBI subproblem is the minimizer of f_1(x)
  xstart = c(ShadowX[,1],0)

  Pareto_Fmat = vector()       # Pareto Optima in F-space
  Pareto_Xmat = vector()       # Pareto Optima in X-space
  X_Near      = vector()

  # solve NBI subproblems
  for(k in 1:num_w){

    w  = Weight[,k]

    # Solve problem only if it is not minimizing one of the individual objectives
    indiv_fn_index = which(w == 1)
    # the boundary solution which has been solved

    if(length(indiv_fn_index) != 0){

      # w has a 1 in indiv_fn_index th component, zero in rest
      # Just read in solution from shadow data
      Pareto_Fmat = cbind(Pareto_Fmat, (-PHI[,indiv_fn_index] + ShadowF))
      Pareto_Xmat = cbind(Pareto_Xmat, ShadowX[,indiv_fn_index])
      X_Near = cbind(X_Near, c(ShadowX[,indiv_fn_index],0))
      # print(Pareto_Fmat)

    }else{

      w = rev(w)

      if(Near[k] > 0){

        xstart = X_Near[,Near[k]]
        #start X is the previous weight-order's X

      }

      #start point in F-space
      g_StartF <<- PHI%*%w + ShadowF

      # SOLVE NBI SUBPROBLEM

      out = nloptr::slsqp(x0 = xstart, fn = myT
                  ,lower = c(VLB,-Inf)
                  ,upper = c(VUB,Inf)
                  ,hin = myCon_ineq
                  ,heq = myTCon_eq)

      x_trial = out$par
      f = out$value
      rm(out)

      # success
      # if(fiasco >= 0){

        Pareto_Fmat = cbind(Pareto_Fmat, -myFM(x_trial[1:nvars]))  # Pareto optima in F-space
        som = 0

        for(k in 2:nvars){som = som + x_trial[k]}

        for(k in 2:nvars){x_trial[k] = x_trial[k]/som}

        Pareto_Xmat = cbind(Pareto_Xmat, x_trial[1:nvars])        # Pareto optima in X-space
        X_Near = cbind(X_Near,x_trial)

#       }else{
#
#         # unsuccess
#
#         num_fiascos = num_fiascos + 1
#         PHI2 = matrix(0,Fnum,Fnum)
#
#         for(i in 1:Fnum){
#
#           PHI2[,i] = myFM(x_trial) - ShadowF
#           PHI2[i,i] = 0
#         }
#
#         g_Normal2  = -PHI2%*%matrix(1,Fnum,1)
#         g_Normal2  = g_Normal2/norm(g_Normal2,type='2') # specifies 2-norm, or Euclidean length  the vector
#         Pareto_Fmat = c(Pareto_Fmat, (g_StartF + xstart[dimFun(xstart)[1]] %*% g_Normal2))  # Pareto optima in F-space
#         X_Near = c(X_Near,x_trial)
#         giveup = readline('Give up ?  (0/1)  ')
#         if(giveup == 1){break}
#
#       }

      }

    }

  #------------------------------Plot Solutions-------------------------------#

#   cat('\n ----Step 6: plot---- \n')

  if(graph==TRUE){plotPareto(Pareto_Fmat, Pareto_Xmat)}

  #------------------------------Output Solutions-------------------------------#

#   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat[2:nrow(Pareto_Xmat),])
  colnames(Pareto_Fmat) = c("AI.ratio","Criterion.Validity")
  colnames(Pareto_Xmat) = c(paste0(rep("P",(nvars-1)),1:(nvars-1)))
  
  if(display_solution == TRUE){
    
    solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
    colnames(solution) = c("AI.ratio","Criterion.Validity", paste0(rep("P",(nvars-1)),1:(nvars-1)))
    cat("\n Pareto-Optimal Solution \n \n")
    print(solution)
  
  }else{
    cat("\n Done. \n \n")
  }


  return(list(Pareto_Fmat = round(Pareto_Fmat, 3),
              Pareto_Xmat = round(Pareto_Xmat, 3)))

}

########################### Supporting Functions (A) ########################

# User-Defined Input for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Input:
## 1) Population correlation matrix (R): criterion & predictor inter-correlation
## 2) Population subgroup difference (d): criterion & predictor mean difference
## between minority and majority subgroups
## 3) Proportion of minority applicants (prop):
## prop = (# of minority applicants)/(total # of applicants)
## 4) Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)

# Related functions:
# myFM
# myCon

###### myFM() ######

#' myFM
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @export
myFM = function(x){

  # Obtain within-package 'global' variables
  d <- d_ParetoR
  R <- R_ParetoR

  R_u = R[-nrow(R),-ncol(R)]
  b = x[-1]

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar = 0
  # mean of majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar = d%*%x[-1]/sigma_p
  # minority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_B = 1 - pnorm(x[1], p_i_bar)
  # majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_W = 1 - pnorm(x[1], p_a_bar)

  # AIratio a_g (DeCorte et al., 2007)
  a_g = SR_B/SR_W

  # Composite Validity R_xy
  R_xy = t(c(t(b),0)%*%R%*%c(t(matrix(0,dimFun(R_u)[1],1)),1))/sqrt(t(b)%*%R_u%*%b) # DeCorte et al., 2007

  f = matrix(1,2,1)
  f[1,] = -a_g
  f[2,] = -R_xy

  return(f)

}

####### myCon_ineq() ######

# Nonlinear inequalities at x

#' myCon_ineq
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @export
myCon_ineq = function(x){return(vector())}

####### myCon_eq() ######

# Nonlinear equalities at x

#' myCon_eq
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @return Equal constraint condition for use in NBI()
#' @export
myCon_eq = function(x){

  # Obtain within-package 'global' variable
  prop <- prop_ParetoR
  sr <- sr_ParetoR
  d <- d_ParetoR
  R <- R_ParetoR

  R_u = R[-nrow(R),-ncol(R)]
  b = x[-1]

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar = 0
  # mean of majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar = d%*%x[-1]/sigma_p
  # p_a_bar = (x[2]*1.00+x[3]*0.23+x[4]*0.09+x[5]*0.33)/sigma_p
  # minority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_B = 1 - pnorm(x[1], p_i_bar)
  # majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_W = 1 - pnorm(x[1], p_a_bar)

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = SR_B*prop + SR_W*(1-prop) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%R_u%*%b) - 1

  return(ceq)

}

########################### Supporting Functions (B) ########################

# Supplementary Functions for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Function List
## assert_col_vec
## dimFun
## WeightsFun
## Weight_Generate
## myLinCom
## myT
## myTCon_eq
## plotPareto

###### assert_col_vec() ######

#' assert_col_vec
#'
#' Support function, refines intermediate variable for use in NBI()
#' @param v Intermediate variable v
#' @return Refined variable v
#' @export
assert_col_vec = function(v){
  if(is.null(dimFun(v))){
    v=v
  }else if(dimFun(v)[1] < dimFun(v)[2]){v = t(t)}
  return(v)}

###### dimFun() ######

#' dimFun
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @return x Checked and refined input predictor weight vector
#' @export
dimFun = function(x){
  if(is.null(dim(x))){
    return(c(0,0))
  }else(return(dim(x)))
}

###### WeightsFun() ######

#' WeightsFun
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
#' @export
WeightsFun = function(n, k){

  # global variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  WeightSub <<- matrix(0,1,n)
  Weights <<- vector()
  # assign("Formers", vector(), envir = .GlobalEnv)
  Formers <<- vector()
  # assign("Layer", n, envir = .GlobalEnv)
  Layer <<- n
  # assign("lastone", vector(), envir = .GlobalEnv)
  lastone <<- -1
  # assign("currentone", -1, envir = .GlobalEnv)
  currentone <<- -1

  Weight_Generate(1, k)

  return(list(Weights = Weights, Formers = Formers))

}

###### Weight_Generate() ######

#' Weight_Generate
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate
#' @export
Weight_Generate = function(n, k){

  # global variables:
  # weight Weights Formers Layer lastone currentone

  # wtgener_test(n,k)
  #
  # Intended to test the weight generation scheme for NBI for > 2 objectives
  # n is the number of objectives
  # 1/k is the uniform spacing between two w_i (k integral)

  if(n == Layer){

    if(currentone >= 0){
      Formers <<- c(Formers,lastone)
      lastone <<- currentone
      currentone <<- -1
    }else{
      num = dimFun(Weights)[2]
      Formers <<- c(Formers,num)
    }

    WeightSub[(Layer - n + 1)] <<- k
    Weights <<- cbind(Weights,t(WeightSub))

  }else{

    for(i in 0:k){
      if(n == (Layer - 2)){
        num = dimFun(Weights)[2]
        currentone <<- num+1
      }

      WeightSub[(Layer - n + 1)] <<- i
      Weight_Generate(n+1, k-i)
    }

  }

}

###### myLinCom() ######

#' myLincom
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @export
myLinCom = function(x){

  # global variable: g_Weight
  F   = myFM(x)
  f = t(g_Weight)%*%F
  return(f)

}

###### myT() ######

#' myT
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
#' @export
myT = function(x_t){

  f = x_t[length(x_t)]
  return(f)

}

###### myTCon_eq() ######

#' myTCon_eq
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @export
myTCon_eq = function(x_t){

  # global variables:
  # g_Normal g_StartF

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = -myFM(x) - g_StartF - t * g_Normal

  # c = myCon_ineq(x)
  ceq1 = myCon_eq(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

###### plotPareto() ######

#' plotPareto
#'
#' Function for plotting Pareto-optimal curve and predictor weights
#' @param CriterionOutput Pareto-Optimal criterion solution
#' @param ParetoWeights Pareto-Optimal predictor weight solution
#' @return Plot of Pareto-optimal curve and plot of predictor weights
#' @export
plotPareto = function(Pareto_Fmat, Pareto_Xmat){

  par(mfrow=c(1,2))

  AIratio = t(Pareto_Fmat[1,])
  Criterion = t(Pareto_Fmat[2,])
  X = t(Pareto_Xmat[2:nrow(Pareto_Xmat),])

  # AI ratio - Composite Validity trade-off

  plot(AIratio, Criterion,
       xlim = c(min(AIratio),max(AIratio)),
       main = "Composite Validity -- AI ratio trade-off",
       xlab = "AI ratio",
       ylab = "Composite Validity",
       type='c',col='blue')

  points(AIratio, Criterion,
         pch=8,col='red')

  # Predictor weights

  plot(AIratio,X[,1],
       xlim=c(min(AIratio),max(AIratio)),ylim=c(0,1),
       main = "Predictor weights trade-off function",
       xlab = "AI ratio",
       ylab = "Predictor weight",
       type='c',col='red')
  points(AIratio,X[,1],pch=8,
         col=rainbow(1))

  for(i in 2:ncol(X)){

    lines(AIratio,X[,i],type='c',
          col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))
    points(AIratio,X[,i],pch=8,
           col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))

  }

  legend('topleft',
         legend=c(paste0('Predictor ',1:ncol(X))),
         lty=c(rep(2,ncol(X))),lwd=c(rep(2,ncol(X))),
         col=rainbow(ncol(X)))

}
