

TruncGaussianLHS = function(samplesize,mean,sd,minval,maxval,varnames=NULL)
{
  N = length(mean)
  if(N < 1)
  {
    cat("error mean should have length > 0\n")
    return()
  }
  if(length(sd)!=N)
  {
    cat("error mean and sd should have equal length\n")
    return()
  }
  if(is.null(varnames))
  {
    varnames=paste("X",1:N,sep="")
  }
  if(length(varnames)!=N)
  {
    cat("error names and mean should have equal length\n")
    return()
  }
  result = data.frame(improvedLHS(samplesize,N))
  names(result) = varnames
  for(k in 1:N)
  {
    result[,k] =  qtruncnorm(result[,k], minval[k], maxval[k], mean[k], sd[k]) 
  }
  return(result)
}


##   The Wageningen Lowland Runoff Simulator (WALRUS): 
##   a lumped rainfall-runoff model for catchments with shallow groundwater
##   
##   Copyright (C) 2014 Claudia Brauer
##   
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##   
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##   
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Calibration function.
#' @description Calibrates model parameters. 
#' You have to run \code{WALRUS_preprocessing} and 
#' \code{WALRUS_preprocessing_calibration} before this function.
#' @param name name for the calibration (not optional; no default)
#' @param samplesize (default: 1000)
#' @param selectionsize (default: 10)
#' @param cW first estimate for cW (default: 250)
#' @param cG first estimate for cG (default: 20e6)
#' @param cQ first estimate for cQ (default: 20)
#' @param cV first estimate for cV (default: 10)
#' @param cS first estimate for cS (default: 4)
#' @param cD first estimate for cD (default: 1500)
#' @param fit_cW should cW be calibrated? (default: TRUE)
#' @param fit_cG should cG be calibrated? (default: TRUE)
#' @param fit_cQ should cQ be calibrated? (default: TRUE)
#' @param fit_cV should cV be calibrated? (default: TRUE)
#' @param fit_cS should cS be calibrated? (default: TRUE)
#' @param fit_cD should cD be calibrated? (default: TRUE)
#' @param aS surface water area fraction (default: 0.01)
#' @param st soil type (not optional; no default)
#' @param Gfrac fraction of discharge originating form quickflow at t=0 (default:0.5)
#' @return a data frame with the model output for all output time steps.
#' @export WALRUS_calibrate
#' @examples
#' x=1
#' 
WALRUS_calibrate = function(name, samplesize=1000, selectionsize=10,
  cW=250, cG=20e6, cQ=20, cV=10, cS=4, cD=1500,
  fit_cW=T, fit_cG=T, fit_cQ=T, fit_cV=F, fit_cS=F, fit_cD=F, 
  aS=0.01, st, Gfrac=0.5)
  {

  
  #########################################################
  ### Making parameter sets with Latin Hypercube Sampling
  #########################################################
  
  # Put parameters to calibrate in vector (order: cW, cG, cQ, cV, cS, cD)
  fit_TF = c(fit_cW,fit_cG,fit_cQ,fit_cV,fit_cS,fit_cD)
  fit_T  = c(1:6)[fit_TF]
  
  # Means and variances for distribution LHS sample
  parID   = c("cW","cG"   ,"cQ","cV","cS" ,"cD" )[fit_T]
  parmean = c(cW  , cG/1e6, cQ , cV , cS  , cD  )[fit_T]
  parsd   = c(100 , 20    , 20 , 10 , 1   , 250 )[fit_T]
  parmin  = c(1   , 0.1   , 1  , 1  , 3   , 1000)[fit_T]
  parmax  = c(450 , 100   , 150, 50 , 6   , 2500)[fit_T]
  
  # make sample
  LHS = TruncGaussianLHS(samplesize, parmean, parsd, parmin, parmax, varnames=parID)
  if(fit_cW==F){LHS$cW = cW}
  if(fit_cG==F){LHS$cG = cG}
  if(fit_cQ==F){LHS$cQ = cQ}
  if(fit_cV==F){LHS$cV = cV}
  if(fit_cS==F){LHS$cS = cS}
  if(fit_cD==F){LHS$cD = cD}
  # reorder columns
  LHS=LHS[c("cW", "cG", "cQ", "cV", "cS", "cD")] 
  
  
  ############
  ### LHS runs
  ############
  
  # Make an empty vector to store mean sums of squares during for-loop.
  LHS$SS = c()  # expand: other objective functions, then select which ones
  
  # Run a for-loop over all parameter sets and run WALRUS in every iteration.
  for(i in 1:samplesize)
  {
    print(paste("Latin Hypercube Sample", i))
    # look up the parameter set for this iteration
    parameters = data.frame(cW=LHS$cW[i], cG=LHS$cG[i]*1e6, cQ=LHS$cQ[i], 
                            cV=LHS$cV[i], cS=LHS$cS[i], cD=LHS$cD[i], 
                            aS=aS, st=st, Gfrac=Gfrac)
    
    # run WALRUS
    modeled    = WALRUS_loop(pars=parameters)

    #print(length(Qobs_forNS))
    #print(length(modeled$Q))
    
    # Compute and store the mean sum of squares.
    LHS$SS[i]  = mean((Qobs_forNS-modeled$Q)^2)
  }
  
  # Write the results to file: parameter values and belonging sum of squares.
  write.table(LHS, paste0("LHS_pars_", name, ".dat"))
  
  # Select best (N_LM) runs
  LHS_best = LHS[order(LHS$SS)[1:selectionsize],]
  
  
  ########################
  # Prepare LM-calibration
  ########################
  
  # Rewrite the model such that the output is a vector of residuals: Qobs-Qmod.
  # Division by parmean is to normalize parameters (better for optimization technique).
  WALRUS_for_optim = function(par)
  {
    nfits=0
    if(fit_cW==F){cW_optim=cW}else{
      nfits=nfits+1
      cW_optim=par[nfits]*parmean[nfits]}
    if(fit_cG==F){cG_optim=cG}else{
      nfits=nfits+1
      cG_optim=par[nfits]*parmean[nfits]*1e6}
    if(fit_cQ==F){cQ_optim=cQ}else{
      nfits=nfits+1
      cQ_optim=par[nfits]*parmean[nfits]}
    if(fit_cV==F){cV_optim=cV}else{
      nfits=nfits+1
      cV_optim=par[nfits]*parmean[nfits]}
    if(fit_cS==F){cS_optim=cS}else{
      nfits=nfits+1
      cS_optim=par[nfits]*parmean[nfits]}
    if(fit_cD==F){cD_optim=cD}else{
      nfits=nfits+1
      cD_optim=par[nfits]*parmean[nfits]}
    
    fit_pars = data.frame(cW=cW_optim, cG=cG_optim, cQ=cQ_optim,
                          cV=cV_optim, cS=cS_optim, cD=cD_optim, 
                          aS=aS, st=st, Gfrac=Gfrac)    
    mod = WALRUS_loop(pars=fit_pars)
    return(Qobs_forNS - mod$Q)
  }
  
  
  #############
  # Calibration
  #############
  
  LM_pars = matrix(nrow=selectionsize, ncol=length(fit_T))
  
  for(i in 1:selectionsize)
  {
    print(paste("Levenberg-Marquardt", i))
    # Define initial and boundary values for the parameters you want to calibrate.
    # You can also calibrate the initial values (if you want).
    par_start = LHS_best[i,fit_T]
    print(par_start)
    
    # This is the Levenberg-Marquardt optimization algorithm
    # which will minimize the sum of squares of the residuals. 
    # Increasing nprint will give more information.
    # Maxiter is the maximum number of iterations.
    cal = nls.lm(par=par_start/parmean, lower=parmin/parmean, upper=parmax/parmean, 
          fn=WALRUS_for_optim, 
          control=nls.lm.control(nprint=1,maxiter=100))
    
    # Retrieve the optimal parameter values found in the calibration procedure.
    LM_pars[i,] = coef(cal) * parmean
  }
  
  
  #########################
  # Run with optimal values
  #########################
  
  for(i in 1:selectionsize)
  {
    calname = paste0("cal_",name,"_", i)
    print(name)
    nfits=0
    if(fit_cW==F){cW_LM=cW}else{
      nfits=nfits+1
      cW_LM=LM_pars[i,nfits]}
    if(fit_cG==F){cG_LM=cG}else{
      nfits=nfits+1
      cG_LM=LM_pars[i,nfits]*1e6}
    if(fit_cQ==F){cQ_LM=cQ}else{
      nfits=nfits+1
      cQ_LM=LM_pars[i,nfits]}
    if(fit_cV==F){cV_LM=cV}else{
      nfits=nfits+1
      cV_LM=LM_pars[i,nfits]}
    if(fit_cS==F){cS_LM=cS}else{
      nfits=nfits+1
      cS_LM=LM_pars[i,nfits]}
    if(fit_cD==F){cD_LM=cD}else{
      nfits=nfits+1
      cD_LM=LM_pars[i,nfits]}
    pars = data.frame(cW=cW_LM, cG=cG_LM, cQ=cQ_LM,
                      cV=cV_LM, cS=cS_LM, cD=cD_LM,
                      aS=aS, st=st, Gfrac=Gfrac)
    mod  = WALRUS_loop(pars=pars)
    WALRUS_postprocessing(o=mod, pars=pars, n=calname)
  }

}








