#' @title Bayesian Inference for a Hierarchical Linear Mixture Model with a Bernoulli Classifier
#'
#' @description
#' Using [rstan] the function fits a hierarchical linear mixture model to a
#' dataset. The model was developed to be able to classify individual fibres from a
#' patient dataset as being healthy or deficient based on a dataset of healthy
#' control fibres. It is assumed that healthy data shows strong
#' linearity and anything that diverges from this is deficient. Individual fibres are
#' classified with a Bernoulli classifier whose parameter is the proportion of
#' deficient fibres within the sample. The parameters and hyper-parameters of
#' prior distributions can be changed to account for different prior beliefs. See
#' details for more a list of parameters whose value can be changed.
#'
#' @details
#' The inference is implemented using [rstan::stan]. The model itself
#' cannot be changed however the parameters and hyper parameters which govern
#' the model can be updated by using the `parameterVals` argument in the
#' function. This allows for a vast array of different prior beliefs and
#' information to be imposed upon the model, which may change the resulting
#' classification and posterior parameter beliefs.
#'
#' Beliefs about the slope are summarised using a normal distribution, however
#' the parameters of this normal distribution, mu_m and tau_m, are unknown and
#' so we impart beliefs about the slope by our choice of priors on mu_m and
#' tau_m. The prior densities for the two parameters are normal and gamma
#' distributions respectively. In the model the hyper-parameters which define
#' these distributions are:
#'  - mean_mu_m : the expected value of mu_m, default = 1.0
#'  - prec_mu_m : the precision of mu_m, default = 100
#'  - shape_mu_m : the shape parameter of the tau_m distribution, default = 27.56372
#'  - rate_mu_m : the rate parameter of the tau_m gamma distribution, default = 1.660233
#'
#' Similarly, the intercept is summarised by a normal distribution whose
#' parameters are not known and follow normal and gamma distributions a priori.
#'  - mean_mu_c : the expected value of mu_c, default = 0.0
#'  - prec_mu_c : the precision of mu_c, default = 25
#'  - shape_mu_c : the shape parameter of the tau_c distribution, default = 1.25
#'  - rate_mu_c : the rate parameter of the tau_c gamma distribution, default = 0.25
#'
#' The model error for healthy patients is considered unknown and follows a
#' gamma distribution a priori.
#'  - shape_tau : the shape parameter for tau, default = 41.97618
#'  - rate_tau : the rate parameter for tau, default = 2.048809
#'
#' The error for the second component, which classifies deficient fibres, is
#' fixed.
#'  - tau_def : the precision of the second component, default = 1e-4
#'
#' The proportion of deficiency follows a beta distribution, with two shape
#' parameters.
#'  - alpha_pi : the first shape parameter, default = 1.0
#'  - beta_pi : the second shape parameter, default = 1.0
#'
#' Any subset of these parameters can be changed for inference by storing their
#' value in a list and passing it to `parameterVals`. If nothing is passed to
#' `parameterVals` the model uses the default parameters, which may in
#' inappropriate for your data set.
#'
#' @param dataMats A list consisting of four elements. The list should be named
#' or the elements should be in the following order.
#'     - ctrl : a matrix of the control subject data.
#'     - pts : a matrix of the patient subject data.
#'     - indexCtrl : a numeric vector indicating which rows of the control subject matrix belong to the same subject.
#'     - indexPat : a numeric vector indicating which rows of the patient subject matrix belong to the same subject.
#' @param parameterVals A list of named values to replace the parameter values used to define prior distributions.
#' @param MCMCout The final number of posterior draws (disregarding burn-in and thinning). Default is 1000.
#' @param MCMCburnin The number of posterior draws discarded for burn-in. Default is 1000.
#' @param MCMCthin The lag at which to the thin the posterior draws. Default is 1 i.e. no thinning.
#'
#' @returns A list of different components of the MCMC output:
#'    - post : a data frame of the posterior draws for all variables, where each column is a different variable
#'    - prior : a data frame of draws from the prior distribution for all variables.
#'    - postpred : a data frame of posterior predictions for the patient subject linear model, with a column for x values and columns for 95% predictive interval and expected values for both the like control and not control models.
#'    - priorpred : a data frame of prior predictions for the patient subject model, similar to postpred.
#'    - classif : a data frame of whose columns correspond to patient fibres and rows are classifications from individual draws from the posterior distirbution
#'
#' @examples
#' exampleData = get_exampleData()
#' #' # the measure of mitochondrial mass - the x-axis of the 2D mito plot
#' mitochan = "VDAC"
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#' sbj = unique(exampleData$sampleID)
#' ctrlID = grep("C", sbj, value=TRUE)
#' pts = grep("C", sbj, value=TRUE, invert=TRUE)
#'
#' chan = channels[1]
#' pat = pts[1]
#'
#' data_mat = getData_mats(exampleData, channels=c(mitochan, chan), ctrlID=ctrlID, pts=pat)
#'
#' tau_mode = 100
#' tau_var = 10
#' rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
#' shape_tau = 1 + tau_mode * rate_tau
#'
#' infOut = stan_inference(data_mat, parameterVals=list(shape_tau=shape_tau, rate_tau=rate_tau))
#'
#'
#' @importFrom rstan stan
#'
#' @export

stan_inference = function(dataMats,
                          MCMCout=1000, MCMCburnin=1000, MCMCthin=1,
                          nChains=10, nCores=1,
                          parameterVals = NULL,
                          save=FALSE,
                          saveOnly=FALSE,
                          saveRoot=""){

  nCtrl = nrow( dataMats$ctrl )
  ctrl_mat = dataMats$ctrl
  ctrl_index = dataMats$indexCtrl
  nPat = nrow( dataMats$pts )
  nCrl = length( unique( dataMats$indexCtrl) )

  nSyn = 1e3
  xSyn = seq(min(c(dataMats$ctrl[,1], dataMats$pts[,1]))-2,
             max(c(dataMats$ctrl[,1], dataMats$pts[,1]))+2,
             length.out=nSyn)

  grad = vector("numeric", length=nCrl)
  inter = grad
  prec = grad

  for( i in 1:nCrl ){
    xCtrl = ctrl_mat[ctrl_index==i, 1]
    yCtrl = ctrl_mat[ctrl_index==i, 2]
    dd = data.frame(mitochan=xCtrl, chan=yCtrl)

    mod = lm(chan ~ mitochan, data=dd)
    grad[i] = mod$coefficients[2]
    inter[i] = mod$coefficients[1]
    prec[i] = 1 / summary(mod)$sigma^2
  }

  tau_mode = mean( prec )
  tau_var = 10
  rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
  shape_tau = 1 + tau_mode*rate_tau

  tau_m_mode = 50
  tau_m_var = 50
  rate_tau_m = 0.5*(tau_m_mode + sqrt(tau_m_mode^2+4*tau_m_var)) / tau_m_var
  shape_tau_m = 1 + tau_m_mode*rate_tau_m

  tau_c_mode = 50
  tau_c_var = 50
  rate_tau_c = 0.5*(tau_c_mode + sqrt(tau_c_mode^2+4*tau_c_var)) / tau_c_var
  shape_tau_c = 1 + tau_c_mode*rate_tau_c

  param_list = list(
    D=2,
    K=2,
    M=nCrl+1,
    ctrl_mat = dataMats$ctrl,
    pat_mat = dataMats$pts,
    nCtrl = nCtrl,
    nCrl = nCrl,
    nPat = nPat,
    nPts = 1,
    ctrlIndex = dataMats$indexCtrl,
    nSyn = nSyn,
    xSyn = xSyn,
    mean_mu_m = mean( grad ), prec_mu_m = 1/0.25^2,
    mean_mu_c = mean( inter ), prec_mu_c = 1/0.25^2,
    shape_tau_m = shape_tau_m, rate_tau_m = rate_tau_m,
    shape_tau_c = shape_tau_c, rate_tau_c = rate_tau_c,
    shape_tau = shape_tau, rate_tau = rate_tau,
    alpha_pi = 1.0, beta_pi = 1.0,
    slope_lb = 0.1,
    tau_def=1e-2
  )

  if (!is.null(parameterVals) && is.list(parameterVals)) {
    for (param in names(parameterVals) ) {
      if (param %in% names(param_list) ) {
        param_list[[param]] = parameterVals[[param]]
      } else {
        message(paste("The parameter `", param, "` is not part of the model."))
      }
    }
  }

  init_list = list()
  for(i in 1:nChains){
    gen_list = list("mu_m"=rnorm(1, param_list[["mean_mu_m"]], 1/sqrt(param_list[["prec_mu_m"]]) ),
                    "mu_c"=rnorm(1, param_list[["mean_mu_c"]], 1/sqrt(param_list[["prec_mu_c"]])),
                    "tau_m"= rgamma(1, param_list[["shape_tau_m"]], param_list[["rate_tau_m"]]),
                    "tau_c"=rgamma(1, param_list[["shape_tau_c"]], param_list[["rate_tau_c"]]),
                    "tau_norm"=rgamma(1, param_list[["shape_tau"]], param_list[["rate_tau"]]),
                    "probdiff"=rbeta(1,param_list[["alpha_pi"]],  param_list[["beta_pi"]] ) )
    m_init = list()
    c_init = list()
    for(j in 1:(nCrl+1)){
      m_init[[paste0("m[",j,"]")]] = rnorm(1, gen_list[["mu_m"]], 1/gen_list[["tau_m"]])
      c_init[[paste0("c[",j,"]")]] = rnorm(1, gen_list[["mu_c"]], 1/gen_list[["tau_c"]])
    }
    init_list[[i]] = c(gen_list, m_init, c_init)
  }

  output = stan("../inst/stan/bhlmm_stan.stan", data=param_list, chains=nChains,
                init=init_list, cores=nCores,
                iter=(MCMCout+MCMCburnin)*MCMCthin, warmup=MCMCburnin, thin=MCMCthin)

  outmat = as.matrix(output)
  outcols = colnames(outmat)

  classifs_mat = outmat[, grepl("classif", outcols)]
  post = outmat[, !( grepl("classif", outcols)|grepl("lp__", outcols)|
                       grepl("probvec", outcols)|grepl("dens", outcols)|
                       grepl("yPred", outcols)|grepl("_tmp", outcols)|
                       grepl("_prior", outcols) ) ]
  postpred_mat = outmat[, grepl("yPred", outcols) & !grepl("_prior", outcols)]

  prior = outmat[, !( grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probvec", outcols)|grepl("dens", outcols)|
                        grepl("yPred", outcols)|grepl("_tmp", outcols) ) & grepl("_prior", outcols) ]
  colnames(prior) = gsub("_prior", "", colnames(prior))

  priorpred_mat = outmat[, grepl("yPred", outcols) & grepl("_prior", outcols)]

  if( nChains == 1 ){
    postpred = apply(postpred_mat, 2, quantile, probs=c(0.025, 0.5, 0.975))
    postpred = cbind(xSyn, t(postpred))
    colnames(postpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")

    priorpred = apply(priorpred_mat, 2, quantile, probs=c(0.025, 0.5, 0.975))
    priorpred = cbind(xSyn, t(priorpred))
    colnames(priorpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
  } else {
    postpred = NULL
    priorpred = NULL

    for( i in 1:nChains ){
      chain_ind = ((i - 1)*MCMCout + 1):(i*MCMCout)
      postpred_chain = apply(postpred_mat[chain_ind, ], 2, quantile, probs=c(0.025, 0.5, 0.975))
      priorpred_chain = apply(priorpred_mat[chain_ind, ], 2, quantile, probs=c(0.025, 0.5, 0.975))
      postpred = rbind(postpred, t(postpred_chain))
      priorpred = rbind(priorpred, t(priorpred_chain))
    }
    postpred = cbind(rep(xSyn, nChains), postpred)
    priorpred = cbind(rep(xSyn, nChains), priorpred)
    colnames(postpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
    colnames(priorpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
  }

  outList = list(POST=post, POSTPRED=postpred, CLASSIF=classifs_mat,
                 PRIOR=prior, PRIORPRED=priorpred)
  if( save ){
    list_saver(outList, root = saveRoot)
    if( saveOnly ){ return( NULL ) }
  }
  return( outList )
}





