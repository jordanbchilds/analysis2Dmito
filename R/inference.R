#' @title Bayesian Inference for a Hierarchical Linear Mixture Model with a Bernoulli Classifier
#'
#' @description
#' Using [rjags] the function fits a hierarchical linear mixture model to a
#' dataset. The was developed to be able to classify individual fibres within a
#' 2Dmito plot as being a set of data control fibre data points. The model
#' fundamentally assumes the control data shows strong linearity and anything
#' diverges from this is not like control.
#'
#' @details
#' The inference is implemented using [rjags::jags.model()]. The model itself
#' cannot be changed however the parameters and hyper parameters which govern
#' the model can be updated by using the `parameterVals` variable in the
#' function. This allows for a vast array of different prior beliefs and
#' information to be imposed upon the model.
#'
#' The function fits a Bayesian hierarchical linear mixture model, classifying
#' the patient fibres as like control or not using a Bernoulli classifier.
#' Let Y be the protein expression level of interest and X be a measure of
#' mitochondrial mass. Suppose that Y\[i,j\] is the j-th fibre from sample i. This
#'notation is purely for convenience here and is not syntax used in the function.
#' Here normal distributions are described using its mean and precision,
#' following the convention of [rjags].
#'
#' The **model**:
#'  Y\[i,j\] ~ Normal( m\[i\]*X\[i,j\]+c\[i\], tau\[i,j\] )
#'
#' Model error: tau\[i,j\]
#'  - tau_def if Z\[i,j\] == 1
#'  - Gamma( shape_tau, rate_tau ) if Z\[i,j\] == 0
#'
#' **Classification**: Z\[i,j\]
#'   - 0 if i is control
#'   - Bernoulli( probdef ) if i is patient
#'
#' **Probability of deficiency**: probdef ~ log-Normal(mu_p, tau_p) T(,1)
#' This is a right-truncated log-normal distribution constraining the parameter
#' space to be \[0,1\]. A truncated log-normal was used instead of Beta
#' distribution (or another density with \[0,1\] support) as this could more
#' accurately captured the show shape of the published data ENTER REFERENCE.
#' The **priors** of the slope and intercept are:
#'  - m\[i\]~ Normal( mu_m0, tau_m0 )
#'  - c\[i\]~ Normal( mu_c0, tau_c0 )
#'
#' The **hyper-priors** for the slop and intercept are:
#'. - mu_m0 ~ Normal( mean_mu_m0, prec_mu_m0 )
#'  - mu_c0 ~ Normal( mean_mu_c0, prec_mu_c0 )
#'  - tau_m0 ~ Gamma( shape_tau_m0, rate_tau_m0 )
#'  - tau_c0 ~ Gamma( shape_tau_c0, rate_tau_c0 )
#'
#' The parameters and hyper-parameters which control the model and can be changed using the `parameterValue` argument of the function are:
#'  - mean_mu_m0 : the expected value of a slope. Default: 1.0.
#'  - prec_mu_m0 : the precision of the expected value of a slope. Default: 16.
#'  - mean_mu_c0 : the expected value of an intercept. Default: 0.0.
#'  - prec_mu_c0 : the precision of the expected value of an intercept. Default: 4/9
#'  - shape_tau_m0 : the shape parameter of a Gamma distribution of the precision of a slope. Default: 1.1956.
#'  - rate_tau_m0 : the rate parameter of a Gamma distribution of the precision of a slope. Default: 0.04889989.
#'  - shape_tau_c0 : the shape parameter of a Gamma distribution of the precision of an intercept. Default: 1.1956.
#'  - rate_tau_c0 : the rate parameter of a Gamma distribution of the precision of an intercept. Default: 0.04889989.
#'  - shape_tau : shape parameter of a Gamma distribution of the precision in the model for healthy, like-control, fibres. Defualt: 41.048809.
#'  - rate_tau : the rate parameter of a Gamma distribution of the precision in the model for healthy, like-control, fibres. Default: 2.048809.
#'  - mu_p : the expected value for the proportion of deficiency in the dataset. Default: -2.549677.
#'  - tau_p : the precision for the proportion of deficiency in the dataset. Default: 0.9549515.
#'  - tau_def : the precision used to classif fibres which are not like control. Default: 0.0001.
#'
#' The model is implemented in JAGS using the following model string:
#' ```{R}
#' "
#' model {
#'  for(i in 1:nCtrl){
#'     Yctrl[i] ~ dnorm(m[index[i]]*Xctrl[i]+c[index[i]], tau_norm)
#'  }
#'  for(j in 1:nPat){
#'    Ypat[j] ~ dnorm(m[nCrl+1]*Xpat[j]+c[nCrl+1], tau_hat[j])
#'    tau_hat[j] = ifelse(class[j]==1, tau_def, tau_norm)
#'    class[j] ~ dbern(probdef)
#'  }
#'  for(k in 1:Nsyn){
#'    Ysyn_norm[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_norm)
#'    Ysyn_def[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_def)
#'  }
#'  tau_norm ~ dgamma(shape_tau, rate_tau)
#'  probdef ~ dlnorm(mu_p, tau_p) T(,1)
#'  for(i in 1:(nCrl+1)){
#'    m[i] ~ dnorm(mu_m0, tau_m0)
#'    c[i] ~ dnorm(mu_c0, tau_c0)
#'  }
#'
#'  # Parent distribution for the slope and intercept
#'  m_pred ~ dnorm(mu_m0, tau_m0)
#'  c_pred ~ dnorm(mu_c0, tau_m0)
#'
#'  # Hyper-priors
#'  mu_m0 ~ dnorm( mean_mu_m0, prec_mu_m0 )
#'  mu_c0 ~ dnorm( mean_mu_c0, prec_mu_c0 )
#'  tau_m0 ~ dgamma( shape_tau_m0, rate_tau_m0 )
#'  tau_c0 ~ dgamma( shape_tau_c0, rate_tau_c0 )
#' }
#' "
#' ```
#' @param dataMats A list consisting of four elements. The list should be named
#' or the elements should be in the following order.
#'     - ctrl : a matrix of the control subject data.
#'     - pts : a matrix of the patient subject data.
#'     - indexCtrl : a numeric vector indicating which rows of the control subject matrix belong to the same subject.
#'     - indexPat : a numeric vector indicating which rows of the patient subject matrix belong to the same subject.
#' @param parameterVals A list or named vector whose names correspond to the parameters of the model and the whose value is the value of the model (hyper-)parameter.
#' @param MCMCout The final number of posterior draws (disregarding burn-in and thinning).
#' @param MCMCburnin The number of posterior draws discarded for burn-in.
#' @param MCMCthin The lag at which to the thin the posterior draws.
#'
#' @returns A list of different components of the MCMC output:
#'    - post : a dataframe of the posterior draws for all variables, where each column is a different variable
#'    - prior : a data frame of draws from the prior distribution for all variables.
#'    - postpred : a dataframe of posterior predictions for the patient subject linear model, with a column for x values and columns for 95% predictive interval and expected values for both the like control and not control models.
#'    - priorpred : a dataframe of prior predictions for the patient subject model, similar to postpred.
#'    - classif : a dataframe of whose columns correspond to patient fibres and rows are classifications from individual draws from the posterior distirbution
#'
#' @examples
#' exampleData = get_exampleData()
#' #' # the measure of mitochondrial mass - the x-axis of the 2D mito plot
#' mitochan = "raw_porin"
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#' sbj = unique(exampleData$sampleID)
#' ctrlid = c("C01", "C02", "C03", "C04", "C05")
#' pts = sbj[ !(sbj %in% ctrlid) ]
#' chan = channels[1]
#' pat = pts[1]
#'
#' data_mat = getData_mats(exampleData, cord=c(mitochan, chan), ctrlID=ctrlid, pts=pat, getIndex=TRUE)
#'
#' tau_mode = 100
#' tau_var = 500
#' rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
#' shape_tau = 1 + tau_mode * rate_tau
#'
#' infOut = inference(data_mat, parameterVals=list(shape_tau=shape_tau, rate_tau=rate_tau))
#'
#'
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#'
#' @export
inference = function(dataMats,
                     parameterVals=NULL,
                     MCMCout = 1000,
                     MCMCburnin = 1000,
                     MCMCthin = 1 ) {
  modelstring = "
       model {
        for(i in 1:nCtrl){
          Yctrl[i] ~ dnorm(m[indexCtrl[i]]*Xctrl[i]+c[indexCtrl[i]], tau_norm)
        }
        for(j in 1:nPat){
            Ypat[j] ~ dnorm(m[nCrl+1]*Xpat[j]+c[nCrl+1], tau_hat[j])
            tau_hat[j] = ifelse(class[j]==1, tau_def, tau_norm)
            class[j] ~ dbern(probdiff)
        }
        for(k in 1:nSyn){
            Ysyn_norm[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_norm)
            Ysyn_def[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_def)
        }
        for(l in 1:(nCrl+1)){
          m[l] ~ dnorm(mu_m, tau_m)
          c[l] ~ dnorm(mu_c, tau_c)
        }
        m_pred ~ dnorm(mu_m, tau_m)
        c_pred ~ dnorm(mu_c, tau_m)

        mu_m ~ dnorm( mean_mu_m, prec_mu_m )
        mu_c ~ dnorm( mean_mu_c, prec_mu_c )
        tau_m ~ dgamma( shape_tau_m, rate_tau_m )
        tau_c ~ dgamma( shape_tau_c, rate_tau_c )
        tau_norm ~ dgamma(shape_tau, rate_tau)
        probdiff ~ dbeta(alpha_pi, beta_pi)
       }
      "
  if (is.null(names(dataMats))) names(dataMats) = c("ctrl", "pts", "indexCtrl", "indexPat")
  ctrl_mat = dataMats$ctrl
  pat_mat = dataMats$pts
  nCtrl = nrow(ctrl_mat)
  nPat = nrow(pat_mat)
  indexCtrl = dataMats$indexCtrl
  indexPat = dataMats$indexPat
  nCrl = length(unique(indexCtrl))

  nSyn = 1e3
  Xsyn = seq(min(0, min(c(ctrl_mat[,1], pat_mat[,1]))),
             max(c(ctrl_mat[, 1], pat_mat[, 1])) * 1.5,
             length.out = nSyn)

  data_list = list(
    Yctrl = ctrl_mat[, 2],
    Xctrl = ctrl_mat[, 1],
    Ypat = pat_mat[, 2],
    Xpat = pat_mat[, 1],
    nCtrl = nCtrl,
    nCrl = nCrl,
    nPat = nPat,
    indexCtrl = indexCtrl,
    nSyn = nSyn,
    Xsyn = Xsyn,
    mean_mu_m = 1.0,
    prec_mu_m = 1/0.1^2, #1 / 0.25 ^ 2,
    mean_mu_c = 0.0,
    prec_mu_c = 1/0.2^2, #1 / 1.5 ^ 2,
    shape_tau_m = 27.56372, #1.1956,
    rate_tau_m = 1.660233, #0.04889989,
    shape_tau_c = 1.25, #1.1956,
    rate_tau_c = 0.25, #0.04889989,
    shape_tau = 41.97618,
    rate_tau = 2.048809,
    alpha_pi = 1,
    beta_pi = 1,
    tau_def=1e-4
  )

  if (!is.null(parameterVals) && is.list(parameterVals)) {
    for (param in names(parameterVals)) {
      if (param %in% names(data_list)) {
        data_list[[param]] = parameterVals[[param]]
      } else {
        message(paste("The parameter `", param, "` is not part of the model."))
      }
    }
  }

  data_priorpred = data_list
  data_priorpred$Yctrl = NULL
  data_priorpred$nCtrl = 0
  data_priorpred$Ypat = NULL
  data_priorpred$nPat = 0

  model_pat = rjags::jags.model(textConnection(modelstring),
                         data = data_list,
                         n.chains = 1)
  model_pat_priorpred = rjags::jags.model(textConnection(modelstring), data =
                                     data_priorpred)

  output_post = rjags::coda.samples(
    model = model_pat,
    n.iter = MCMCout * MCMCthin,
    thin = MCMCthin,
    variable.names = c(
      "m",
      "c",
      "m_pred",
      "c_pred",
      "mu_m",
      "tau_m",
      "mu_c",
      "tau_c",
      "tau_norm",
      "Ysyn_norm",
      "Ysyn_def",
      "class",
      "probdiff"
    )
  )
  output_prior = rjags::coda.samples(
    model = model_pat_priorpred,
    n.iter = MCMCout,
    thin = 1,
    variable.names = c(
      "m",
      "c",
      "m_pred",
      "c_pred",
      'mu_m',
      "tau_m",
      "mu_c",
      "tau_c",
      "tau_norm",
      "Ysyn_norm",
      "Ysyn_def",
      "probdiff"
    )
  )

  post_out = as.data.frame(output_post[[1]])
  prior_out = as.data.frame(output_prior[[1]])

  summ_pat = summary(output_post)
  classifs = post_out[, grepl("class", colnames(post_out))]
  post_names = colnames(post_out)

  post = post_out[, c(
    paste0("m[", 1:(nCrl + 1), "]"),
    paste0("c[", 1:(nCrl + 1), "]"),
    "m_pred",
    "c_pred",
    "mu_m",
    "tau_m",
    "mu_c",
    "tau_c",
    "tau_norm",
    "probdiff"
  )]
  postpred_norm = apply(post_out[, paste0("Ysyn_norm[", 1:nSyn, "]")], 2, quantile, probs =
                          c(0.025, 0.5, 0.975))
  postpred_def = apply(post_out[, paste0("Ysyn_def[", 1:nSyn, "]")], 2, quantile, probs =
                         c(0.025, 0.5, 0.975))
  postpred = cbind(Xsyn, t(postpred_norm), t(postpred_def))
  colnames(postpred) = c("mitochan",
                         "lwrNorm",
                         "medNorm",
                         "uprNorm",
                         "lwrDef",
                         "medDef",
                         "uprDef")

  prior_names = colnames(prior_out)
  prior = prior_out[, c(
    paste0("m[", 1:(nCrl + 1), "]"),
    paste0("c[", 1:(nCrl + 1), "]"),
    "m_pred",
    "c_pred",
    "mu_m",
    "tau_m",
    "mu_c",
    "tau_c",
    "tau_norm",
    "probdiff"
  )]

  priorpred_norm = apply(prior_out[, paste0("Ysyn_norm[", 1:nSyn, "]")], 2, quantile, probs =
                           c(0.025, 0.5, 0.975))
  priorpred_def = apply(prior_out[, paste0("Ysyn_def[", 1:nSyn, "]")], 2, quantile, probs =
                          c(0.025, 0.5, 0.975))
  priorpred = cbind(Xsyn, t(priorpred_norm), t(priorpred_def))
  colnames(priorpred) = c("mitochan",
                          "lwrNorm",
                          "medNorm",
                          "uprNorm",
                          "lwrDef",
                          "medDef",
                          "uprDef")

  out_list = list(
    POST = post,
    POSTPRED = postpred,
    PRIOR = prior,
    PRIORPRED = priorpred,
    CLASSIF = classifs
  )
  return(out_list)
}
