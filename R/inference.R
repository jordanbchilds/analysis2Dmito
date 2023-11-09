#' @title Bayesian Inference for a Hierarchical Linear Mixture Model with a Bernoulli Classifier
#'
#' @description
#' Using [rjags] the function fits a hierarchical linear mixture model to a
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
#' The inference is implemented using [rjags::jags.model]. The model itself
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
        c_pred ~ dnorm(mu_c, tau_c)

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
    prec_mu_m = 1/0.1^2,
    mean_mu_c = 0.0,
    prec_mu_c = 1/0.2^2,
    shape_tau_m = 27.56372,
    rate_tau_m = 1.660233,
    shape_tau_c = 1.25,
    rate_tau_c = 0.25,
    shape_tau = 41.97618,
    rate_tau = 2.048809,
    alpha_pi = 1,
    beta_pi = 1,
    tau_def = 1e-4
  )

  if (!is.null(parameterVals) && is.list(parameterVals)) {
    for (param in names(parameterVals)) {
      if (param %in% names(data_list)) {
        data_list[[param]] = parameterVals[[param]]
      } else {
        message(paste("The parameter `", param, "` does not exist. See function documentation for a list of possible parameters."))
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
