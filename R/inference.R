#' @title Bayesian Inference for a Hierarchical Linear Mixture Model with a Bernoulli Classifier
#'
#' @description
#' Using [rjags] the function fits a hierarchical linear mixture model to a
#' dataset. The was developed to be able to classify individual fibres within a
#' 2Dmito plot as being a set of data control fibre datapoints. The model
#' fundamentally assumes the contrl data shows strong linearity and anything
#' diverges from this is not like control.
#'

inference = function(dataMats,
                     parameterVals = NULL,
                     MCMCout = 1000,
                     MCMCburnin = 1000,
                     MCMCthin = 1,
                     MCMCchain = 1) {
  modelstring = "
       model {
       # Fit the model to the control subjects
        for(i in 1:nCtrl){
          Yctrl[i] ~ dnorm(m[indexCtrl[i]]*Xctrl[i]+c[indexCtrl[i]], tau_norm)
        }

        for(j in 1:nPat){
            Ypat[j] ~ dnorm(m[nCrl+1]*Xpat[j]+c[nCrl+1], tau_hat[j])
            tau_hat[j] = ifelse(class[j]==1, tau_def, tau_norm)
            class[j] ~ dbern(probdef)
        }

        # Predictive for the patient subject
        for(k in 1:nSyn){
            Ysyn_norm[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_norm)
            Ysyn_def[k] ~ dnorm(m[nCrl+1]*Xsyn[k]+c[nCrl+1], tau_def)
        }

        # Prior beliefs about parameters
        for(l in 1:(nCrl+nPt)){
          m[l] ~ dnorm(mu_m0, tau_m0)
          c[l] ~ dnorm(mu_c0, tau_c0)
        }

        m_pred ~ dnorm(mu_m0, tau_m0)
        c_pred ~ dnorm(mu_c0, tau_m0)

        mu_m0 ~ dnorm( mean_mu_m0, prec_mu_m0 )
        mu_c0 ~ dnorm( mean_mu_c0, prec_mu_c0 )
        tau_m0 ~ dgamma( shape_tau_m0, rate_tau_m0 )
        tau_c0 ~ dgamma( shape_tau_c0, rate_tau_c0 )
        tau_norm ~ dgamma(shape_tau, rate_tau)
        # tau_def ~ dgamma(shape_tauDef, rate_tauDef)
        probdef ~ dlnorm(mu_p, tau_p)
       }
      "

  ctrl_mat = dataMats$ctrl
  pat_mat = dataMats$pts
  nCtrl = nrow(ctrl_mat)
  nPat = nrow(pat_mat)
  indexCtrl = dataMats$indexCtrl
  indexPat = dataMats$indexPat
  nPt = length(unique(indexPat))
  nCrl = length(unique(indexCtrl))


  # prior parameters for control data
  mean_mu_m0 = 1.0
  prec_mu_m0 = 1 / 0.25 ^ 2
  mean_mu_c0 = 0.0
  prec_mu_c0 = 1 / 1.5 ^ 2

  tau_m0_mode = 1 / 0.5 ^ 2 # expected sd of 0.5
  tau_m0_var = 500
  rate_tau_m0 = 0.5 * (tau_m0_mode + sqrt(tau_m0_mode ^ 2 + 4 * tau_m0_var)) / tau_m0_var
  shape_tau_m0 = 1 + tau_m0_mode * rate_tau_m0

  tau_c0_mode = 1 / 0.5 ^ 2 # expected sd of 0.5
  tau_c0_var = 500
  rate_tau_c0 = 0.5 * (tau_c0_mode + sqrt(tau_c0_mode ^ 2 + 4 * tau_c0_var)) / tau_c0_var
  shape_tau_c0 = 1 + tau_c0_mode * rate_tau_c0

  tau_mode = 1 / sqrt(0.05) ^ 2
  tau_var = 10
  rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
  shape_tau = 1 + tau_mode * rate_tau

  mu_p = -2.549677 # values taken from the Ahmed_2022 data
  tau_p = 1 / 1.023315 ^ 2 # values taken from the Ahmed_2022 data

  tau_def = 0.1
  tauDef_mode = 1 / 15 ^ 2 # expected sd of 5
  tauDef_var = 1
  rate_tauDef = 0.5 * (tauDef_mode + sqrt(tauDef_mode ^ 2 + 4 * tauDef_var)) / tauDef_var
  shape_tauDef = 1 + tauDef_mode * rate_tauDef

  # curve(dgamma(x, shape_tauDef, rate_tauDef), to=10)

  nSyn = 1e3
  Xsyn = seq(0, max(c(ctrl_mat[, 1], pat_mat[, 1])) * 1.5, length.out =
               nSyn)
  tau_def = 0.001

  data_list = list(
    Yctrl = ctrl_mat[, 2],
    Xctrl = ctrl_mat[, 1],
    Ypat = pat_mat[, 2],
    Xpat = pat_mat[, 1],
    nCtrl = nCtrl,
    nCrl = nCrl,
    nPt = nPt,
    nPat = nPat,
    indexCtrl = indexCtrl,
    nSyn = nSyn,
    Xsyn = Xsyn,
    mean_mu_m0 = mean_mu_m0,
    prec_mu_m0 = prec_mu_m0,
    mean_mu_c0 = mean_mu_c0,
    prec_mu_c0 = prec_mu_c0,
    shape_tau_m0 = shape_tau_m0,
    rate_tau_m0 = rate_tau_m0,
    shape_tau_c0 = shape_tau_c0,
    rate_tau_c0 = rate_tau_c0,
    shape_tau = shape_tau,
    rate_tau = rate_tau,
    mu_p = mu_p,
    tau_p = tau_p,
    # shape_tauDef = shape_tauDef,
    # rate_tauDef = rate_tauDef
    tau_def = tau_def
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
                         n.chains = MCMCchain)
  model_pat_priorpred = rjags::jags.model(textConnection(modelstring), data =
                                     data_priorpred)
  update(model_pat, n.iter = MCMCburnin)
  #nconverge_pat=coda.samples(model=model_pat,variable.names=c("m","c","tau_par","class","probdef"),n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
  output_post = rjags::coda.samples(
    model = model_pat,
    n.iter = MCMCout * MCMCthin,
    thin = MCMCthin,
    variable.names = c(
      "m",
      "c",
      "m_pred",
      "c_pred",
      "mu_m0",
      "tau_m0",
      "mu_c0",
      "tau_c0",
      "tau_norm",
      "tau_def",
      "Ysyn_norm",
      "Ysyn_def",
      "class",
      "probdef"
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
      'mu_m0',
      "tau_m0",
      "mu_c0",
      "tau_c0",
      "tau_norm",
      "tau_def",
      "Ysyn_norm",
      "Ysyn_def",
      "probdef"
    )
  )

  post_out = as.data.frame(output_post[[1]])
  prior_out = as.data.frame(output_prior[[1]])

  summ_pat = summary(output_post)
  # save every output from the Bernoulli dist for TGo
  classifs = post_out[, grepl("class", colnames(post_out))]
  # classifs_pat = summ_pat$statistics[grepl("class",rownames(summ_pat$statistics)),"Mean"]

  post_names = colnames(post_out)

  post = post_out[, c(
    paste0("m[", 1:(nCrl + nPt), "]"),
    paste0("c[", 1:(nCrl + nPt), "]"),
    "m_pred",
    "c_pred",
    "mu_m0",
    "tau_m0",
    "mu_c0",
    "tau_c0",
    "tau_norm",
    "tau_def",
    "probdef"
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
    paste0("m[", 1:(nCrl + nPt), "]"),
    paste0("c[", 1:(nCrl + nPt), "]"),
    "m_pred",
    "c_pred",
    "mu_m0",
    "tau_m0",
    "mu_c0",
    "tau_c0",
    "tau_norm",
    "tau_def",
    "probdef"
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
    post = post,
    postpred = postpred,
    prior = prior,
    priorpred = priorpred,
    classif = classifs
  )

  # output_saver(outroot, out_list, folder, pat_only=TRUE)
  return(out_list)
}
