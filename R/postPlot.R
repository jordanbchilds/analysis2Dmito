#' @title Plot Posteriors, Predictions and Classifications for [analysis2Dmito::inference()]
#'
#' @description
#' For plotting the model output of the [analysis2Dmito::inference()]
#' function. Plotted are the priors/parent distributions and posteriors of
#' variables in the model as well as the 95\% predictive posterior for the
#' like-control linear model and the classifications.
#'
#' @param post Either a [data.frame] containing the posterior draws of the model or the a filepath to a saved file, whose columns are the posterior draws from the output. Both should be output of the [analysis2Dmito::inference()].
#' @param prior Either a [data.frame] containing draws from the prior beliefs about model parameters or the a filepath to a saved file, whose columns are the prior draws. Both should be output from the output of the [analysis2Dmito::inference()]. The default is NULL. If this is true, the prior beliefs are not added to the plots.
#' @param postpred A [data.frame] containing posterior predictions from the model, as outputed from [analysis2Dmito::inference()], or a filepath to a saved file of this output.
#' @param priorpred A [data.frame] containing prior predictions from the model, as outputed from [analysis2Dmito::inference()], or a filepath to a saved file of this output. The default is NULL. If this is true, no prior predictions will be plotted.
#' @param dataMats A character vector describing the names of the variables to be plotted as passed to [rjags] during inference. The default list is the full list of parameters in the model.
#' @param classifs A numeric vector of the probabilies that each fibre in the patient datasets are not like control.
#' @param var.names A named character vector the name of the variables to be plotted on the density plots. The default is NULL, if this is TRUE then names as given in 'var.names' will be plotted.
#' @param xlabs A character vector of names of variables as they should be added to the plot.
#'
#' @return NULL.
#'
#' @export
postPlot = function(post,
                    prior = NULL,
                    postpred,
                    priorpred = NULL,
                    dataMats,
                    classifs,
                    var.names = c("mu_m0",
                                  "tau_m0",
                                  "tau_norm",
                                  "mu_c0",
                                  "tau_c0",
                                  "probdef",
                                  "m",
                                  "c"),
                    xlabs = NULL,
                    ...) {
  if (!is.null(xlabs)) {
    if (is.null(names(xlabs)))
      names(xlabs) = var.names
  } else {
    xlabs = var.names
    names(xlabs) = var.names
  }

  if (is.character(post)) {
    if (file.exists(post)) {
      post = data.table::fread(post, header = TRUE)
      post = as.data.frame(post)
    } else {
      stop(paste0("The file, ", post, ", does not exist."))
    }
  }

  if (is.character(postpred)) {
    if (file.exists(postpred)) {
      postpred = data.table::fread(postpred, header = TRUE)
      postpred = as.data.frame(postpred)
    } else {
      stop(paste0("The file, ", postpred, ", does not exist."))
    }
  }

  if (!is.null(prior)) {
    if (is.character(prior)) {
      if (file.exists(prior)) {
        prior = data.table::fread(prior, header = TRUE)
        prior = as.data.frame(prior)
      } else {
        message(paste0("The file, ", prior, ", does not exist."))
      }
    }
  }

  if (!is.null(priorpred)) {
    if (is.character(priorpred)) {
      if (file.exists(priorpred)) {
        priorpred = data.table::fread(priorpred, header = TRUE)
        priorpred = as.data.frame(priorpred)
      } else {
        message(paste0("The file, ", priorpred, ", does not exist."))
      }
    }
  }

  colnames = colnames(post)
  for (var in var.names) {
    if (sum(grepl(paste0(var, "\\["), colnames)) > 0) {
      post_var = post[, grepl(paste0(var, "\\["), colnames)]
      post_dens = list()

      priorParent = prior[, min(grep(paste0(var, "\\["), colnames))]
      priorParent_dens = density(priorParent)

      postParent = rnorm(nrow(post_var), post[, paste0("mu_", var, "0")], 1 /
                           sqrt(post[, paste0("tau_", var, "0")]))
      postParent_dens = density(postParent)

      xlim_varMin = min(postParent_dens$x)
      xlim_varMax = max(postParent_dens$x)
      ylim_var = max(postParent_dens$y)

      for (param in colnames(post_var)) {
        post_dens[[paste(param)]] = density(post_var[[param]])

        ylim_var = max(c(ylim_var, post_dens[[param]]$y))
        xlim_varMax = max(c(xlim_varMax, post_dens[[param]]$x))
        xlim_varMin = min(c(xlim_varMin, post_dens[[param]]$x))
      }
      if (is.null(xlabs)) {
        plot(
          NA,
          main = "",
          xlab = paste(var),
          ylab = "",
          cex.lab = 3,
          xlim = c(xlim_varMin, xlim_varMax),
          ylim = c(0, ylim_var)
        )
      } else {
        plot(
          NA,
          main = "",
          xlab = xlabs[var],
          ylab = "",
          cex.lab = 3,
          xlim = c(xlim_varMin, xlim_varMax),
          ylim = c(0, ylim_var)
        )
      }
      lines(
        priorParent_dens,
        lwd = 5,
        col = alphaPink(1.0),
        lty = 4
      )
      lines(
        postParent_dens,
        lwd = 5,
        col = alphaGreen(1.0),
        lty = 4
      )

      ind = grep(paste0(var, "\\["), colnames, value = TRUE)
      ind = gsub("\\]", "", gsub(paste0(var, "\\["), "", ind))
      max_ind = max(as.numeric(ind))

      for (param in colnames(post_var)[paste0(var, "[", max_ind, "]") != colnames(post_var)]) {
        lines(post_dens[[param]], lwd = 4, col = alphaGreen(0.4))
      }
      lines(post_dens[[paste0(var, "[", max_ind, "]")]], lwd = 5, col =
              alphaGreen(1.0))

    } else {
      post_den = density(post[, paste(var)])
      prior_den = density(prior[, paste(var)])
      xlims = range(c(prior_den$x, post_den$x))
      yMax = max(c(prior_den$y, post_den$y))

      plot(
        NA,
        xlim = xlims,
        ylim = c(0, yMax),
        xlab = xlabs[var],
        cex.lab = 3,
        ylab = "",
        main = ""
      )
      lines(prior_den, lwd = 4, col = alphaPink(1.0))
      lines(post_den, lwd = 4, col = alphaGreen(1.0))
    }
  }

  xlims = range((c(dataMats$ctrl[, 1], dataMats$pts[, 1])))
  ylims = range((c(dataMats$ctrl[, 2], dataMats$pts[, 2])))

  plot(
    NULL,
    xlab = mitochan,
    ylab = chan,
    main = "",
    xlim = xlims,
    ylim = ylims
  )
  # rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray93")
  points(dataMats$ctrl,
         pch = 20,
         cex = 1.0,
         col = alphaBlack(0.05))
  points(dataMats$pts,
         pch = 20,
         cex = 1.0,
         col = classcols(classifs))

  lines(
    postpred[,"mitochan"],
    postpred[,"lwrNorm"],
    lty = 2,
    col = alphaGreen(1.0),
    lwd = 2
  )
  lines(
    postpred[,"mitochan"],
    postpred[,"medNorm"],
    lty = 1,
    col = alphaGreen(1.0),
    lwd = 4
  )
  lines(
    postpred[,"mitochan"],
    postpred[,"uprNorm"],
    lty = 2,
    col = alphaGreen(1.0),
    lwd = 2
  )

  return(NULL)
}
