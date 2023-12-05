#' @title Plot Posteriors, Predictions and Classifications for [analysis2Dmito::inference()]
#'
#' @description
#' For plotting the model output of the [analysis2Dmito::inference]
#' function. Plotted are the priors/parent distributions and posteriors of
#' variables in the model as well as the 95\% predictive posterior for the
#' like-control linear model and the classifications.
#'
#' @param post Either a [data.frame] containing the posterior draws of the model
#' or the a filepath to a saved file, whose columns are the posterior draws from
#' the output. Both should be output of the [analysis2Dmito::inference].
#' @param prior Either a [data.frame] containing draws from the prior beliefs
#' about model parameters or the a filepath to a saved file, whose columns are
#' the prior draws. Both should be output from the output of the
#' [analysis2Dmito::inference]. The default is NULL. If this is true, the
#' prior beliefs are not added to the plots.
#' @param postpred A [data.frame] containing posterior predictions from the
#' model, as outputed from [analysis2Dmito::inference], or a filepath to a
#' saved file of this output.
#' @param priorpred A [data.frame] containing prior predictions from the model,
#' as outputed from [analysis2Dmito::inference], or a filepath to a saved file
#' of this output. The default is NULL. If this is true, no prior predictions
#' will be plotted.
#' @param dataMats A character vector describing the names of the variables to
#' be plotted as passed to [rjags] during inference. The default list is the
#' full list of parameters in the model.
#' @param classifs A numeric vector of the probabilies that each fibre in the
#' patient datasets are not like control.
#' @param var.names A named character vector the name of the variables to be
#' plotted on the density plots. The default is NULL, if this is TRUE then names
#' as given in 'var.names' will be plotted.
#' @param xlabs A character vector of names of variables as they should be added
#' to the plot.
#' @param ... Any additional parameters to be passed to the plotting functions
#' e.g. lwd, cex, etc.
#' @param mitoPlot_xlab x axis label to be placed on the final classification
#' and posterior predictive model for the patient sample. Default is "".
#' @param mitoPlot_ylab y axis label to be placed on the final classification
#' and posterior predictive model for the patient sample. Default is "".
#' @param main_title Title to be placed over at the top of the plotting window.
#' Default is "".
#' @param chains The number of chains present in the posterior draws. The
#' default is 1.
#'
#' @details
#' If multiple chains are present in the posterior files, then the density plots
#' ignore this and plot all draws as though they were one - assuming that if the
#' the chains do not converge to the same stationary distribution this would be
#' evident and produce a multi-modal plot. The posterior predictive plot does
#' not ignore this and will plot each posterior predictive over the
#' classification plot. This will indicate if there are differences in the
#' chains as the posterior predictives will different.
#'
#' @return NULL.
#'
#' @examples
#' exampleData = get_exampleData()
#' # the measure of mitochondrial mass - the x-axis of the 2D mito plot
#' mitochan = "VDAC"
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#' sbj = unique(exampleData$sampleID)
#' ctrlid = c("C01", "C02", "C03", "C04")
#' pts = sbj[ !(sbj %in% ctrlid) ]
#'
#' chan = channels[1]
#' pat = pts[1]
#'
#' data_mat = getData_mats(exampleData,
#'                         cord=c(mitochan, chan),
#'                         ctrlID=ctrlid,
#'                         pts=pat,
#'                         getIndex=TRUE)
#'
#' infOut = stan_inference( data_mat )
#' class = apply(infOut$classif, 2, mean)
#'
#' postPlot(post=infOut$POST, prior=infOut$PRIOR, classifs=class, postpred=infOut$POSTPRED,
#'          dataMats=data_mat)
#'
#' @importFrom data.table fread
#' @importFrom stats density
#' @importFrom graphics lines
#' @importFrom graphics points
#'
#' @export
postPlot = function(post,
                    prior = NULL,
                    postpred,
                    priorpred = NULL,
                    dataMats,
                    classifs,
                    var.names = c("mu_m",
                                  "tau_m",
                                  "tau_norm",
                                  "mu_c",
                                  "tau_c",
                                  "probdiff",
                                  "m",
                                  "c"),
                    xlabs = NULL,
                    mitoPlot_xlab="",
                    mitoPlot_ylab="",
                    main_title="",
                    chains=1,
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

      priorParent = prior[, paste0(var, "_pred")]
      priorParent_dens = stats::density(priorParent)

      postParent = post[, paste0(var, "_pred")]
      postParent_dens = stats::density(postParent)

      xlim_varMin = min(postParent_dens$x)
      xlim_varMax = max(postParent_dens$x)
      ylim_var = max(postParent_dens$y)

      for (param in colnames(post_var)) {
        post_dens[[paste(param)]] = stats::density(post_var[,param])

        ylim_var = max(c(ylim_var, post_dens[[param]]$y))
        xlim_varMax = max(c(xlim_varMax, post_dens[[param]]$x))
        xlim_varMin = min(c(xlim_varMin, post_dens[[param]]$x))
      }
      plot(
        NA,
        main = "",
        xlab = xlabs[var],
        ylab = "",
        xlim = c(xlim_varMin, xlim_varMax),
        ylim = c(0, ylim_var),
        ...
      )
      graphics::lines(
        priorParent_dens,
        col = analysis2Dmito::alphaPink(1.0),
        lty = 4,
        ...
      )
      graphics::lines(
        postParent_dens,
        col = analysis2Dmito::alphaGreen(1.0),
        lty = 4,
        ...
      )

      ind = grep(paste0(var, "\\["), colnames, value = TRUE)
      ind = gsub("\\]", "", gsub(paste0(var, "\\["), "", ind))
      max_ind = max(as.numeric(ind))

      for (param in colnames(post_var)[paste0(var, "[", max_ind, "]") != colnames(post_var)]) {
        graphics::lines(
          post_dens[[param]],
          col = analysis2Dmito::alphaGreen(0.4),
          ...
        )
      }
      graphics::lines(
        post_dens[[paste0(var, "[", max_ind, "]")]],
        col = analysis2Dmito::alphaGreen(1.0),
        ...
      )
    } else {
      post_den = stats::density(post[, paste(var)])
      prior_den = stats::density(prior[, paste(var)])
      xlims = range(c(prior_den$x, post_den$x))
      yMax = max(c(prior_den$y, post_den$y))

      plot(
        NA,
        xlim = xlims,
        ylim = c(0, yMax),
        xlab = xlabs[var],
        ylab = "",
        main = "",
        ...
      )
      graphics::lines(
        prior_den,
        col = analysis2Dmito::alphaPink(1.0),
        ...
      )
      graphics::lines(
        post_den,
        col = analysis2Dmito::alphaGreen(1.0),
        ...
      )
    }
  }

  xlims = range((c(dataMats$ctrl[, 1], dataMats$pts[, 1])))
  ylims = range((c(dataMats$ctrl[, 2], dataMats$pts[, 2])))

  plot(
    NULL,
    xlab = mitoPlot_xlab,
    ylab = mitoPlot_ylab,
    main = "",
    xlim = xlims,
    ylim = ylims,
    ...
  )
  graphics::points(
    dataMats$ctrl,
    pch = 20,
    col = analysis2Dmito::alphaBlack(0.05),
    ...
  )
  graphics::points(
    dataMats$pts,
    pch = 20,
    col = analysis2Dmito::classcols(classifs),
    ...
  )

  for( i in 1:chains ){
    chain_ind = ( (i-1)*1000 + 1 ):( i*1000 )
    graphics::lines(
      postpred[chain_ind,"mitochan"],
      postpred[chain_ind,"lwrNorm"],
      lty = 2,
      col = analysis2Dmito::alphaGreen(1.0),
      ...
    )
    graphics::lines(
      postpred[chain_ind,"mitochan"],
      postpred[chain_ind,"medNorm"],
      lty = 1,
      col = analysis2Dmito::alphaGreen(1.0),
      ...
    )
    graphics::lines(
      postpred[chain_ind,"mitochan"],
      postpred[chain_ind,"uprNorm"],
      lty = 2,
      col = analysis2Dmito::alphaGreen(1.0),
      ...
    )
  }
}
