#' @title Diagnostic MCMC Output Plots
#'
#' @description
#' Plot the trace, autocorrelation and kernel density estimate for MCMC output.
#' For each variable present the three are plotted in a single row. The number
#' of rows in the plotting environment can be user defined.
#'
#' @details
#' If a file path is passed to the function it is read using the
#' [data.table::fread()] function and then coverted to a [data.frame]. The plotting environment is created within
#' the function using the [par()] command, it is also ended within the function.
#' The kernel density estimates (KDEs) of the priors and posteriors are
#' calculated using the default [density()] function.
#'
#'
#' @param post A file path to saved output file or a [data.frame] containing the MCMC output. Where each column represents a different variable.
#' @param prior A file path to a saved output file or a [data.frame] constraining prior draws. The kernel density estimate of this will be added to the plot of the KDE for the posterior of the variable, if nothing is passed nothing is added to the plot.
#' @param lag The maximum lag index used in the autocorrelation plot.
#' @param nRow The number of rows in the plotting environment, the default is NULL. If this is true then the number of variables will be used.
#' @param colPost An appropriate colour definition to be used for the kernel density estimate of the posterior distribution for all variables.
#' @param colPrior An appropriate colour definition to be used for the kernel density estimate of the prior distribution.
#' @param ... Any parameters passed to the plotting environment created by the [par()].
#'
#' @return NULL.
#'
#' @export
#' @importFrom data.table fread
#'
MCMCplot = function(post,
                    prior = NULL,
                    lag = 20,
                    nRow = NULL,
                    colPost = alphaPink(1.0),
                    colPrior = alphaGreen(1.0),
                    ...) {
  if (is.character(post)) {
    if (file.exists(post)) {
      post = data.table::fread(post, header = TRUE, stringsAsFactors = FALSE)
      post = as.data.frame(post)
    } else {
      stop(
        paste(
          "The file,",
          post,
          ", does not exist, check the file path and working directory."
        )
      )
    }
  }

  if (is.null(prior)) {
    message("No prior draws given - only posteriors will be plotted")
  }

  if (is.character(prior)) {
    if (file.exists(prior)) {
      prior = data.table::fread(prior, header = TRUE, stringsAsFactors = FALSE)
      prior = as.data.frame(prior)
    } else {
      stop(
        paste(
          "The file,",
          prior,
          ", does not exist, check the file path and working directory."
        )
      )
    }
  }

  col.names = colnames(post)
  if (is.null(nRow))
    nRow = length(col.names)
  op = par(mfrow = c(nRow, 3), ...)
  for (param in col.names) {
    post_vec = post[, param]
    plot(
      ts(post_vec),
      xlab = "Iteration",
      ylab = paste(param),
      main = "",
      cex.lab = 1.2
    )
    if (sum(post_vec == post_vec[1]) != length(post_vec)) {
      acf(
        post[, param],
        xlab = "lag index",
        ylab = "ACF",
        main = "",
        cex.lab = 1.2,
        lag = lag
      )
    } else {
      plot(
        NA,
        type = 'n',
        xlim = c(0, lag),
        ylim = c(0, 1),
        xlab = "lag index",
        ylab = "ACF",
        main = ""
      )
    }
    xlims = range(post[, param])
    if (!is.null(prior)) {
      xlims = range(c(post[, param], prior[, param]))
    }
    plot(
      density(post[, param]),
      xlim = xlims,
      lwd = 2,
      col = alphaGreen(1.0),
      xlab = paste(param),
      ylab = "Density",
      main = ""
    )
    if (!is.null(prior))
      lines(density(prior[, param]), lwd = 2, col = alphaPink(1.0))
  }
  par(op)

  return(NULL)
}
