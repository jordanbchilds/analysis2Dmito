#' @title 2Dmito Plot with Classifications
#'
#' @description
#' Plots control and patient data with patient data classifications on a scale
#' from blue to red (healthy to deficient).
#'
#' @param dataMats A list of matrices containing the control data, named `ctrl`,
#' and the patient data, `pts`. Where the first column represents the data
#' for protein along the x-axis of the 2Dmito plot and the second column is the
#' y-axis data.
#' @param classifs A numeric vector where the i-th element is the probability
#' that the fibre in the i-th observation of the patient data matrix is
#' deficient The default is `NULL`, if this is the case the fibres are plotted as
#' green.
#' @param postpred A [data.frame] whose columns contain the 95-th percentile
#' posterior predictive interval and the corresponding x-axis values for a
#' linear regression model for this dataset. The columns should be labelled;
#' `mtiochan`, `lwrNorm`, `medNorm` and `uprNorm`. This is the form of the
#' output given by [analysis2Dmito::inference], in the `POSTPRED`.
#' @param colAlpha Transparency value of the plotted points. A numeric between
#' 0.0 and 1.0.
#'
#' @return NULL.
#'
#' @examples
#' exampleData = get_exampleData()
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
#' pat = "P01"
#'
#' data_mat = getData_mats(exampleData, channels=c(mitochan, chan), ctrlID=ctrlid, pts=pat, getIndex=TRUE)
#' data_mat$ctrl = data_mat$ctrl
#' data_mat$pts = data_mat$pts
#'
#' infOut = inference( data_mat )
#' class = apply(infOut$classif, 2, mean)
#'
#' classif_plot(dataMats=data_mat, classifs=class, postpred=infOut$postpred, xlab=mitochan, ylab=chan)
#'
#' @export
classif_plot = function(dataMats,
                        classifs = NULL,
                        postpred = NULL,
                        colAlpha = 1.0,
                        ...) {
  xlims = range(c(dataMats$ctrl[, 1], dataMats$pts[, 1]))
  ylims = range(c(dataMats$ctrl[, 2], dataMats$pts[, 2]))

  plot(NULL, xlim = xlims, ylim = ylims, ...)
  points(dataMats$ctrl,
         pch = 20,
         col = alphaBlack(colAlpha*0.5),
         ...)
  if (is.null(classifs)) {
    cols = alphaGreen(colAlpha)
  } else {
    cols = classcols(classifs, alphaLevel=colAlpha)
  }
  points(dataMats$pts,
         pch = 20,
         col = cols,
         ...)
  lines(
    postpred[, "mitochan"],
    postpred[, "lwrNorm"],
    col = alphaGreen(1.0),
    lty = 2,
    ...
  )
  lines(
    postpred[, "mitochan"],
    postpred[, "medNorm"],
    col = alphaGreen(1.0),
    lty = 1,
    ...
  )
  lines(
    postpred[, "mitochan"],
    postpred[, "uprNorm"],
    col = alphaGreen(1.0),
    lty = 2,
    ...
  )
}
