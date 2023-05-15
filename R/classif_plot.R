#' @title 2Dmito Plot with Classifications
#'
#' @description
#' Plots control and patient data with patient data classifications on a scale from blue to red.
#'
#' @param dataMats A list of matrices containing the control data, named 'ctrl', and the patient data, named 'pts'. Where the first column represents the data for protein along the x-axis of the 2Dmito plot and the second column is the y-axis data.
#' @param classifs A numeric vector where the i-th element is the probability that the fibre in the i-th row of the patient fibre matrix is not like control. The default is NULL, if this is the case the fibres are plotted as green.
#' @param postpred A [data.frame] outputted from whose columns contain the 95\% posterior predictive interval and the corresponding x-axis values for a linear regression model for this dataset. The columns should be labelled; 'mitochan', 'lwrNorm', 'medNorm' and 'uprNorm'. This is the form of the output given by [analysis2Dmito::inference()].
#'
#' @return NULL.
#'
#' @examples
#' data(exampleData)
#' mitochan = "raw_porin"
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#' sbj = unique(exampleData$sampleID)
#' ctrlid = c("C01", "C02", "C03", "C04", "C05")
#' pts = sbj[ !(sbj %in% ctrlid) ]
#' chan = channels[1]
#' pat = "P01"
#' data_mat = getData_mats(exampleData, channels=c(mitochan, chan), ctrlID=ctrlid, pts=pat, getIndex=TRUE)
#' data_mat$ctrl = log(data_mat$ctrl)
#' data_mat$pts = log(data_mat$pts)
#' infOut = inference(data_mat)
#' class = apply(infOut$classif, 2, mean)
#' classif_plot(dataMats=data_mat, classifs=class, postpred=infOut$postpred, xlab=paste0("log(", mitochan, ")"), ylab=paste0("log(",chan,")"))
#'
#' @export
classif_plot = function(dataMats,
                        classifs = NULL,
                        postpred = NULL,
                        ...) {
  xlims = range(c(dataMats$ctrl[, 1], dataMats$pts[, 1]))
  ylims = range(c(dataMats$ctrl[, 2], dataMats$pts[, 2]))

  plot(NULL, xlim = xlims, ylim = ylims, ...)
  points(dataMats$ctrl,
         pch = 20,
         col = alphaBlack(0.1))
  if (is.null(classifs)) {
    cols = alphaGreen(0.7)
  } else {
    cols = classcols(classifs)
  }
  points(dataMats$pts,
         pch = 20,
         col = cols)
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
