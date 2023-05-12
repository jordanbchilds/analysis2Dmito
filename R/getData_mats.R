#' @title Data preparation for [analysis2Dmito::inference()].
#'
#' @description The function aggregates the data into the form needed for [analysis2Dmito::inference]. To be able to do this the data passed to the function must in a specific form, see details for more info.
#' The data.frame object passed to the function must be in long form and have the following columns; 'value', 'channel', 'sampleID' and 'fibreID'. Where 'value' is the protein expression value, 'channel' is the protein or channel that the value expresses, 'sampleID' is the identifying name associated with the sample the expression is from and'fibreID' is the fibre identification from that sample. Other columns can be present but are not needed.
#' To be able to transform the data into long form, we suggest using the [tidyverse] package and the [tidyr::pivot_longer()] function.
#'
#' @param data A data.frame object of the protein expression data. See details for the format the data should be in to be able to parse to the data.
#' @param channels A vector of strings or a single string containing the channels which will form the columns of the returned matrix.
#' @param ctrlID A character vector denoting the ID for the control samples in the 'sampleID' column of data.
#' @param pts A vector of subject names or a single subject name, whose expression levels are wished to be returned. The defualt value is NULL, if this is the case the protein expression levels for all subjects within the dataset are returned.
#' @param ctrl_only A boolean variable indicating whether to return only patient data, the default is FALSE. If this is TRUE and getIndex=FALSE the function will return a single matrix otherwise it will return a list containing all outputs.
#' @param getINDEX A boolean parameter, if TRUE two vectors of indexes, one for the control data and one for the patient data, are returned. Each indicating which observations came from the same sample by grouping them from 1 upto the number of samples in the data requested. The default is FALSE.
#'
#' @returns If ctrl_only=TRUE and getIndex=FALSE then a matrix is returned of just the control sample fibre expressiosn. The columns of the matrix correspond to the channels requested. Otherwise a list is returned containing the requested elements.
#' - ctrl : The control sample data matrix.
#' - pts : The patient sample data matrix.
#' - indexCtrl : A vector of control sample indexes.
#' - indexPts : A vector of patient sample indexes.
#'
#' @export
#'
#' @examples
#' data(exampleData)
#' # the measure of mitochondrial mass - the x-axis of the 2D mito plot
#' mitochan = "raw_porin"
#'
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#' sbj = unique(exampleData$sampleID)
#' ctrlid = c("C01", "C02", "C03", "C04", "C05")
#' pts = sbj[ !(sbj %in% ctrlid) ]
#' chan = channels[1]
#' pat = pts[1]
#' data_mat = getData_mats(exampleData, cord=c(mitochan, chan), ctrlID=ctrlid, pts=pat, getIndex=TRUE)
#'


getData_mats = function(data,
                        channels,
                        ctrlID,
                        pts = NULL,
                        ctrl_only = FALSE,
                        getIndex = FALSE) {
  sbj = sort(unique(data$sampleID))
  if (is.null(pts)) {
    pts = sbj[!(sbj %in% ctrlID)]
  }
  ctrl_mat = NULL
  indexCtrl = NULL
  indexPat = NULL
  ind = 1

  for (cr in ctrlID) {
    ctrl_data = data[data$sampleID == cr,]
    mat = NULL
    for (chan in channels) {
      mat = cbind(mat, ctrl_data[ctrl_data$channel == chan, "value"])
    }
    ctrl_mat = rbind(ctrl_mat, mat)
    indexCtrl = c(indexCtrl, rep(ind, nrow(mat)))
    ind = ind + 1
  }

  if (ctrl_only) {
    if (getIndex) {
      return(list(ctrl = ctrl_mat, index = index))
    } else {
      return(ctrl_mat)
    }
  } else {
    pts_mat = NULL
    for (pat in pts) {
      pat_data = data[data$sampleID == pat,]
      mat = NULL
      for (chan in channels) {
        mat = cbind(mat, pat_data$value[pat_data$channel == chan])
      }
      pts_mat = rbind(pts_mat, mat)
      indexPat = c(indexPat, rep(ind, nrow(mat)))
      ind = ind + 1
    }
  }

  if (getIndex) {
    return(list(
      ctrl = ctrl_mat,
      pts = pts_mat,
      indexCtrl = indexCtrl,
      indexPat = indexPat
    ))
  } else {
    return(list(ctrl = ctrl_mat, pts = pts_mat))
  }
}
