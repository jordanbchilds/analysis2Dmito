#' @title Data preparation for [analysis2Dmito::inference].
#'
#' @description The function aggregates the data into the form needed for
#' [analysis2Dmito::stan_inference]. To be able to do this the data passed to the
#' function must be in a specific form, see details for more info.
#'
#' @details
#' The data frame passed to the function must be in long form and have
#' the following columns; `Value`, `Channel`, `sampleID` and `fibreID`. Where
#' 'Value' is the protein expression level, 'Channel' is the protein or channel
#' that the value expresses, `sampleID` is the identifying name associated with
#' the tissue sample on which the measurement was made and `fibreID` is the
#' fibre identification from that sample. Other columns can be present but are
#' not needed. One column which is helpful is a column identifying which samples
#' are from control subjects and which are from patients. For this to be used
#' the column must be called `sbjType` and the each value must be labelled
#' `control` or `patient`. To be able to transform the data into long form, we
#' suggest using the [tidyverse] package and the [tidyr::pivot_longer] function.
#'
#' @param data A data.frame object of the protein expression data. See details
#' for the format the data should be in.
#' @param channels A vector of strings or a single string containing the
#' channels which will form the columns of the returned matrices.
#' @param ctrlID A character vector denoting the IDs for the control samples in
#' the `sampleID` column of data. If this is not given the data must have a
#' `sbjType` column specifying which samples are `control` or `patient`.
#' @param pts A vector of subject names or a single subject name, whose
#' expression levels are to be returned. The default value is `NULL`, if this is
#' the case the protein expression levels for all samples within the dataset are
#' returned.
#' @param ctrlOnly A boolean variable indicating whether to return only control
#' subject data, the default is `FALSE.`
#' @param getINDEX A boolean parameter. If `TRUE`, a vector is returned for both
#' the control and patient data matrices indicating which observations come the
#' same sample.
#'
#'
#' @returns If `ctrlOnly=TRUE` and `getIndex=FALSE` then a matrix is returned of
#' just the control sample fibre expression otherwise a list is returned
#' containing data matrices and index vectors. The columns of the data matrices
#' correspond to the channels requested, in the order they are given in `channels`
#' argument. The possible elements in the list are:
#' - `ctrl` : The control sample data matrix.
#' - `pts` : The patient sample data matrix.
#' - `indexCtrl` : A vector of indexes for the control data matrix, where the i-th
#' element corresponds to the i-th row of the matrix, indicating which
#' observations belong to the same sample.
#' - `indexPts` : A vector of indexes for the patient data matrix, where the i-th
#' element corresponds to the i-th row of the matrix, indicating which
#' observations belong to the same sample.
#'
#' The indexes, in the `indexPts` and `indexCtrl` vectors, range from 1 to the total
#' number of control samples in the data . For example, if four control samples are in the dataset (and their
#' sample IDs are passed to the `ctrlID` argument) then observations with an index
#' of 1 will be associated with the first element in `ctrlID`, an index of 2 will
#' be observations from the second element in `ctrlID` etc. Continuing this example, the patient indexes
#' will start at 5 (as there are four control sample) and will increase
#' similarly depending on the number of patient samples being outputted, which
#' can be controlled by the `pts` argument. If the data is then going to be used
#' in with the [analysis2Dmito::stan_inference] function, then there should only
#' be one patient sample.
#'
#' @export
#'
#' @examples
#' exampleData = get_exampleData()
#' # the measure of mitochondrial mass - the x-axis of the 2D mito plot
#' mitochan = "VDAC"
#' # all channels available in the dataset
#' channelsAll = unique(exampleData[,"Channel"])
#' # remove mitochan from the channels of interest
#' channels = channelsAll[ channelsAll!=mitochan ]
#'
#' sbj = unique(exampleData$sampleID)
#' ctrlID = grep("C", sbj, value=TRUE )
#' pts = grep("C", sbj, value=TRUE, invert=TRUE)
#'
#' chan = channels[1]
#' pat = pts[1]
#'
#' data_mat = getData_mats(exampleData, channels=c(mitochan, chan), ctrlID=ctrlID, pts=pat)
#'


getData_mats = function(data,
                        channels,
                        ctrlID=NULL,
                        pts = NULL,
                        ctrlOnly = FALSE,
                        getIndex = TRUE) {
  if( is.null(ctrlID) & !is.null(data$sbjType) ){
    ctrlID = unique( data[data$sbjType=="control", "sampleID"] )
  } else if ( is.null(ctrlID) & is.null(data$sbjType)){
    stop("Must pass `ctrlID` or have a column in the data called `sbjType` which indicates which sample IDs are `control` or `patient`.")
  }
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
      mat = cbind(mat, ctrl_data[ctrl_data$Channel == chan, "Value"])
    }
    ctrl_mat = rbind(ctrl_mat, mat)
    indexCtrl = c(indexCtrl, rep(ind, nrow(mat)))
    ind = ind + 1
  }

  if (ctrlOnly) {
    if (getIndex) {
      return(list(ctrl = ctrl_mat, index = indexCtrl))
    } else {
      return(ctrl_mat)
    }
  } else {
    pts_mat = NULL
    for (pat in pts) {
      pat_data = data[data$sampleID == pat,]
      mat = NULL
      for (chan in channels) {
        mat = cbind(mat, pat_data$Value[pat_data$Channel == chan])
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
