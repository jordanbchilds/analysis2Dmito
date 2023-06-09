#' @title
#' Retrive IMC or QIF Data
#'
#' @description
#' The function aggregates the data into the form needed for [analysis2Dmito::inference()]. To be able to do this the data passed to the function must in a specific form, see details for more info.
#' The data.frame object passed to the function must be in long form and have the following columns; 'value', 'channel', 'sampleID' and 'fibreID'. Where 'value' is the protein expression value, 'channel' is the protein or channel that the value expresses, 'sampleID' is the identifying name associated with the sample the expression is from and'fibreID' is the fibre identification from that sample. Other columns can be present but are not needed.
#' To be able to transform the data into long form, we suggest using the [tidyr] package and the [tidyr::pivot_longer()] function.
#'
#' @details
#' The file is read using the [data.table::fread] function. If the file does not
#' exist here then an error is returned stating as such.
#'
#' @details
#' For NPC correction to take place the `filenames` column in the dataset must
#' be apprropriately labelled as "<sampleID> NPC" or "<sampleID> OXPHOS"
#'
#'
#' @param fname A character describing the filepath to the dataset from the current wokring directory.
#' @param cord A character vector of the channels of interest in the dataset.
#' @param correctnpc A boolean. If true the expression values are corrected NPC.
#' @param noNPCchannels If `correctnpc=TRUE` then the
#'
#' @importFrom data.table fread
#' @importFrom stats aggregate
#'
#' @return
#'
#' @export

getData = function(fname,
                   cord,
                   correctnpc = TRUE,
                   noNPCchannels = NULL) {
  if (file.exists(fname)) {
    dat = data.table::fread(
      fname,
      sep = "\t",
      stringsAsFactors = FALSE,
      header = TRUE
    )
  } else {
    stop(paste(fname, "does not exist in current or parent directory."))
  }

  dat = dat[dat$Channel %in% cord, ]
  dat$fn = gsub(" NPC", "", dat$Filename)
  dat$fn = gsub(" OXPHOS", "", dat$fn)
  dat$pch = paste(dat$fn, dat$ch)

  if (correctnpc) {
    npc = dat[grepl("NPC", dat$Filename),]
    oxphos = dat[grepl("OXPHOS", dat$Filename),]
    agg = stats::aggregate(npc$Value, by = list(npc$pch), mean)
    lu = agg$x
    names(lu) = gsub("_L._C._S._R.", "", agg$Group.1)

    oxphos$Filename = gsub(" OXPHOS", "", oxphos$Filename)

    newvals = pmax(1.0, oxphos$Value - lu[gsub("_L._C._S._R.", "", oxphos$pch)])
    if (!is.null(noNPCchannels))
      newvals[oxphos$chan %in% noNPCchannels] = oxphos$Value[oxphos$chan %in% noNPCchannels]
    oxphos$Value = newvals

    oxphos$Filename = oxphos$fn
    oxphos$fn = NULL
    oxphos$pch = NULL
    dat = oxphos
  } else {
    dat = dat[!grepl("NPC", dat$Filename), ]
    dat$Filename = gsub(" OXPHOS", "", dat$Filename)
  }
  dat$type = "Mean intensity"
  return(data.frame(dat, stringsAsFactors = FALSE))
}
