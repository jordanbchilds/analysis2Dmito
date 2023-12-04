#' @title Save List of Dataframes
#'
#' @description
#' Saves all data frames stored within a list by suffixing the name of the
#' dataframe, as written in the list,  onto the end of a root file path. The
#' files are saved as .txt files defaulted to the tab separator.
#'
#' @details
#' The function uses [write.table] to save the files. By default files are saved
#' with a `.txt` extension and overwrite any file with the same name. The name
#' of the data frame is combined with the root file path by a `_` by default
#' before adding the extension e.g. if the name of the data frame was
#' "datasetOne" and the root file path was "data/processed" then the complete
#' file path, under which the file is saved, would be
#' "data/processed_datasetOne.txt".
#'
#' @param df_list A list of data.frame objects to be saved.
#' @param root A character vector of the root path to which the names of the
#' list will be added. If this is a single character then the same root will be
#' used for files.
#' @param replace A boolean variable. If `TRUE` the function will overwrite any file that exists under the same name, default = `TRUE`.
#' @param nameSep A string which is used to join the data name and the file path, its default is "_". A vector of different joiners is allowed.
#' @param sep A character of an appropriate separator value, default=`\t`.
#' @param col.names,row.names A boolean, passed directly to [write.table], default `TRUE` for column names and `FALSE` for row names.
#' @param ... Any additional parameters to be passed to [write.table].
#'
#' @examples
#' data_list = list(fileOne = data.frame(a=1:20, b=21:40), fileTwo=data.frame(x=1:1000, y=rnorm(1000)))
#' dir.create("dir")
#' rootFile = "data/sourceOne"
#' output_saver(data_list, rootFile, replace=FALSE, nameSep="__")
#'
#' @importFrom utils write.table
#'
#' @export

list_saver = function(df_list,
                      root = "",
                      ext = ".txt",
                      replace = TRUE,
                      rootSep = "",
                      sep = "\t",
                      col.names = TRUE,
                      row.names = FALSE,
                      ...) {
  if (root==""){ rootSep = "" }
  if (length(rootSep) == 1) {
    nameSep = rep(rootSep, length(df_list))
  }
  if (length(root) == 1) {
    root = rep(root, length(df_list))
  }

  file_names = names(df_list)
  for (ii in seq_along(file_names)) {
    filePath = paste0(paste(root[ii], file_names[ii], sep = rootSep[ii]), ext)
    if (file.exists(filePath) && !replace) {
      message("File exists and has NOT been replaced.")
    } else {
      utils::write.table(
        df_list[[ii]],
        filePath,
        sep = sep,
        row.names = row.names,
        col.names = col.names,
        ...
      )
    }
  }
}

