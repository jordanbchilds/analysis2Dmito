#' @title Save List of Dataframes
#'
#' @description
#' Saves all data frames stored within a list by suffixing the name of the
#' dataframe, as written in the list,  onto the end of a root file path. The
#' files are saved as .txt files defaulted to the tab separator.
#'
#' @details
#' The function uses the [write.table()] function to save the files.
#'
#' @param df_list A list of dataframes to be saved.
#' @param root A character vector of the root paths to which the names of the list will be added. If this is a single character then the same root will be used for files.
#' @param replace A boolean variable. If TRUE the function will overwrite any file that exists under the same name.
#' @param sep A character of an appropriate separator value.
#' @param col.names,row.names A boolean, passed directly to [write.table()].
#' @param ... Any additional parameters to be passed to write.table.
#'
#' @examples
#' data_list = list(fileOne = data.frame(a=1:20, b=21:40), fileTwo=data.frame(x=1:1000, y=rnorm(1000)))
#' dir.create("dir")
#' rootFile = "data/sourceOne"
#' output_saver(data_list, rootFile, replace=FALSE)
#'
#' @importFrom utils write.table
#'
#' @export

list_saver = function(df_list,
                        root,
                        replace = TRUE,
                        sep = "\t",
                        col.names = TRUE,
                        row.names = FALSE,
                        ...) {
  if (length(root) == 1) {
    root = rep(root, length(df_list))
  }
  names(root) = names(df_list)
  for (df in names(df_list)) {
    filePath = paste0(paste(root[df], df, sep = "_"), ".txt")
    if (file.exist(filePath) && !replace) {
      message("File exists and has NOT been replaced.")
    } else {
      utils::write.table(
        df_list[[df]],
        filePath,
        sep = sep,
        row.names = row.names,
        col.names = col.names,
        ...
      )
    }
  }
}
