#' @title Transparent Colour Creator
#'
#' @description Defining a transparent version of a given colour.
#'
#' @details
#' Either a valid colour name or RGB values must be passed to the function. If
#' both `name` and `rgbVals` are passed to the function, `name` is used to create
#' the colour.
#'
#' @param col A valid colour name, as listed by [colors] function. OR a positive
#' integer indexing the colour listed in the [palette] function. OR a string of
#' the form "#rrggbb", see [rgb].
#' @param alpha The level of transparency, a `numeric` in the range [0,1].
#' @param maxValue The maximum colour value, default to 255.
#'
#' @return A character vector of the hexidecimal value defining the colour.
#'
#' @export
#'
#' @importFrom grDevices rgb
#' @importFrom grDevices col2rgb
#'
#' @examples
#' colOne = alphaCol(col="blue", alpha=0.5)
#' colTwo = alphaCol(col=6, alpha=0.1)
#' colThree = alphaCol(rgb(0.6, 0.1, 0.6), alpha=0.8)
#'
#'
alphaCol = function(col, alpha=0.5) {
  if (alpha<0.0 | alpha>1.0 ) {
    stop("`alpha` must be in the range [0,1].")
  }

  if (is.character(col) && (col%in%colors() | str_sub(col,1,1)=='#')) {
    if (str_sub(col,1,1)=='#' & nchar(col)==9) {
      col = str_sub(col, 1, 7)
    }
    rgbVals = t( col2rgb(col, alpha=FALSE) )
    return (rgb(red=rgbVals[1], green=rgbVals[2], blue=rgbVals[3], alpha=alpha*255, maxColorValue=255))
  } else if (is.numeric(col) && (ceiling(col)==floor(col) & col>0)) {
    rgbVals = t( col2rgb(col, alpha=FALSE) )
    return (rgb(red=rgbVals[1], green=rgbVals[2], blue=rgbVals[3], alpha=alpha*255, maxColorValue=255))
  } else {
    stop ("`col` was not appropriately defined, see documentation.")
  }
}
