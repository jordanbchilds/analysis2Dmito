#' @title Transparent Blue
#'
#' @description Create a blue with a specified level of transparency.
#'
#' @details The RGB values used for alphaBlue are 0,0 and 128, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque. Default = 0.5.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @importFrom grDevices rgb
#'
#' @export
alphaBlue = function(alpha=0.5) rgb(0,0,128/255, alpha=alpha)
