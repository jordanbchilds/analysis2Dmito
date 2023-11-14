#' @title Transparent Black
#'
#' @description Create a black with a specified level of transparency.
#'
#' @details The RGB values for alphaBlack are 0,0,0.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely
#' transparent and 1.0 being opaque. Default = 0.5.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @importFrom grDevices rgb
#'
#' @export
alphaBlack = function(alpha=0.5) rgb(0,0,0, alpha=alpha)
