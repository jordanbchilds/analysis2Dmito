#' @title Transparent Green
#'
#' @description Create a green with a specified level of opacity.
#'
#' @details The RGB values used for alphaGreen are 0,150 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @importFrom grDevices rgb
#'
#' @export
alphaGreen = function(alpha) rgb(0,150/255,0, alpha)
