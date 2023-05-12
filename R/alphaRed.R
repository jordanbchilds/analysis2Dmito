#' @title Transparent Red
#'
#' @description Create a red with a specified level of opacity.
#'
#' @details The RGB values for alphaRed are 255,0 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaRed = function(alpha) rgb(1,0,0, alpha)
