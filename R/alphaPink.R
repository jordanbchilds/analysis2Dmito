#' @title Transparent Pink
#'
#' @description Create a pink with a specified level of opacity.
#'
#' @details The RGB values for alphaPink are 255,62 and 150, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaPink = function(alpha) rgb(255/255,62/255,150/255, alpha)
