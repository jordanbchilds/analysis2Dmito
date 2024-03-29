#' @title Transparent Pink
#'
#' @description Create a pink with a specified level of transparency.
#'
#' @details The RGB values for alphaPink are 255,62 and 150, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque. Default = 0.5.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @importFrom grDevices rgb
#'
#' @export
alphaPink = function(alpha=0.5) rgb(255/255,62/255,150/255, alpha=alpha)
