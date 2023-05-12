
#' Transparent Black
#'
#' Create a black with a specified level of opacity.
#'
#' The RGB values for alphaBlack are 0,0,0.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely
#' transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
alphaBlack = function(alpha) rgb(0,0,0, alpha)


#' Transparent Blue
#'
#' Create a blue with a specified level of opacity.
#'
#' The RGB values used for alphaBlue are 0,0 and 128, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#' @returns A character object of the hex value colour name created.
alphaBlue = function(alpha) rgb(0,0,128/255, alpha)


#' Transparent Red
#'
#' Create a red with a specified level of opacity.
#'
#' The RGB values for alphaRed are 255,0 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#' @returns A character object of the hex value colour name created.
alphaRed = function(alpha) rgb(1,0,0, alpha)


#' Transparent Green
#'
#' Create a green with a specified level of opacity.
#'
#' The RGB values used for alphaGreen are 0,150 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#' @returns A character object of the hex value colour name created.
alphaGreen = function(alpha) rgb(0,150/255,0, alpha)


#' Transparent Pink
#'
#' Create a pink with a specified level of opacity.
#'
#' The RGB values for alphaPink are 255,62 and 150, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#' @returns A character object of the hex value colour name created.
alphaPink = function(alpha) rgb(255/255,62/255,150/255, alpha)


#' Blue-Red Colour Ramp
#'
#' For creating transitions colours between red and blue, with a given transparency.
#'
#' The colour gradient lies between [analysis2Dmito::alphaBlue()] and [analysis2Dmito::alphaRed()]
#'
#' @param alphaLevel The level of opacity of the two end-point colours. A number between 0.0 and 1.0, where 0.0 is completely transparent and 1.0 being opaque.
#' @returns A character object of the hex value colour name created.
ramp_BlueRed = function(alphaLevel){
  colorRamp(c(alphaBlue(alphaLevel), alphaRed(alphaLevel)), alpha=TRUE)
}


#' Classification Colour
#'
#' The colour given to a fibre given its probability of being deficient. The
#' colours are on a scale of red to blue where red corresponds to deficiency i.e.
#' a probability of 1.0.
#'
#' The colour gradient is created using [analysis2Dmito::ramp_BlueRed()].
#'
#' @param classif A numeric vector for the deficiency probabilities. All elements must be between 0.0 and 1.0.
#' @param alphaLevel A number between 0.0 and 1.0 to control the level of the opacity of the colours.
#' @returns A character vector hex value RGB names.
classcols = function(classif, alphaLevel=1.0){
  alphaCramp = ramp_BlueRed(alphaLevel)
  rgbvals = alphaCramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}

