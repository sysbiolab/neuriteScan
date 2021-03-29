#' Image Pre-Processing
#'
#' Function that performs image pre-processing with the aim of minimize the
#' noise and contrast problems that might be associated with a poor staining or
#' bad image quality.
#'
#' This function uses skimage and numpy functions. Skimage and numpy are
#' modules imported from python via reticulate package.
#'
#' @param img An image object as input to be pre-processed.
#'
#' @return An improved image with adjusted contrast and less noise.
#'
#' @seealso \code{\link[EBImage{clahe}]}
#'
#' @examples
#' image <- system.file("images", "cells.jpg", package = "neuriteScan")
#' img <- EBImage::readImage(image)
#' img <- preProcessing(img)
#' EBImage::display(img)
#'
#' @importFrom reticulate import
#' @export
#'
#'


preProcessing <- function(img){
  img <- EBImage::clahe(img)
  # numpy <- reticulate::import("numpy")
  # skimage <- reticulate::import("skimage")
  # skirest <- reticulate::import("skimage.restoration")
  #
  # sigma_est <- numpy$mean(skirest$estimate_sigma(img, multichannel = TRUE))
  # nlm <- skirest$denoise_nl_means(img, h = as.integer(1.15 * sigma_est),
  #                                 fast_mode = TRUE, patch_size = 5,
  #                                 patch_distance = 6, multichannel = TRUE)
  # img <- Image(nlm, colormode = Color)
  return(img)
}
