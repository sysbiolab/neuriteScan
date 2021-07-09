#' Image Pre-Processing
#'
#' Function that performs image pre-processing with the aim of minimize the
#' noise and contrast problems that might be associated with a poor staining or
#' bad image quality.
#'
#' This function utilizes skimage and numpy functions. Skimage and numpy are
#' modules imported from Python via reticulate package.
#'
#' @param img An image object as input to be pre-processed.
#' @param fast_mode TRUE or FALSE. When set to FALSE a spatial Gaussian
#' weighting is applied to the patches when computing patch distances. When set
#' to TRUE a faster algorithm employing uniform spatial weighting on the
#' patches is applied (default = TRUE).
#' @param patch_size An integer value for patch size (default = 5).
#' @param patch_distance An integer value for patch distance (default = 6).
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


preProcessing <- function(img, fast_mode = TRUE,
                          patch_size = 5, patch_distance = 6){
  img <- EBImage::clahe(img)
  numpy <- reticulate::import("numpy")
  skimage <- reticulate::import("skimage")
  skirest <- reticulate::import("skimage.restoration")

  sigma_est <- numpy$mean(skirest$estimate_sigma(img, multichannel = TRUE))
  nlm <- skirest$denoise_nl_means(img, h = (1.15 * sigma_est),
                                  fast_mode = fast_mode,
                                  patch_size = as.integer(patch_size),
                                  patch_distance = as.integer(patch_distance),
                                  multichannel = TRUE)
  img <- Image(nlm, colormode = Color)
  return(img)
}
