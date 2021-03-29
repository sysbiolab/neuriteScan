#' Merge image files into RGB
#'
#' Function that performs files merging to RGB format.
#'
#' @param red A path to the file containing an image for the red channel.
#' @param green A path to the file containing an image for the green channel.
#' @param blue A path to the file containing an image for the blue channel.
#'
#' @return The image in RGB format.
#'
#' @seealso \code{\link{channelSegmentation}}
#' @seealso \code{\link[EBImage{channel}]}
#' @seealso \code{\link[EBImage{readImage}]}
#' @seealso \code{\link[EBImage{rgbImage}]}
#'
#' @examples
#' blue <- system.file("images", "nucleus.jpg", package = "neuriteScan")
#' green <- system.file("images", "cytoplasm.jgp", package = "neuriteScan")
#' img <- readRgbChannels(green = green, blue = blue)
#' EBImage::display(img)
#'
#' @importFrom EBImage channel readImage rgbImage
#' @export
#'


readRgbChannels<-function(red=NULL,green=NULL,blue=NULL){
  img<-NA
  tryReadRgb<-function(red,green,blue){
    if(!is.null(red)){
      red<-EBImage::readImage(red)
      red<-EBImage::channel(red,"red")
      red<-EBImage::channel(red,"grey")
    }
    if(!is.null(green)){
      green<-EBImage::readImage(green)
      green<-EBImage::channel(green,"green")
      green<-EBImage::channel(green,"grey")
    }
    if(!is.null(blue)){
      blue<-EBImage::readImage(blue)
      blue<-EBImage::channel(blue,"blue")
      blue<-EBImage::channel(blue,"grey")
    }
    img<-EBImage::rgbImage(red,green,blue)
  }
  try(img<-tryReadRgb(red,green,blue),TRUE)
  img
}
