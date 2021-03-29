#' Returning some grids in the resulting image
#'
#' Function that performs displaying grid lines in an image object.
#'
#' @param mask A binary image used as mask, such as mask contained in the
#' resulting list of the channelSegmentation function.
#' @param gridscan An integer value to define the size of the grid that will be
#' applicated to the segmented image.
#' @param thickness An integer value to define the thickness of the grid lines.
#'
#' @return An binary image containing grid lines.
#'
#' @seealso \code{\link{channelSegmentation}}
#' @seealso \code{\link{neuriteScan}}
#'
#' @examples
#' grid <- getGridLines(segImage$cytomask)
#' paintedImage <- segImage$paintedImage
#' paintedImage[grid==1] <- 1
#' EBImage::display(paintedImage)
#'
#' @export
#'


getGridLines<-function(mask, gridscan=10, thickness=1){
  if(length(gridscan)==1)gridscan<-c(gridscan,gridscan)
  if(thickness<0)thickness=0
  nc<-ncol(mask)
  nr<-nrow(mask)
  ic<-as.integer(c(seq.int(1,nc,by=nc/gridscan[1])))[-1]
  ir<-as.integer(c(seq(1,nr,by=nr/gridscan[2])))[-1]
  gridmask<-array(0,dim=dim(mask))
  gridmask[ir,]<-1L
  gridmask[,ic]<-1L
  for(i in c(1:thickness)){
    gridmask[ir+i,]<-1L
    gridmask[,ic+i]<-1L
    gridmask[ir-i,]<-1L
    gridmask[,ic-i]<-1L
  }
  return(gridmask)
}
