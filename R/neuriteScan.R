#' Scan objects in a segmented image
#'
#' Function that performs neurite scanning in a segmented image.
#'
#' @param segImage A segmented image obtained using the channelSegmentation
#' function. That image will be used as input for neurite scanning.
#' @param gridscan An integer value to define the size of the grid that will be
#' applicated to the segmented image.
#' @param lengthLimit An integer value to define the maximum size to be consider
#' as a neurite length. The value is in number of pixels.
#'
#' @return A list containing the final results. This list contains
#' the results for the total count of neurites, the total count of cells, the
#' number of neurites per cytoplasm area and the number of neurites per cell.
#'
#' @seealso \code{\link{channelSegmentation}}
#'
#' @examples
#' results <- neuriteScan(segImage,gridscan=10, lengthLimit = 30)
#'
#' @export
#'


neuriteScan <- function(segImage, gridscan=10, lengthLimit = 30){
  cytomask <- segImage$cytomask
  cytomask [cytomask>0] <- 1
  nucmask <- segImage$nucmask
  if(length(gridscan)==1)gridscan<-c(gridscan,gridscan)
  nObj <- length(unique(as.numeric(nucmask)))-1
  nc <- ncol(cytomask)
  nr <- nrow(cytomask)
  ic <- as.integer(c(seq.int(1,nc,by=nc/gridscan[1]),nc))
  ir <- as.integer(c(seq(1,nr,by=nr/gridscan[2]),nr))
  sumObjCol <- sapply(ic,function(i){
    tp <- cytomask[1:nr,i]
    sum(.bitcounts(tp)<lengthLimit)
  })
  sumObjRow <- sapply(ir,function(i){
    tp <- cytomask[i,1:nc]
    sum(.bitcounts(tp)<lengthLimit)
  })

  neuriteTotal <- mean(c(sum(sumObjCol),sum(sumObjRow)))
  kpix <- 1000
  areacyto <- (sum(cytomask==1))/kpix
  resnuc <- neuriteTotal/nObj
  rescyto <- neuriteTotal/areacyto

  res_neurite <- list(neuriteCount = neuriteTotal,celCount = nObj,
                      neuriteCyto = rescyto,neuriteCel = resnuc)
  neuriteResults <- data.frame(t(sapply(res_neurite,c)))
  return(neuriteResults)
}


.bitcounts <- function(x){
  g <- rep(0, length(x))
  lab <- 1
  for(i in 2:length(x)){
    if(x[i]==1){
      if(x[i-1]==0){
        g[i] <- lab; lab <- lab + 1
      } else {
        g[i]=g[i-1]
      }
    }
  }
  gg <- split(x, g)[-1]
  nobj <- as.numeric(lapply(gg, length))
  length(nobj)
  return(nobj)
}

