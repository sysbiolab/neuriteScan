#' Scanning images in batch
#'
#' Function that scans images in batch mode.
#'
#' @param path The file directory containing the images for analysis.
#' @param fflag The pattern observed in the name of the files.
#' @param resdir The name of the directory that will be created for the storage
#' of the resulting files.
#' @param bandwidth An integer value to define the width of the band to be
#' considered. As higher the bandwidth higher the band to be considered during
#' segmentation.
#' @param fpicks Defines the precision of the extensions to be considered. As
#' higher the fpicks less extensions are identified in the segmentation.
#' @param lengthLimit An integer value to define the minimum size to be consider
#' as a neurite length. The value is in number of pixels.
#' @param gridscan An integer value to define the size of the grid that will be
#' applicated to the segmented image.
#'
#' @return A file directory containing the resulting segmented images along with
#' a .csv file including the quantitative results obtained by the function
#' neuriteScan.
#'
#' @seealso \code{\link{channelSegmentation}}
#' @seealso \code{\link{neuriteScan}}
#'
#'
#' @importFrom EBImage readImage writeImage distmap
#' @export
#'

batchScanning <- function(path, fflag="controle", resdir="resdir", bandwidth,
                          fpicks, lengthLimit = 30, gridscan = 10){
  resdir<-paste(path.expand(path),resdir,"/", sep="")
  dir.create(resdir)
  fi<-list.files(path=path,include.dirs=FALSE,recursive=FALSE, pattern=fflag,
                 ignore.case=TRUE, no..=TRUE)
  fo<-paste("res_",fi,sep="")
  fo<-sub(".png",".jpg",fo)
  fo<-sub(" +$", "", fo)
  fo<-sub("\\s", "", fo)
  labs<-sub(".jpg","",fo)
  labs<-tolower(labs)
  fi<-paste(path,fi,sep="")
  fo<-paste(resdir,fo,sep="")
  #---
  #run segmentation
  res<-NULL
  for(i in 1:length(fi)){
    fls<-fi[[i]]
    fls <- EBImage::readImage(fls)
    fls <- preProcessing(fls)
    segImage<-channelSegmentation(img = fls,bandwidth=bandwidth, fpicks = fpicks)
    EBImage::writeImage(segImage$paintedImage, files=fo[i], type='jpg')
    res<-rbind(res,neuriteScan(segImage, lengthLimit = lengthLimit,
                               gridscan = gridscan))
  }
  rownames(res)<-labs
  #---
  resfile<-"neuritescan"
  tp<-data.frame(round(res,3),check.rows=FALSE, check.names=FALSE)
  write.csv(tp,file=paste(resdir,resfile,".csv",sep=""))
  #---
  invisible(res)
}
