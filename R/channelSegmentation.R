#' Two-channel segmentation of an image object
#'
#' Function that runs two-channel segmentation for one image.
#'
#' @param img An image object as input to be segmented.
#' @param nwin describe
#' @param smooth describe
#' @param minShape describe
#' @param maxShape describe
#' @param channelnuc A character variable to define the channel color for
#' the nucleus staining.
#' @param channelcyto A character variable to define the channel color for
#' the cytoplasm staining.
#' @param bandwidth describe
#' @param fpicks describe
#' @param colnuc A character variable to define the color of the nucleus in the
#' resulting painted image.
#' @param colcyto A character variable to define the color of the cytoplasm in
#' the resulting painted image.
#' @param colscheme An integer to define the scheme of colors in the resulting
#' painted image.
#'
#' @return A list containing the original image, the resulting grayscale images,
#' the binary masks and the segmented final image.
#'
#' @seealso \code{\link[EBImage{channel}]}
#' @seealso \code{\link[EBImage{paintObjects}]}
#' @seealso \code{\link{preProcessing}}
#'
#' @examples
#' image <- system.file("images", "cells.jpg", package = "neuriteScan")
#' img <- EBImage::readImage(image)
#' segImage <- channelSegmentation(img=img, bandwidth=10, fpicks = 0.01)
#' EBImage::display(segImage$paintedImage)
#'
#' @importFrom EBImage channel paintObjects distmap
#' @export
#'


channelSegmentation <- function(img, nwin=3, smooth=3, minShape=40,
                                maxShape=800,channelnuc="blue",
                                channelcyto="green", bandwidth=20,
                                fpicks=0.01, colnuc='cyan', colcyto='white',
                                colscheme=3)
{
  if(class(img)!="Image")stop("'img' should be an 'Image' class object.")
  imagecyto <- EBImage::channel(img,channelcyto)
  imagenuc <- EBImage::channel(img,channelnuc)

  #segmentation
  nucmask <- .nucSegmentation(imagenuc,nwin,smooth,minShape,maxShape)
  cytomask <- .cytoSegmentation(imagecyto, nucmask, nwin, bandwidth, fpicks)

  #relabel masks
  cellmasks <- .relabel2(nucmask,cytomask)

  #image for merging
  imgc <- imagecyto
  imgc[cellmasks$nucmask>0] <- 0
  imgn <- imagenuc
  msk <- cellmasks$cytomask>0 & cellmasks$nucmask==0
  imgn[msk] <- 0

  #paint image
  if(colscheme==1){
    paintedImage = EBImage::rgbImage(green=imgc, blue=imgn)
    paintedImage = EBImage::paintObjects(cellmasks$cytomask, paintedImage,
                                         col=colcyto,thick=FALSE)
    paintedImage = EBImage::paintObjects(cellmasks$nucmask, paintedImage,
                                         col=colnuc,thick=TRUE)
  } else if(colscheme==2){
    paintedImage = EBImage::rgbImage(green=imgc, blue=imgn)
    paintedImage = EBImage::paintObjects(cellmasks$cytomask, paintedImage,
                                         col=colcyto,thick=FALSE)
    paintedImage = EBImage::paintObjects(cellmasks$nucmask, paintedImage,
                                         col=colnuc,thick=TRUE)
    paintedImage@.Data[,,1][cellmasks$cytomask>0 & cellmasks$nucmask==0]<-0.9
  } else {
    imgc[cellmasks$cytomask==0]<-imgc[cellmasks$cytomask==0]*0.5
    imgc[msk]<-imgc[msk]+0.2*(1-imgc[msk])
    imgn[cellmasks$nucmask>0]<-imgn[cellmasks$nucmask>0]+
      0.2*(1-imgn[cellmasks$nucmask>0])
    paintedImage = EBImage::rgbImage(green=imgc, blue=imgn)
    paintedImage = EBImage::paintObjects(cellmasks$cytomask, paintedImage,
                                         col=colcyto,thick=FALSE)
    paintedImage = EBImage::paintObjects(cellmasks$nucmask, paintedImage,
                                         col=colnuc,thick=TRUE)
  }

  #get res and return
  res <- list(image=img,imagenuc=imagenuc, imagecyto=imagecyto,
              nucmask=cellmasks$nucmask,cytomask=cellmasks$cytomask,
              paintedImage=paintedImage)
  message("Segmentation done!")
  return(res)

}

#re-label cells using top-down ordering
.relabel1<-function(imgmask){
  regions<-unique(as.vector(imgmask))
  regions<-regions[regions!=0]
  lb<-cbind(regions,1:length(regions))
  labs<-rep(0,max(regions))
  labs[lb[,1]]<-lb[,2]
  v1<-as.vector(imgmask)
  v2<-v1[v1>0]
  idx<-which(v1>0)
  v2<-labs[v2]
  v1[idx]<-v2
  v1<-array(v1,dim(imgmask))
  imgmask<-v1
  imgmask
}
#re-label cells using top-down ordering (2 mask)
.relabel2<-function(nucmask,cytomask){
  nucmask[cytomask==0] <- 0
  regions<-unique(as.vector(nucmask))
  regions<-regions[regions!=0]
  lb<-cbind(regions,1:length(regions))
  labs<-rep(0,max(regions))
  labs[lb[,1]]<-lb[,2]
  #re-label nucmask
  v1<-as.vector(nucmask)
  v2<-v1[v1>0]
  idx<-which(v1>0)
  v2<-labs[v2]
  v1[idx]<-v2
  v1<-array(v1,dim(nucmask))
  nucmask<-v1
  #re-label cytomask
  v1<-as.vector(cytomask)
  v2<-v1[v1>0]
  idx<-which(v1>0)
  v2<-labs[v2]
  v1[idx]<-v2
  v1<-array(v1,dim(cytomask))
  cytomask<-v1
  return(list(nucmask=nucmask,cytomask=cytomask))
}

###################################################
### Nucleus segmentation
###################################################
.nucSegmentation <- function(imgN,nwin,smooth,minShape,maxShape,threshold){
  message("Start segmentation...")
  #reverse
  imgN[,] <- 1-imgN[,]
  #--------------------------
  #local thresholding
  if(missing(threshold)){
    whitepx<-imgN==1
    imgB<-.computeBinaryImage1(imgN,nwin=nwin,whitepx=whitepx)
  } else {
    imgB<-.computeBinaryImage1(imgN,nwin=nwin,threshold=threshold)
  }
  #--------------------------
  #smooth shape and set labels
  imgB <- .removetips1(imgB,th=smooth)
  imgLB <- EBImage::bwlabel(imgB)
  #set min shape
  imgD <- EBImage::imageData(imgLB)+1
  tb <- tabulate(as.vector(imgD))
  idx <- array(tb[imgD]<minShape,dim(imgD))
  imgD[idx] <- 1;imgD=imgD-1
  #watershed segmentation
  imgWS <- EBImage::watershed(distmap(imgD))
  #set max shape
  imgD <- EBImage::imageData(imgWS)+1
  hc <- tabulate(as.vector(imgD))
  idx <- array(hc[imgD]<maxShape,dim(imgD))
  imgD[idx] <- 1;imgD<-imgD-1
  imgWS[imgD>0] <- 0
  #resegment large regions
  imgL <- imgN
  imgL[imgD==0] <- 0
  imgd <- as.vector(imgD)
  imgl <- as.vector(imgL)
  imgtp <- imgd
  segs <- unique(imgd)
  segs <- segs[segs!=0]
  if(length(segs)>0){
    sapply(segs,function(i){
      indexgpx <- which(imgd==i)
      gpx <- imgl[indexgpx]
      th <- .localth1(gv=as.vector(gpx[gpx!=1]))
      thPixel <- gpx
      thPixel[thPixel<th]<--1
      thPixel[thPixel!=-1]<-0
      thPixel[thPixel==-1]<-1
      imgtp[indexgpx]<<-thPixel
      NULL
    })
  }
  imgtp <- array(imgtp,dim(imgN))
  #watershed segmentation
  imgtp <- EBImage::watershed(distmap(imgtp))
  imgtp <- imgtp+max(imgWS)
  imgtp[imgtp==max(imgWS)]<-0
  #large segments resegmented
  imgWS <- imgWS+imgtp
  #re-set min shape
  imgD <- EBImage::imageData(imgWS)+1
  tb <- tabulate(as.vector(imgD))
  idx <- array(tb[imgD]<minShape,dim(imgD))
  imgD[idx] <- 1;imgD <- imgD-1
  imgWS <- imgD
  #fill hull
  imgWS <- EBImage::fillHull(imgWS)
  #re-label using top-down ordering
  regions <- unique(as.vector(imgWS))
  regions <- regions[regions!=0]
  if(length(regions)>0){
    imgWS <- .relabel1(imgWS)
  } else {
    message("Segmentation fault! no object has been identified.")
  }
  imgWS
}

###################################################
### Compute binary images
###################################################
.computeBinaryImage1 <- function (imgG, threshold=NULL, nwin=1, whitepx){
  message("Thresholding...")
  if (is.null(threshold)){
    imgB<- .computeThreshold1(imgG, nwin, whitepx)
  } else {
    imgB<-array(0,dim(imgG))
    imgB[which(imgG>threshold)]<-1
  }
  imgB
}

###################################################
### Compute image threshold
###################################################
.computeThreshold1 <- function(imgG, nwin, whitepx){
  if(missing(whitepx)){
    whitepx <- array(FALSE,dim(imgG))
  }
  nr <- dim(imgG)[1]
  nc <- dim(imgG)[2]
  ix <- as.integer(seq(from=1, to=nr, by=(nr-1)/(nwin)))
  jx <- as.integer(seq(from=1, to=nc, by=(nc-1)/(nwin)))
  imgB <- array(0,dim(imgG))
  sapply(1:nwin,function(i){
    sapply(1:nwin,function(j){
      vi<-ix[i]:(ix[i+1])
      vj<-jx[j]:(jx[j+1])
      gw<-imgG[vi,vj]
      ww<-whitepx[vi,vj]
      if(length(as.vector(gw[ww==FALSE]))>0){
        bthresh<-.localth1(gv=as.vector(gw[ww==FALSE]))
        bw<-array(0,dim(gw))
        bw[which(gw<bthresh)]<-1
        imgB[vi,vj]<<-bw
      }
    })
  })
  imgB
}

###################################################
### Compute local threshold
###################################################
.localth1 <- function(gv){
  gv<-round(gv*255)
  pi<-(tabulate(gv)+0.001)/length(gv)
  len<-length(pi)
  if(len>1){
    ui<-sum((1:len)*(pi[1:len]))
    pk<-function(p,k)sum(p[1:k])
    uk<-function(p,k)sum((1:k)*(p[1:k]))
    maxth<-function(p,k,u){(u*pk(p,k)-uk(p,k))^2/(pk(p,k)*(1-pk(p,k)))}
    th<-optimize(maxth,c(1,len),maximum=TRUE,p=pi,u=ui)
    th<-th$maximum/255
  } else {
    th<-0
  }
  th
}

###################################################
### Filter background on edges
###################################################
.removetips1 <- function(imgB, th=2){
  x<-EBImage::distmap(imgB,"manhattan")
  x[x>1]=2
  nr<-nrow(x);nc<-ncol(x)
  x[1,]=0;x[,1]=0;x[nr,]=0;x[,nc]=0
  nr=nr-1;nc=nc-1
  ii<-jj<-c(-1,0,1)
  sapply(2:nr,function(i){
    sapply(which(x[i,]==1),function(j){
      nn<-sum(x[ii+i,jj+j]>1)
      imgB[i,j]<<-ifelse(nn>=th,1,0)
    })
  })
  imgB
}

###################################################
### Cytoplasm segmentation
###################################################
.cytoSegmentation <- function(imagecyto,nucmask,nwin,bandwidth,
                              fpicks,fillhull=TRUE){
  cytomask <- .computeBinaryImage2(imagecyto,nucmask=nucmask,nwin=nwin,
                                   bandwidth=bandwidth, fpicks=fpicks)
  cytomask <- EBImage::propagate(cytomask, seeds=nucmask, mask=cytomask)
  if(fillhull){
    cytomask <- EBImage::fillHull(cytomask)>0
    cytomask <- EBImage::propagate(cytomask, seeds=nucmask, mask=cytomask)
  }
  cytomask
}

###################################################
### Compute binary images
###################################################
.computeBinaryImage2 <- function(imgC, nucmask, nwin, bandwidth, fpicks){

  #---add cell density information to bandwidth
  bw <- sum(nucmask!=0)/prod(dim(nucmask))+bandwidth

  #---set a ceiling for bandwidth
  bw <- min(c(bw,100))

  #---set window index
  nr <- dim(imgC)[1];nc<-dim(imgC)[2]
  ix <- as.integer(seq(from=1, to=nr, by=(nr-1)/(nwin)))
  jx <- as.integer(seq(from=1, to=nc, by=(nc-1)/(nwin)))

  #---set cell body
  imgB1 <- array(0,dim(imgC))
  for(i in 1:nwin){
    for(j in 1:nwin){
      vi<-ix[i]:(ix[i+1])
      vj<-jx[j]:(jx[j+1])
      imgB1[vi,vj]<-.findbody(img=imgC[vi,vj],bandwidth=bw)
    }
  }

  #---enhance signal intensity of distant pixels
  imgC <- .enhanceDistSignal(imgC,nucmask)

  #---apply a filter to smooth signal based on local intensity
  #---this filter is only used to find picks
  imgF <- .nfilter(imgC,100)
  #---find picks
  imgB2 <- array(0,dim(imgC))
  for(i in 1:nwin){
    for(j in 1:nwin){
      vi<-ix[i]:(ix[i+1])
      vj<-jx[j]:(jx[j+1])
      imgB2[vi,vj]<- .findpicks(img=imgC[vi,vj],imgf=imgF[vi,vj], fpicks=fpicks)
    }
  }
  imgB <- imgB1 | imgB2
  return(imgB)
}

###################################################
### kernel adjusts
###################################################
#find body (one-signal assumption)
.findbody <- function(img,bandwidth){
  dt <- as.numeric(img)
  qt <- quantile(dt,probs=seq(0, 1, 0.01))
  img>qt[101-bandwidth]
}

#find picks
.findpicks <- function(img,imgf,fpicks=0.01){
  img > imgf+fpicks
}

#find edges
.findedges <- function(img){
  distmap(!img)==1
}

#fft on img
.nfilter <- function(img,df=50){
  df <- as.integer(dim(img)/df)
  idx <- df%%2==0
  df[idx] <- df[idx]+1
  filter <- array(1,df)
  filter <- filter/sum(filter)
  dx <- dim(img)
  df <- dim(filter)
  cx <- dx%/%2
  cf <- df%/%2
  wf <- matrix(0, nrow = dx[1], ncol = dx[2])
  wf[(cx[1] - cf[1]):(cx[1] + cf[1]), (cx[2] - cf[2]):(cx[2] + cf[2])] <- filter
  wf <- fft(wf)
  dim(img) <- c(dx[1:2], prod(dx)/prod(dx[1:2]))
  index1 <- c(cx[1]:dx[1], 1:(cx[1] - 1))
  index2 <- c(cx[2]:dx[2], 1:(cx[2] - 1))
  pdx <- prod(dim(img)[1:2])
  imgf <- apply(img, 3, function(xx) {
    dim(xx) <- dx[1:2]
    Re(fft(fft(xx) * wf, inverse = TRUE)/pdx)[index1, index2]
  })
  dim(imgf) <- dx
  return(imgf)
}

# enhance signal intensity of distant pixels
.enhanceDistSignal <- function(imgG, nucmask){
  #---for every pixel, compute the distance to the closest cell body
  dm <- EBImage::distmap(nucmask==0)
  dm <- dm/max(dm)
  #---based on dm matrix, enhance signal intensity of distant pixels
  idx <- imgG>quantile(imgG,probs=0.5)
  imgG[idx] <- (imgG*(1+dm^2))[idx]
  return(imgG)
}

