# blob.r
# program to do the analysis of the UV photo when threshold is changed
# or width of central analysed part is changed
# Written July 05

#  first clear workspace,close all graphics windows 
rm(list=ls())
times <- Sys.time()
#graphics.off()

##### INCLUDE FUNCTIONS FIRST  #########
########################  imcrop  ##################################
imcrop <- function(p,toplim,leftlim,newheight,newwidth) {
#imcrop - to crop an image matrix
#toplim = first row to be in new cropped file
#leftlim = first column to be in new cropped file
#newheight = number of rows in new cropped file
#newwidth = number of columns in new cropped file
#
	bottomlim <- toplim+newheight-1
	rightlim <- leftlim+newwidth-1
	q <- p[toplim:bottomlim,leftlim:rightlim]
	return(q)
}
############################################################################
##########################   fliplr   ######################################
fliplr <- function(p) {
	if(is.vector(p)) {  #check if is a vector, if so treat as single row
		nc <- length(p)
		q <- p[nc:1]
	}
		if (is.matrix(p)) { #if a matrix, reverse each row
			nc <- ncol(p)
			q <- p[,nc:1]
		}
			return(q)
}
############################################################################
##########################   flipud   ######################################
flipud <- function(p) {
	if(is.vector(p)) {  #check if is a vector, if so treat as single row
		nc <- length(p)
		q <- p[nc:1]
	}
		if (is.matrix(p)) { #if a matrix, reverse each column
			nr <- nrow(p)
			q <- p[nr:1,]
		}
			return(q)
}
##########################################################################
###################### bwlabel2 ##################################
	bwlabel2 <- function(imgin){
		w <- dim(imgin)[2]
		h <- dim(imgin)[1]
		matrix(.C("bwlabel2",as.double(imgin),as.integer(w),as.integer(h),imgout=double(w*h))$imgout,nrow=h,ncol=w)
	}
###########################  myimage  #####################################
################################################################################
	myimage <- function(a,xpixel=1,ypixel=1){
		nrows0 <- nrow(a)
		ncols0 <- ncol(a)
		fact <- round(nrows0/120)
		nrp.seq <- seq(1,nrows0,fact)
		ncp.seq <- seq(1,ncols0,fact)
		aplot <- a[nrp.seq,ncp.seq]
		nrows <- nrow(aplot)
		ncols <- ncol(aplot)
		xplot <- seq(1,ncols)*xpixel*fact
		yplot <- seq(1,nrows)*ypixel*fact
		maxaplot <- max(aplot,na.rm=TRUE)
		maxaplot <- round(maxaplot/10)*10
		image(xplot,yplot,fliplr(t(aplot)),asp=1,,xlab="",ylab="")
	}
########################################################################
######################### start main program  ##############################
# dynamic library loading for laptop
#  dyn.load("~/desktop/c_for_backanal/blobs/build/blobs.dylib",local=FALSE)
	
# dynamic library loading for G4
	dyn.load("~/xcode_files/blobs/build/blobs.dylib",local=FALSE)
	
	library(rimage,lib.loc="/Users/fiona/Library/R/library") # thresholding and maxImg
	load("/Users/fiona/desktop/supercardr/mycol.Rdata") # MATLAB colours in contour plots
	time1 <- Sys.time()
	
# read text file with size of gruvred binary file
	temp <- read.table("~/Desktop/SuperCardR/rdata/uvdata.txt")
	imh <- temp$V3
	imw <- temp$V4
	len <- imh*imw
	flt <- c(1:temp$V5)*0.0+1 #create filter coefficients, length fringepix
	time2 <- Sys.time()
	
# binary read of uv file
	zz <- file("~/Desktop/SuperCardR/rdata/gruvred.raw","rb")
	temp <- readBin(zz,integer(),n=len,size=1,signed=FALSE)
	close(zz)
	gruvred <- matrix(temp,imh,imw)
	time3 <- Sys.time()
	
# read text file with threshold and width data from SuperCard
	temp <- read.table("~/Desktop/SuperCardR/rdata/threshwidth.txt")
	thresh <- temp$V1
	width <- temp$V2
	start0 <- (1-width)/2
	time4 <- Sys.time()
	
#crop to centre section to reduce processing
	newimw <- round(imw*width) 
	startpix <- round(imw*start0)  
	uvreduced2 <- imcrop(gruvred,1,startpix,imh,newimw)
	uvreduced <- filter(uvreduced2,flt,method="convolution",sides=2,circular=TRUE)/length(flt)
 	test1uv <- thresholding(uvreduced,mode="fixed",thresh*255)  
	time5 <- Sys.time()
	
#now locate blobs
	test3uv <- bwlabel2(test1uv)
	num <- max(test3uv)
	time6 <- Sys.time()
	
#eliminate blobs of less than 40 pixels
	npixels <- matrix(0,1,num)
	rbar <- matrix(0,1,num)
	cbar <- matrix(0,1,num)
	for (j in 1:num){
		temp <- (which(test3uv==j))/imh
		c <- ceiling(temp)
		r <- round((temp+1-c)*imh) 
		npixels[j] <- length(c) 
		rbar[j] <- sum(r)/npixels[j] 
		cbar[j] <- sum(c)/npixels[j] 
	}
	skiplog <- matrix(0,1,num)
	for (j in 1:num){
		skiplog[j] <- 0 
		if (npixels[j] > 40){
			skiplog[j] <- 1 
		}
	}
	num2 <- sum(skiplog) #number of blobs found
	time7 <- Sys.time()
	
# find pixel locations only of blobs > 40 pixels		
	offset1 <- 0 
	cbar2 <- matrix(0,1,num2) 
	rbar2 <- matrix(0,1,num2) 
	for (j in 1:num){
		if (skiplog[j]){
			rbar2[j-offset1] <- rbar[j] 
			cbar2[j-offset1] <- cbar[j] 
		}
		else {
			offset1 <- offset1+1 
		}
	} 
# correct so pixel locations are relative to full cropped image not
# narrower UV analysed image or original uncropped image
	cbar2 <- cbar2+startpix-1

# write a raw data file to disk containing the thresholded blobs
# for reading and display in SuperCard 
# transpose in writeBin gets the image into portrait format
	zz <- file("~/Desktop/SuperCardR/rdata/blobs.raw","wb")
	writeBin(as.integer(as.vector(t(test1uv*255))),zz,size=1)
	close(zz)
	time8 <- Sys.time()

# write a text file containing the size of the blob file
# so SuperCard knows the size of the file to be read
	zz <- file("~/desktop/supercardr/rdata/blobdata.txt","wt")
	outputdata <- c(imh,newimw,num2,thresh) 
	write(outputdata,zz)
	close(zz)
	time9 <- Sys.time()

# write a text file containing pixel locations that will be needed for clinical analysis
# file has to be same format as one that could be written from SuperCard
	zz <- file("~/desktop/supercardr/rdata/pixels.txt","wt")
	outputdata <- cbind(t(cbar2),t(rbar2))
	write.table(outputdata,zz,col.names=FALSE,row.names=FALSE)
	close(zz)
	time10 <- Sys.time()
	
# finish up
	deltat <- (time10-times)
	cat(sprintf('Processing time is %g\n',deltat))
	cat(sprintf('Run completed successfully\n'))
#END


