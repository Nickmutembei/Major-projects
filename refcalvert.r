# refcal.r
#
# Date of program : Nov 05
# R file to read ref plane photograph and store it as Rdata for
# using in surfcolour.r
# also calculates deltapix0 and f0
# This calibration data is written to disk as a text file. It will be needed in 
# surfcolour.r
#
#
#  first clear workspace,close all graphics windows 
rm(list=ls())
times <- Sys.time()
graphics.off()
#
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
		image(xplot,yplot,fliplr(t(aplot)),asp=1,xlab="",ylab="")
	}
####################  find_end  ##################################
# finds the location of the first frame of data points less than a set level
	find_end <- function(x, frame=1,level=0) {
		n <- length(x)
		for (i in 1:(n-frame+1)) {
			if(all(x[i:(i+frame-1)]< level)) {return(i)}
		}
		return(NULL)
	}
##################################################################
####################  find_end2  ##################################
# finds the location of the first frame of data points greater than a set level
	find_end2 <- function(x, frame=1,level=0) {
		n <- length(x)
		for (i in 1:(n-frame+1)) {
			if(all(x[i:(i+frame-1)]> level)) {return(i)}
		}
		return(NULL)
	}
#####################################################################
############## start main program ###################################
#
	library(jpeg) # thresholding and maxImg
# read ref plane
	grgb0 <- readJPEG("refplane.jpg")
	gr0 <- grgb0[,,2] # use green channel for FTP

	imh0 <- nrow(gr0)
	imw0 <- ncol(gr0)

# get calibration data
	
# get first estimate of horizontal location of left and right of circle
	temp1 <- apply(grgb0[(imh0/2-7):(imh0/2+7),,1],2,sum)
	v1raw <- which.min(temp1[1:round(imw0/2)])
	v2raw <- which.min(temp1[(round(imw0/2)+1):imh0])+round(imw0/2)
	av <- sum(temp1[(v1raw-20):(v1raw+20)])/41
	av2 <- sum(temp1[(v2raw-20):(v2raw+20)])/41
# find width of line and take mid-pixel as better estimate of centre of line
	tempf <- ifelse(temp1 < av,1,0)
	v1l <- find_end2(tempf[(v1raw-30):(v1raw+30)],4,0.5)+v1raw-30-1
	v1r <- v1raw+30+1-find_end2(tempf[(v1raw+30):(v1raw-30)],4,0.5)
	tempf <- ifelse(temp1 < av2,1,0)
	v2l <- find_end2(tempf[(v2raw-30):(v2raw+30)],4,0.5)+v2raw-30-1
	v2r <- v2raw+30+1-find_end2(tempf[(v2raw+30):(v2raw-30)],4,0.5)
	v1h <- (v1r+v1l)/2
	v2h <- (v2r+v2l)/2
	
	#deltapix0 <- 2*393/(v2r+v2l-v1r-v1l) #Nuffield system
	deltapix0 <- 2*316/(v2r+v2l-v1r-v1l)
	
# translate v1 and v2 into vertical readings for getting f0
	v1 <- round((imh0-v2h+v1h)/2)
	v2 <- round((imh0+v2h-v1h)/2)
	temp <- fft(gr0[v1:v2,round(imw0/2)]) #just look at central column
	klim <- round((v2-v1)/2)
	k <- which.max(abs(temp[10:klim]))
	k <- 9+k-1
	# do a barycentric fit to get better estimate of f0
	pminus <- (abs(temp[k-1]))**2
	pplus <- (abs(temp[k+1]))**2
	k0 <- k + (pplus-pminus)/(2*(pplus+pminus))
	f0 <- k0/((v2-v1)*deltapix0)
	cat(sprintf('Raw results %g %g %g %g\n',v1,v2,deltapix0,f0))

# write an R data file with deltapix0, f0 and gr0
save(deltapix0,f0,gr0,file="refdata.Rdata")
time6 <- Sys.time()

# write a text file with deltapix0 and f0 that can be read by SuperCard
# not needed for the analysis so could be deleted
	zz <- file("calib.txt","wt")
	outputdata <- c(round(deltapix0,5),round(f0,5))
	write(outputdata,zz)
	close(zz)

# finish up
	cat(sprintf('Run completed successfully\n'))
#END
	
