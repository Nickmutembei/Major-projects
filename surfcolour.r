#  surfcolour.r
#  
#  Date of program : Nov 2005
# Edited May 2009 for use on Windows
#
#  Sept 06 changed the tolerance in the unwrapping to 3.6 instead of 4
#
#  R file to do analysis of set of back photographs
#  Called from SuperCard (can be run standalone in R)
#  First it autocrops to back region, creates a mask, does FTP to give 3D back surface.
#  Then it processes UV image and writes the blob raw data to disk for SuperCard to display.
#
#  Object illuminated from above so image with horizontal fringes.
#  Camera mounted in portrait orientation.
#  h is the height of the surface, phase is the phase angle of the surface
#  d0 is distance between projector and camera
#  l0 is distance between camera and reference plane
#  f0 is grating frequency (i.e. lines/mm) at centre of ref plane
#  includes Cobb calculation on estimated vertebrae positions
#
#  Changed width of filter from plus/minus 0.4 to 0.7 July 2007
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
############################  meshgrid  ######################################
###############################################################################
meshgrid <- function(a,b){
	list(
		 x=outer(b*0,a,FUN="+"),
		 y=outer(b,a*0,FUN="+")
		 )
}
#########################  meshplot  ########################################
##############################################################################
#does a persp plot with about 60 points in the rows direction
meshplot <- function(a){
	nrows0 <- nrow(a)
	ncols0 <- ncol(a)
	fact <- round(nrows0/60)
	nrp.seq <- seq(1,nrows0,fact)
	ncp.seq <- seq(1,ncols0,fact)
	aplot <- a[nrp.seq,ncp.seq]
	nrows <- nrow(aplot)
	ncols <- ncol(aplot)
	xplot <- seq(1,ncols)
	yplot <- seq(1,nrows)
	persp(xplot,yplot,fliplr(t(aplot)),theta=0,phi=40,expand=0.5,asp=1,axes=FALSE,box=FALSE)
}
###########################  mycontourplot  #####################################
################################################################################
#mycontourplot <- function(a,xpixel=1,ypixel=1){
	mycontourplot <- function(a,col1=mycol,xpixel=1,ypixel=1){
# make sure mycol is loaded in main program if you want blue(lowest) to red(highest)
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
#        image(xplot,yplot,fliplr(t(aplot)),asp=1,,xlab="",ylab="",axes=FALSE)
		image(xplot,yplot,fliplr(t(aplot)),asp=1,,xlab="",ylab="",axes=FALSE,col=col1)
		contour(xplot,yplot,fliplr(t(aplot)),levels=seq((maxaplot-100),maxaplot,by=5),asp=1,add=TRUE,axes=FALSE)
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
###########################  bamplot  #####################################
################################################################################
	bamplot <- function(a,col1=mycol,limit=c(0,50)){
# make sure mycol is loaded in main program if you want blue(lowest) to red(highest)
		nrows0 <- nrow(a)
		ncols0 <- ncol(a)
		fact <- round(nrows0/120)
		nrp.seq <- seq(1,nrows0,fact)
		ncp.seq <- seq(1,ncols0,fact)
		aplot <- a[nrp.seq,ncp.seq]
		nrows <- nrow(aplot)
		ncols <- ncol(aplot)
		xplot <- seq(1,ncols)*fact
		yplot <- seq(1,nrows)*fact
		maxaplot <- max(aplot,na.rm=TRUE)
		maxaplot <- round(maxaplot/10)*10
		
		filled.contour(xplot,yplot,fliplr(t(aplot)),asp=1,zlim=limit,xlab="",ylab="",axes=FALSE,col=col1)
	}
#################################################################
#################### mypoly  ####################################
	mypoly <- function(p,myx,myy){
		yy <- p[1]+p[2]*myx+p[3]*myx^2+p[4]*myx^3+p[5]*myx^4+p[6]*myx^5
		ss <- sum((yy-myy)^2)
		return(ss)
	}
###################################################################
################## polyval  #######################################
	polyval <- function(p,x){
#use Horner's method which needs descending coefficients'
		nc <-length(p)
		p <- fliplr(p) #make coefficients descending
		y <-seq(1,length(x))*0
		if (length(p)>0){
			y[] <- p[1]
		}
		for (ii in 2:nc){
			y <- x*y+p[ii]
		}
		return(y)
	}
#######################################################################
#####################  excise  #########################################
	excise <- function(maskin){
		w <- dim(maskin)[2]
		h <- dim(maskin)[1]
		matrix(.C("excise",as.double(maskin),as.integer(w),as.integer(h),maskout=double(w*h))$maskout,nrow=h,ncol=w)
	}
##########################################################################
###################### bwlabel2 ##################################
	bwlabel2 <- function(imgin){
		w <- dim(imgin)[2]
		h <- dim(imgin)[1]
		matrix(.C("bwlabel2",as.double(imgin),as.integer(w),as.integer(h),imgout=double(w*h))$imgout,nrow=h,ncol=w)
	}
######################### hanning ######################################
#######################################################################
	hanning <- function(n){
#produces a SYMMETRIC hanning window without the first and last zero values
		if (n%%2){
#odd length window
			half <- (n+1)/2
			w <- 0.5*(1-cos(2*pi*(1:half)/(n+1)))
			temp <- fliplr(w[1:(length(w)-1)])
			w <- c(w,temp)
		}
		else {
			half <- n/2
			w <- 0.5*(1-cos(2*pi*(1:half)/(n+1)))
			temp <- fliplr(w)
			w <- c(w,temp)
		}
	}
###############################################################
################### rgl plot ##################################
	backrglplot <- function(surf){
		sizeh <- nrow(surf)
		sizew <- ncol(surf)
		nr.seq <- seq(1,sizeh,2)
		nc.seq <- seq(sizew,1,-2) #this flips surface so left is at the left
		htemp <- surf[nr.seq,nc.seq]
		colorlut <- heat.colors(256)
		col <- colorlut[256-255*(htemp-min(htemp))/(max(htemp)-min(htemp))]
		rgl.surface(1:length(nr.seq),1:length(nc.seq),htemp/2,col)
	}# rotate 3D function, time demanding, not neccesary
##################################################################
####################  find_end  ##################################
# finds the location of the first frame of data points less than a set level
	find_end <- function(x, frame=1,level=0) {
		n <- length(x)
		for (i in 1:(n-frame+1)) {
			if(all(x[i:(i+frame-1)]< level)) {return(i)}
		}
		return(0)
	}
##################################################################
####################  find_end2  ##################################
# finds the location of the first frame of data points greater than a set level
	find_end2 <- function(x, frame=1,level=0) {
		n <- length(x)
		for (i in 1:(n-frame+1)) {
			if(all(x[i:(i+frame-1)]> level)) {return(i)}
		}
		return(0)
	}
############################################################################
###################### unwrap   ##############################################
unwrap<-function(p,cutoff){
#unwraps a phase function p, using discontinuity tolerance cutoff
# check for vector
vecflag <- 0
if (is.vector(p))  {
    p <- (as.matrix(p)) # produces a column matrix, nrow=length(p)
    vecflag <- 1
}
nrows <- nrow(p) 
ncols <- ncol(p)
# Unwrap phase angles 
b<-matrix(0,nrow=nrows,ncol=ncols)
b[2:nrows,] <- diff(p)    # Differentiate phases but make b equal in length to p
c <- -(b > cutoff) 
d <- (b < -cutoff)		        # Locations of jumps.
e <- (c + d) * 2 * pi;          # Array of 2*pi jumps.
f <- apply(e,2,cumsum)			# Integrate to get corrections.
# take initial jump (may be more than 2pi) into account
jump <- round(b[2,]/(2*pi))
jump1 <- ifelse(jump<2,0,jump-1)
jump2 <- ifelse(jump>-2,0,jump+1)
jump3 <- (jump1+jump2)*-2*pi
jump4 <- matrix(jump3,nrow=nrows,ncol=ncols,byrow=TRUE)
jump4[1,] <- 0.0
q <- p + f + jump4;    # Reshape output if originally transposed
if (vecflag) { q <- as.vector(q)
    }
return(q)
}

######################### start main program  ##############################
#
# dynamic library loading for G4
	dyn.load("OxtersDLL/Debug/oxtersDLL.dll",local=FALSE)
	dyn.load("blobsDLL/Debug/blobsDLL.dll",local=FALSE)
#	dyn.load("~/xcode_files/unwrapc/build/unwrapc.dylib",local=FALSE)
	
	library("rimage") # jpeg read, thresholding and maxImg
	load("mycol.Rdata") # MATLAB colours in contour plots
	
#	l0 <- 2000 #Nuffield system
#	d0 <- 505
	l0 <- 2200 #ROH system
	d0 <- 468
	theta <- atan(d0/l0)
	
	time1 <- Sys.time()
	
#read ref plane and calibration data to get gr0,deltapix0 and f0
# written to R data file earlier by refcal.r
	load("refdata.Rdata")

# read data file
	grgb <- read.jpeg("back.jpg")
	gr <- grgb[,,2] #green channel
	time2 <- Sys.time()
	 
# get cropping limits
	imh0 <- nrow(gr)
	imw0 <- ncol(gr)
	centreh0 <- round(imh0/2)
	centrew0 <- round(imw0/2)
	fringepix <- round(1/(deltapix0*f0)) #number of pixels in fringe at centre
#	flt <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
	flt <- c(1:fringepix)*0.0+1 #create filter coefficients, length fringepix
	# filter whole data file down the columns and threshold (binarize)
	tmp1 <- filter(gr,flt,method="convolution",sides=2,circular=TRUE)/length(flt)
	mask0 <- thresholding(tmp1, mode="fixed", th=0.24)
	#	
	tmp <- mask0[,centrew0]
	toplim1 <- centreh0+50+1-which.min(tmp[(centreh0+50):1])
	bottomlim1 <- centreh0+50-1+which.min(tmp[(centreh0+50):imh0])
#	rlimlev <- round(toplim+0.135*(bottomlim-toplim)) # for dummy
	rlimlev <- round(toplim1+0.25*(bottomlim1-toplim1)) # for people
	tmp2 <- mask0[rlimlev,]
	leftlim1 <- centrew0-which.min(tmp2[(centrew0-1):1]) # searching to left for first < limlim
	rightlim1 <- centrew0+which.min(tmp2[centrew0:imw0])-1
	
	# redo top and bottom limits at centre column of back (rather than photo)
	# this corrects limits if large imbalance means neck of patient is offcentre
	tmp3 <- mask0[,round((leftlim1+rightlim1)/2)]
	toplim <- centreh0+50+1-which.min(tmp3[(centreh0+50):1])
	bottomlim <- centreh0+50-1+which.min(tmp3[(centreh0+50):imh0])
#	rlimlev <- round(toplim+0.135*(bottomlim-toplim)) # for dummy
	rlimlev <- round(toplim+0.25*(bottomlim-toplim)) # for people
	tmp4 <- mask0[rlimlev,]
	leftlim <- centrew0-which.min(tmp4[(centrew0-1):1]) # searching to left for first < limlim
	rightlim <- centrew0+which.min(tmp4[centrew0:imw0])-1
	time3 <- Sys.time()

# select suitable width for FFT
	deltan0 <- bottomlim-toplim 
	fftwidth <- c(512, 640, 768, 896, 960, 1024, 1152, 1280, 1408, 1536, 1664) 
	deltan1 <- deltan0+10-fftwidth 
	negposcode <- ifelse(deltan1 > 0,1,0)
	pos <- which.min(negposcode) #first location of negative deltan1
	fftw <- fftwidth[pos] 
	imw <- rightlim-leftlim+1 
	diftb <- deltan0-fftw 
	toplim0 <- toplim+trunc(diftb/2)
	time4 <- Sys.time()

# crop image
	g <- imcrop(gr,toplim0,leftlim,fftw,imw)
	time5 <- Sys.time()
	
# crop mask
	mask <- imcrop(mask0,toplim0,leftlim,fftw,imw)
	time6 <- Sys.time()
	
# zero bits outside the back included to have a 'good' FFT length
# make sure no end effects remain from the filtering operation
	zerolim <- trunc(-diftb/2)
	if (zerolim < trunc(length(flt)/2)) zerolim <- trunc(length(flt)/2) 
	mask[1:zerolim,] <- 0.0
	mask[(fftw-zerolim+1):fftw,] <- 0.0
	time7 <- Sys.time()
		
	time8 <- Sys.time
			
# crop ref plane
	g0 <- imcrop(gr0,toplim0,leftlim,fftw,imw) 
	time9 <- Sys.time()
	
# improve mask by dilation and erosion
	#mask2 <- maxImg(mask)
	#mask4 <- maxImg(mask2)
	#mask4 <- maxImg(mask4)
	mask4 <- mask
	mask5 <- mask4 
	time11 <- Sys.time()
			
# now find positions of oxters and excise arms
# first find oxters
# crop to just the outer thirds (where oxters will be)
	subh <- round(fftw/4)
	subw <- round(imw/3)
	left <- imcrop(mask4,subh,1,(2*subh),subw)
	right <- imcrop(mask4,subh,imw-subw+1,(2*subh),subw)
	test1 <- diff(t(left))
	test1a <- diff(t(right)) 
	test4 <- apply((abs(test1)),2,sum) 
	test5 <- apply((abs(test1a)),2,sum)  
# get the height of the left oxter
#	posloxy <- subh-1+find_end2(test4,10,1) # searching down for first > 1
	ytemp <- find_end(fliplr(test4),2,1) # searching up for first < 1
	if (ytemp==0) ytemp <- length(test4)
	posloxy <- 3 *subh-ytemp+1
#then height of right oxter
#	posroxy <- subh-1+find_end2(test5,10,1)
	ytemp <- find_end(fliplr(test5),2,1)
	if (ytemp==0) ytemp <- length(test5)
	posroxy <- 3 *subh-ytemp+1
#now horiz position of left oxter
	test6 <- mask4[posloxy,6:round(imw/2)] # avoid edges
	ytemp <- find_end(test6,1,1)
	posloxx <- ytemp+5
#now horiz position of right oxter
	test7 <- fliplr(mask4[posroxy,round(imw/2):(imw-6)])
	ytemp <- find_end(test7,1,1)
	posroxx <- imw-ytemp-5
	
# now get rid of upper arms directly across from oxters and down for 100 pixels
#first left
	mask5[posloxy:(posloxy+100),1:(posloxx+5)] <- 0.0
#then right
	mask5[posroxy:(posroxy+100),(posroxx-5):imw] <- 0.0
	
#get rid of lower arms, making sure hips stay
#first left
	tempmask <- mask4[(posloxy+101):fftw,1:(posloxx+5)]
	tempmask <- fliplr(tempmask)
	tempmask2 <- excise(tempmask)
	tempmask2 <- fliplr(tempmask2)
	mask5[(posloxy+101):fftw,1:(posloxx+5)] <- tempmask2

#then right
	tempmask <- mask4[(posroxy+101):fftw,(posroxx-5):imw]
	tempmask2 <- excise(tempmask) 
	mask5[(posroxy+101):fftw,(posroxx-5):imw] <- tempmask2 
	time12 <- Sys.time()
	
# get rid of upper arms towards shoulders, going upwards and 
# outwards at an angle of 20 or 30 degrees from the vertical from oxters
	myang <- 30/57.2957795
	
# first left
	myy <- 1:posloxy
	myx1 <- fliplr(round(posloxx-myy*tan(myang)))
	myx <- ifelse(myx1 > 0, myx1,1)
	
	for (i in 1:posloxy){
		mask5[i,1:myx[i]] <- 0.0
	}

# then right
	myy <- 1:posroxy
	myx1 <- fliplr(round(posroxx+myy*tan(myang)))
	myx <- ifelse(myx1 > imw, imw,myx1)
	
	for (i in 1:posroxy){
		mask5[i,myx[i]:imw] <- 0.0
	}
	
	
#now do FTP processing
#get size of new arrays
	imh <- fftw  #height of image (number of rows), imw already defined
#cat(sprintf('Doing FFT processing\n'))
	centrew <- trunc(imw/2) 
	centreh <- trunc(imh/2) 
#do FFT for ref plane
	fref <- mvfft(g0) 
	compfref <- apply(Mod(fref),1,max) 
	k1 <- which.max(compfref[10:(imh/2)]) 
	k1 <- k1+9 #compensate for position of peak being relative to line 10
	time13 <- Sys.time()
		
#do FFT for object
	f <- mvfft(g)
	time14 <- Sys.time()
			
#create a comb filter, pass band around fundamental
#zero elsewhere
	filtfjb <- matrix(0,imh,imw)  #filt for back
#calculate filter bandwidth 
	bw <- round(0.7*k1) #half bandwidth for filter
	klow <- k1-bw #k1 is position of fringe frequency in spec
	khigh <- k1+bw 
#apply shaped filter
	filtwind <- hanning(2*bw+1)
	for (k in klow:khigh) { 
			filtfjb[k,] <- filtwind[k-klow+1] 
	} 
	time15 <- Sys.time()
								
#apply filter to back and ref plane
	ffilt <- f*filtfjb 
	ffiltref <- fref*filtfjb 
	time16 <- Sys.time()
						
#now do IFT for both object and ref plane
	invf <- mvfft(ffilt,inverse=TRUE) 
	invf0 <- mvfft(ffiltref,inverse=TRUE) 
	time17 <- Sys.time()
						
#calculate phase
	phase <- -Im(log(invf*Conj(invf0))) 
	time18 <- Sys.time()
						
#unwrap
	centreh <- round(fftw/2)
	centrew <- round(imw/2)

#we will unwrap across the rows from the centre out,
#but first unwrap central vertical column from 
#centre point up and from centre point down
	line1c <- flipud(phase[1:centreh,centrew]) 
	line2c <- phase[centreh:imh,centrew] 
#unwrap upwards and downwards from same central point
	line3c <- unwrap(line1c,3.6) 
	line6c <- unwrap(line2c,3.6) 
	line4c <- flipud(line3c) 
#remove centre point from this half column (top half)
	line5c <- line4c[1:(centreh-1)] 
	line7c <- c(line5c, line6c) 
#set central line of phase array to this unwrapped line
	phase[,centrew] <- line7c 

#now split matrix into the two halves, each including central line
# flip the left half and transpose (unwrap works down the columns)
# transpose right; unwrap both; transpose both; flip left back
# finally, re-combine to get whole unwrapped matrix
	phasel <- t(phase[,centrew:1])
	phaser <- t(phase[,centrew:imw])
	uphasel <- fliplr(t(unwrap(phasel,3.6)))
	uphaser <- t(unwrap(phaser,3.6))
	uphase <- cbind(uphasel,uphaser[,2:ncol(uphaser)])
#	uphase <- unwrapc(phase,4)
	time19 <- Sys.time()
						
#adjust phase so min value corresponds roughly with reality
#so that calculation of h is more accurate
	toptemp <- round(imh/8) 
	bottomtemp <- round(imh/2) 
	lefttemp <- round(imw*3/16) 
	righttemp <- round(imw*13/16) 
	tempphase <- uphase[toptemp:bottomtemp,lefttemp:righttemp] 
#figure,mesh(tempphase) 
	mintemp <- min(tempphase) 
	deltatemp <- -33-mintemp # -33 corresponds to h=120 mm
	ntemp <- round(deltatemp/(2*pi)) 
	shiftphase <- ntemp*2*pi 
	uphase <- uphase+shiftphase
	time20 <- Sys.time() 
	
#calculate height from unwrapped phase
	h <- l0*uphase/(uphase-2*pi*f0*d0)
	time21 <- Sys.time()
	
# perhaps add smooth results later ####
	flt2 <- c(1:40)*0.0+1
	h1 <- filter(h,flt2,method="convolution",sides=2,circular=TRUE)/length(flt2)
	time23 <- Sys.time()
	
# apply zero mask
	h2 <- h1*mask5 #h2 has zero mask
	time22 <- Sys.time() 
	
# reduce to every second/third/fourth point for processing
	fact <- 2 
	if (fftw < 600) fact <- 1
	if (fftw > 1600) fact <- 3 
	deltapix <- deltapix0*fact 
	nr.seq <- seq(1,imh,fact)
	nc.seq <- seq(1,imw,fact)
	h3 <- h2[nr.seq,nc.seq] #h3 reduced dataset for processing,filtered, zero mask
	sizeh <- nrow(h3) 
	sizew <- ncol(h3) 
	x1 <- seq(deltapix,deltapix*sizew,deltapix)
	y1 <- seq(deltapix,deltapix*sizeh,deltapix)
	time24 <- Sys.time()
		
#reduce mask too and make it have NaNs off the area of interest
	maskred <- mask5[nr.seq,nc.seq]
#	maskred <- ifelse (maskred < 0.1, NaN, maskred) #mask NaNs if maskred < 0.1
	maskred[maskred < 0.1] <- NaN
	time25 <- Sys.time()
													
#now find marker positions
# get suitable mix of channels to emphasise markers, minimise fringes	
	#gruv <- grgbuv[,,1]-0.3*(grgbuv[,,3]) # R and B channels for UV
	gruv <- grgb[,,3]-0.9*(grgb[,,1]) # B and R channels for turquoise/blue
	#gruv <- grgb[,,2]-0.9*(grgb[,,1]) # G and R channels for green
	#gruv <- grgb[,,2]-1.1*(grgb[,,3]) # G and B channels for yellow
	time26 <- Sys.time()
	
#crop marker image to same as final back image size
	gruvred <- imcrop(gruv,toplim0,leftlim,fftw,imw) 
	gruvred[gruvred < 0] <- 0
	
#crop further to centre section to reduce processing
	newimw <- round(imw*0.5) # medium width
	startpix <- round(imw*0.25)  
#	uvreduced <- imcrop(gruv,toplim0,leftlim+startpix,imh,newimw)
	uvreduced2 <- gruvred[,startpix:(startpix+newimw-1)]
	uvreduced <- filter(uvreduced2,flt,method="convolution",sides=2,circular=TRUE)/length(flt)
	time27 <- Sys.time() 

#now locate blobs
	thresh <- 0.1
	test1uv <- thresholding(uvreduced,mode="fixed",thresh) # standard threshold  
	if (sum(test1uv) > 6000) {
		thresh <- 0.2
		test1uv <- thresholding(uvreduced,mode="fixed",thresh)
		if (sum(test1uv) > 6000) {
			thresh <- 0.3
			test1uv <- thresholding(uvreduced,mode="fixed",thresh)
			if (sum(test1uv) > 6000) {
				thresh <- 0.4
				test1uv <- thresholding(uvreduced,mode="fixed",thresh)
			}
		}
	}
	test3uv <- bwlabel2(test1uv)
	num <- max(test3uv)
	time28 <- Sys.time()
	
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
	time28b <- Sys.time()

# write a raw data file to disk containing the thresholded blobs
# for reading and display in SuperCard 
# transpose in writeBin gets the image into portrait format
	zz <- file("blobs.raw","wb")
	writeBin(as.integer(as.vector(t(test1uv*255))),zz,size=1)
	close(zz)
	time29 <- Sys.time()

# write a text file containing the size of the blob file
# so SuperCard knows the size of the file to be read
# also the number of blobs detected and what threshold was
	zz <- file("blobdata.txt","wt")
	outputdata <- c(imh,newimw,num2,thresh)
	write(outputdata,zz)
	close(zz)
	time30 <- Sys.time()

# write a raw data file to disk containing the cropped uv image
# for reading and further analysis in blob.r if necessary 
	zz <- file("gruvred.raw","wb")
	writeBin(as.integer(as.vector(gruvred*255)),zz,size=1)
	close(zz)
	time31 <- Sys.time()
			
# write a text file containing the size of the uv cropped file
# for reading and further analysis in blob.r if necessary
	zz <- file("uvdata.txt","wt")
	outputdata <- c(leftlim,toplim0,imh,imw,length(flt))
	write(outputdata,zz)
	close(zz)
	time32 <- Sys.time()
	
# save data about back surface that will be needed for clinical analysis
	save(h3,x1,y1,maskred,imw,imh,fact,deltapix,sizew,sizeh,leftlim,toplim0,file="surfdata.Rdata")
	time33 <- Sys.time()
	
# write a text file containing pixel locations that will be needed for clinical analysis
# file has to be same format as one that could be written from SuperCard
	zz <- file("pixels.txt","wt")
	outputdata <- cbind(t(cbar2),t(rbar2))
	write.table(outputdata,zz,col.names=FALSE,row.names=FALSE)
	close(zz)
	time34 <- Sys.time()

# finish up
	deltat <- (time34-times)
	cat(sprintf('Processing time is %g\n',deltat))
	cat(sprintf('Run completed successfully\n'))
#END
