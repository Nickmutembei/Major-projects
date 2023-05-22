# clinic2.r
# program to do the rotate and tilt based on good blobs data and then all
# the clinical processing and pdf generation

# 
#  first clear workspace,close all graphics windows 
rm(list=ls())
times <- Sys.time()
graphics.off()
##### INCLUDE FUNCTIONS  ###################################################
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
		if(fact < 1) fact <- 1
		nrp.seq <- seq(1,nrows0,fact)
		ncp.seq <- seq(1,ncols0,fact)
		aplot <- a[nrp.seq,ncp.seq]
		aplot2 <- ifelse(aplot>80,80,aplot)
		nrows <- nrow(aplot)
		ncols <- ncol(aplot)
		xplot <- seq(1,ncols)*xpixel*fact
		yplot <- seq(1,nrows)*ypixel*fact
		maxaplot <- max(aplot,na.rm=TRUE)
		maxaplot <- round(maxaplot/10)*10
		image(xplot,yplot,fliplr(t(aplot2)),asp=1,zlim=c(-80,80),xlab="",ylab="",axes=FALSE,col=col1)
#		contour(xplot,yplot,fliplr(t(aplot)),levels=seq((maxaplot-100),maxaplot,by=5),asp=1,add=TRUE,axes=FALSE)
		contour(xplot,yplot,fliplr(t(aplot)),levels=seq(-60,maxaplot,by=5),asp=1,add=TRUE,axes=FALSE)
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
##################################################################
######################### start main program  ##############################
#	load some files from disc
	load("mycol.Rdata") # MATLAB colours in contour plots	
	load("surfdata.Rdata") # all data stored by surfcolour.r 
	# the following parameters are loaded  h3,x1,y1,maskred,imw,imh,fact,deltapix,sizew,sizeh,leftlim,toplim0
	# h3 = back surface heights without NaNs reduced by factor fact on original cropped image size
	# x1, y1 = vectors of x and y pixel spacings in mm (after fact taken into account) ~1 mm each
	# maskred = reduced size mask with NaNs instead of zeros
	# fact = factor used to reduce cropped image size to data for further processing
	# deltapix = pixel size in mm (after reduction by factor fact i.e. increased by fact)
	# sizew, sizeh = size of reduced data files (number of columns, number of rows)
	# leftlim, toplim0 = column number and row number for original cropping based on back location (not needed)
#	read text file containing pixel locations of stickers relative to cropped non-reduced image 
#	pixels.txt may have been written by surfcolour.r or by SuperCard
	temp <- read.table("pixels.txt") # matrix of sticker columns and rows
# make sure they are in correct order VP first then all other spinal markers, then PSIS points
	nblobs <- dim(temp)[1]
	o <- order(temp$V2)
	spine <- cbind(temp$V1[o], temp$V2[o]) 
	yspine <- spine[1:(nblobs-2),] # markers on spinous processes, first is VP, exact locations fractions of pixels
	psisl <- spine[(nblobs-1),] # left psis point
	psisr <- spine[nblobs,]     # right psis point
	if (psisr[1] < psisl[1]) {  # swop over if psis points have been stored r then l, rather than l then r
		temp2 <- psisl
		psisl <- psisr
		psisr <- temp2
	}
#   calculate integer pixel numbers for markers in reduced dataset
	dvlx <- round(psisl[1]/fact) #closest pixel column for reduced dataset
	dvly <- round(psisl[2]/fact) #closest pixel row
	dvrx <- round(psisr[1]/fact) 
	dvry <- round(psisr[2]/fact) 
	vpx <- round(yspine[1,1]/fact) 
	vpy <- round(yspine[1,2]/fact) 
	sacx <- round((dvlx+dvrx)/2) 
	sacy <- round((dvly+dvry)/2) 
	yspinepix <- round(yspine/fact) 
#yspinepix[1] will be vertebra prominens
#set up x and y grids (matrices)
	ggrid <- meshgrid(x1,y1)
	x1g <- ggrid[[1]]
	y1g <- ggrid[[2]]
#########################################################################
## now do tilt and rotate  ##############################################
	nh0 <- h3
#Step 1:get x,y coordinates in pixels and change to mm
#cat(sprintf('Doing rotate and tilt\n')) 
#cat(sprintf('Step 1 - Set x,y coordinates of VP and DVs and change to mm')) 
#cat(sprintf('\n')) 
#cat(sprintf('\n')) 
# change to mm ready for rotation
	dvl <- c(dvlx*deltapix,dvly*deltapix,nh0[dvly,dvlx]) #all in mm
	dvr <- c(dvrx*deltapix,dvry*deltapix,nh0[dvry,dvrx]) 
	vp <- c(vpx*deltapix,vpy*deltapix,nh0[vpy,vpx]) 
	sac <- c((dvl[1]+dvr[1])/2,(dvl[2]+dvr[2])/2,nh0[sacy,sacx]) 
#cat(sprintf('DVL1 is %7.2f %7.2f %7.2f\n',dvl[1],dvl[2],dvl[3])) 
#cat(sprintf('DVR1 is %7.2f %7.2f %7.2f\n',dvr[1],dvr[2],dvr[3])) 
#cat(sprintf('VP1 is %7.2f %7.2f %7.2f\n',vp[1],vp[2],vp[3])) 
#cat(sprintf('SAC1 is %7.2f %7.2f %7.2f\n',sac[1],sac[2],sac[3])) 
#cat(sprintf('\n')) 
#print locations of spinal markers in mm
	numspine <- nblobs-2 
	sp1 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp1[i,] <- c(x1[yspinepix[i,1]],y1[yspinepix[i,2]],nh0[yspinepix[i,2],yspinepix[i,1]]) #x1 still 1D here
	}
#cat(sprintf('Locations of spinal markers in mm are: \n'))
#for (i in 1:numspine){
#    cat(sprintf('%8.2f %8.2f %8.2f\n',sp1[i,1],sp1[i,2],sp1[i,3])) 
#}
#cat(sprintf('\n')) 
#Step 2 translate origin to dvl
	x2 <- x1g-dvl[1] 
	y2 <- y1g-dvl[2] 
	nh2 <- nh0-dvl[3] 
	dvl2 <- c(x2[dvly,dvlx],y2[dvly,dvlx],nh2[dvly,dvlx]) 
	dvr2 <- c(x2[dvry,dvrx],y2[dvry,dvrx],nh2[dvry,dvrx]) 
	vp2 <- c(x2[vpy,vpx],y2[vpy,vpx],nh2[vpy,vpx]) 
	sac2 <- c(x2[sacy,sacx],y2[sacy,sacx],nh2[sacy,sacx]) 
	sp2 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp2[i,] <- c(x2[yspinepix[i,2],yspinepix[i,1]],y2[yspinepix[i,2],yspinepix[i,1]],nh2[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 3 rotate about y-axis to bring dvr z value to zero
	s3 <- dvr2[3]/sqrt(dvr2[1]^2+dvr2[3]^2) 
	c3 <- dvr2[1]/sqrt(dvr2[1]^2+dvr2[3]^2) 
	rot <- atan(s3/c3)*360/(2*pi) 
	x3 <- c3*x2+s3*nh2 
	y3 <- y2 
	nh3 <- -s3*x2+c3*nh2 
	dvl3 <- c(x3[dvly,dvlx],y3[dvly,dvlx],nh3[dvly,dvlx]) 
	dvr3 <- c(x3[dvry,dvrx],y3[dvry,dvrx],nh3[dvry,dvrx]) 
	vp3 <- c(x3[vpy,vpx],y3[vpy,vpx],nh3[vpy,vpx]) 
	sac3 <- c(x3[sacy,sacx],y3[sacy,sacx],nh3[sacy,sacx]) 
	sp3 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp3[i,] <- c(x3[yspinepix[i,2],yspinepix[i,1]],y3[yspinepix[i,2],yspinepix[i,1]],nh3[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 4 rotate about z-axis to bring dvrs y value to zero
	s4 <- dvr2[2]/sqrt(dvr2[1]^2+dvr2[2]^2) 
	c4 <- dvr2[1]/sqrt(dvr2[1]^2+dvr2[2]^2) 
	x4 <- x3*c4+y3*s4 
	y4 <- -x3*s4+y3*c4 
	nh4 <- nh3 
	dvl4 <- c(x4[dvly,dvlx],y4[dvly,dvlx],nh4[dvly,dvlx]) 
	dvr4 <- c(x4[dvry,dvrx],y4[dvry,dvrx],nh4[dvry,dvrx]) 
	vp4 <- c(x4[vpy,vpx],y4[vpy,vpx],nh4[vpy,vpx]) 
	sac4 <- c(x4[sacy,sacx],y4[sacy,sacx],nh4[sacy,sacx]) 
	sp4 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp4[i,] <- c(x4[yspinepix[i,2],yspinepix[i,1]],y4[yspinepix[i,2],yspinepix[i,1]],nh4[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 5 Translate x-axis to going through sacrum point
#i.e. z <- 0 plane moves to sacrum
	x5 <- x4 
	y5 <- y4 
	nh5 <- nh4-sac4[3] 
	dvl5 <- c(x5[dvly,dvlx],y5[dvly,dvlx],nh5[dvly,dvlx]) 
	dvr5 <- c(x5[dvry,dvrx],y5[dvry,dvrx],nh5[dvry,dvrx]) 
	vp5 <- c(x5[vpy,vpx],y5[vpy,vpx],nh5[vpy,vpx]) 
	sac5 <- c(x5[sacy,sacx],y5[sacy,vpy],nh5[sacy,sacx]) 
	sp5 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp5[i,] <- c(x5[yspinepix[i,2],yspinepix[i,1]],y5[yspinepix[i,2],yspinepix[i,1]],nh5[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 6 rotate about x axis until vertebra prominens 
#is over z <- 0 plane
	s6 <- -vp5[3]/sqrt(vp5[2]^2+vp5[3]^2) 
	c6 <- -vp5[2]/sqrt(vp5[2]^2+vp5[3]^2) 
	tilt <- atan(s6/c6)*360/(2*pi) 
#cat(sprintf('Tilt angle is %4.1f degrees\n',tilt)) 
	x6 <- x5 
	y6 <- y5*c6+nh5*s6 
	nh6 <- -y5*s6+nh5*c6 
	dvl6 <- c(x6[dvly,dvlx],y6[dvly,dvlx],nh6[dvly,dvlx]) 
	dvr6 <- c(x6[dvry,dvrx],y6[dvry,dvrx],nh6[dvry,dvrx]) 
	vp6 <- c(x6[vpy,vpx],y6[vpy,vpx],nh6[vpy,vpx]) 
	sac6 <- c(x6[sacy,sacx],y6[sacy,sacx],nh6[sacy,sacx]) 
	sp6 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp6[i,] <- c(x6[yspinepix[i,2],yspinepix[i,1]],y6[yspinepix[i,2],yspinepix[i,1]],nh6[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 7 rotate about z axis again until dvr is back to old y value
	s7 <- -s4 
	c7 <- c4 
	x7 <- x6*c7+y6*s7 
	y7 <- -x6*s7+y6*c7 
	nh7 <- nh6 
	dvl7 <- c(x7[dvly,dvlx],y7[dvly,dvlx],nh7[dvly,dvlx]) 
	dvr7 <- c(x7[dvry,dvrx],y7[dvry,dvrx],nh7[dvry,dvrx]) 
	vp7 <- c(x7[vpy,vpx],y7[vpy,vpx],nh7[vpy,vpx]) 
	sac7 <- c(x7[sacy,sacx],y7[sacy,sacx],nh7[sacy,sacx]) 
	sp7 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp7[i,] <- c(x7[yspinepix[i,2],yspinepix[i,1]],y7[yspinepix[i,2],yspinepix[i,1]],nh7[yspinepix[i,2],yspinepix[i,1]]) 
	}
#Step 8 shift origin to vertebra prominens
	x8 <- x7-vp7[1] 
	y8 <- y7-vp7[2] 
	nh8 <- nh7-vp7[3] 
	dvl8 <- c(x8[dvly,dvlx],y8[dvly,dvlx],nh8[dvly,dvlx]) 
	dvr8 <- c(x8[dvry,dvrx],y8[dvry,dvrx],nh8[dvry,dvrx]) 
	vp8 <- c(x8[vpy,vpx],y8[vpy,vpx],nh8[vpy,vpx]) 
	sac8 <- c(x8[sacy,sacx],y8[sacy,sacx],nh8[sacy,sacx]) 
	sp8 <- matrix(0,numspine,3)
	for (i in 1:numspine){
		sp8[i,] <- c(x8[yspinepix[i,2],yspinepix[i,1]],y8[yspinepix[i,2],yspinepix[i,1]],nh8[yspinepix[i,2],yspinepix[i,1]]) 
	}
#
# Step 9 - create new xy grid and adjust so origin is at VP
	deltax <- deltapix/cos(rot*2*pi/360) #(xmax-xmin)/(sizew-1)
	deltay <- deltapix/cos(tilt*2*pi/360) #(ymax-ymin)/(sizeh-1)
	xng <- (1:sizew)*deltax  #seq(xmin,xmax,length=sizew) 
	yng <- (1:sizeh)*deltay  #seq(ymin,ymax,length=sizeh)
	xng <- xng-xng[vpx] #shift grid zero to VP
	yng <- yng-yng[vpy]
	nggrid <-meshgrid(xng,yng)
	xng2 <- nggrid[[1]]
	yng2 <- nggrid[[2]]
	nh9 <- nh8
	x9 <- xng2
	y9 <- yng2
	xmin <- min(x9)
	ymin <- min(y9)
	vp9rownum <- vpy 
	vp9colnum <- vpx 
	sac9rownum <- sacy 
	sac9colnum <- sacx 
	dvl9 <- c(x9[dvly,dvlx],y9[dvly,dvlx],nh9[dvly,dvlx]) 
	dvr9 <- c(x9[dvry,dvrx],y9[dvry,dvrx],nh9[dvry,dvrx]) 
	vp9 <- c(x9[vp9rownum,vp9colnum],y9[vp9rownum,vp9colnum],nh9[vp9rownum,vp9colnum])
	sac9 <- c(x9[sac9rownum,sac9colnum],y9[sac9rownum,sac9colnum],nh9[sac9rownum,sac9colnum])
	nh9d <- nh9*maskred
	sp9 <- matrix(0,(numspine+1),3) #include sacrum as well as all spinous process markers
	for (i in 1:numspine){
		sp9[i,] <- c(x9[yspinepix[i,2],yspinepix[i,1]],y9[yspinepix[i,2],yspinepix[i,1]],nh9[yspinepix[i,2],yspinepix[i,1]])
	}
	sp9[numspine+1,1] <- sac9[1]
	sp9[numspine+1,2] <- sac9[2]
	sp9[numspine+1,3] <- sac9[3]
###################################################################
#cat(sprintf('Calculating clinical parameters and plotting pdf\n'))
###### START PDF PLOTTING  ########################################
###################################################################
#	oldpar <- par(no.readonly=TRUE) #save current parameters
	pdffile <- "results.pdf"
	pdf(pdffile,width=7.5,height=7.5,bg="white",onefile=TRUE,family="Times")
	layout(matrix(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,5),2,14,byrow=TRUE))
# ########HEIGHT MAP ########################
	par(mai=c(0.0,0.0,0.15,0.0))
	meshplot(nh9d)
	title("Height map")
# #########CONTOUR PLOT #######################
	par(mai=c(0.4,0.3,0.3,0.16))
	mycontourplot(nh9d)
	title("Contour plot")
# add spinal markers to the plot
for (j in 1:(nblobs-2)){
	points(yspinepix[j,1],sizeh+1-yspinepix[j,2],pch=19,col="blue")
#	points(yspinepix[j,1]+9,sizeh+1-yspinepix[j,2],pch=21,bg="blue",col="blue",cex=3)
}
dvl9rownum <- dvly 
dvr9rownum <- dvry 
dvl9colnum <- dvlx 
dvr9colnum <- dvrx
points(dvr9colnum,sizeh-dvr9rownum,pch=19,col="blue")
points(dvl9colnum,sizeh-dvl9rownum,pch=19,col="blue")
#points(dvr9colnum,sizeh-dvr9rownum,pch=21,bg="blue",col="blue",cex=3)
#points(dvl9colnum,sizeh-dvl9rownum,pch=21,bg="blue",col="blue",cex=3)
#
# add max hump markers
colliml1 <- round(vp9colnum - 0.3 * (sac9rownum - vp9rownum))
if (colliml1 < 1) colliml1 <- 1 # in case patient is very narrow for height
colliml2 <- vp9colnum
rowliml1 <- round(vp9rownum + 0.15 * (sac9rownum - vp9rownum))
rowliml2 <- round(vp9rownum + 0.55 * (sac9rownum - vp9rownum))
left <- nh9d[rowliml1:rowliml2,colliml1:colliml2]
maxleft <- max(left,na.rm=TRUE)
maxleftloc <- which.max(left)
nrowleft <- nrow(left)
ncolleft <- ncol(left)

leftremain <- maxleftloc/nrowleft - floor(maxleftloc/nrowleft)
maxleftrowloc <- nrowleft*leftremain
maxleftcolloc <- floor(maxleftloc/nrowleft)+1
if (leftremain==0) {
	maxleftrowloc <- nrowleft
	maxleftcolloc <- floor(maxleftloc/nrowleft)
}
maxleftrowloc <- maxleftrowloc+rowliml1-1
maxleftcolloc <- maxleftcolloc+colliml1-1

collimr1 <- round(vp9colnum + 0.3 * (sac9rownum - vp9rownum))
if (collimr1 > ncol(nh9d)) collimr1 <- ncol(nh9d) # in case patient is very narrow for height
collimr2 <- vp9colnum
rowlimr1 <- round(vp9rownum + 0.15 * (sac9rownum - vp9rownum))
rowlimr2 <- round(vp9rownum + 0.55 * (sac9rownum - vp9rownum))
right <- nh9d[rowlimr1:rowlimr2,collimr2:collimr1]
maxright <- max(right,na.rm=TRUE)
maxrightloc <- which.max(right)
nrowright <- nrow(right)
ncolright <- ncol(right)

rightremain <- maxrightloc/nrowright - floor(maxrightloc/nrowright)
maxrightrowloc <- nrowright*rightremain
maxrightcolloc <- floor(maxrightloc/nrowright)+1
if (rightremain==0) {
	maxrightrowloc <- nrowright
	maxrightcolloc <- floor(maxrightloc/nrowright)
}
maxrightrowloc <- maxrightrowloc+rowlimr1-1
maxrightcolloc <- maxrightcolloc+collimr2-1
points(maxleftcolloc,sizeh-maxleftrowloc,pch=17,col="blue")
points(maxrightcolloc,sizeh-maxrightrowloc,pch=17,col="blue")

scapdiffheight <- round(-(maxrightrowloc-maxleftrowloc)*deltay)  #if +ve, right higher up back than left

############################################################################
# calculate  fifth order polynomial fitted to x and y coordinates of spinal markers
# add one extra point above and below spine markers to force spline
# gradient to be almost vertical at top and bottom of spine
	x0spine <- c(sp9[1,1],sp9[,1],sp9[numspine+1,1]) #just 6 or 8 values +2 extra in mm
	y0spine <- c((sp9[1,2]-40),sp9[,2],(sp9[numspine+1,2]+40)) #just 6 or 8 values +2 extra in mm
# normalise first
	backlength <- sac9[2]
	x0spine <- x0spine/backlength
	y0spine <- -y0spine/backlength
	xmin <- xmin/backlength
	yng <- -yng/backlength
	xng <- xng/backlength
	deltaxnorm <- deltax/backlength
# calculate fifth order polynomial coefficients
	pc <- c(1,1,1,1,1,1) #starting guesses for coefficients
	non <- nlm(mypoly,pc,gradtol=1e-8,myx=y0spine,myy=x0spine)
	ps <- non$estimate #coefficients of fifth order polynomial, ascending order 
	x1spine <- polyval(ps,yng) #x1spine is x values for the fitted 5th order poly for all yng
	x1sppix <- round((x1spine-xmin)/deltaxnorm+1) #pixel position of fitted spine line
	z1spine <- seq(1,sizeh)*0
	for (i in 1:sizeh){
		z1spine[i] <- nh9d[i,x1sppix[i]] #get z value associated with each pair
	}
	
# get distance (in mm) from spinous processes line to max points on scapulae
	maxleftoff <- round((x1sppix[maxleftrowloc]-maxleftcolloc)*deltax)
	maxrightoff <- round((maxrightcolloc-x1sppix[maxrightrowloc])*deltax) 
	
##################### TRANSVERSE PLOTS ###########################
# calculate the row numbers of 19 levels distributed
# between VP and SAC (including them)
	numlev <- 19
	lev <- matrix(0,1,numlev)
	lev[1] <- vp9rownum
	lev[numlev] <- sac9rownum
	deltalev <- (sac9rownum-vp9rownum)/(numlev-1)
	for (i in 2:(numlev-1)){
		lev[i] <- round(vp9rownum+(i-1)*deltalev)
	}
# put the sections into a matrix
	paramed <- backlength/10
	plotwidth <- round(paramed/deltax) #no of points between spine and paramedian lines
	section <- matrix(0,numlev,sizew)
	spsectpix <- seq(1,numlev)*0
	for ( i in 1:numlev){
		section[i,] <- nh9d[lev[i],]
		spsectpix[i] <- x1sppix[lev[i]]
	}
# first define width of plot wanted for each level
#const <- [0.7 0.95 1.2 1.95 2.7 2.7 2.7 2.7 2.6 2.5 2.25 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0]
	const <- c(0.7, 1.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.5, 2.25, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
# then define min and max to be plotted for each level
	minw <- spsectpix-round(plotwidth*const)
	maxw <- spsectpix+round(plotwidth*const)
# check for minw and max w being outside the data (minw < 1, maxw > sizew)
	for (i in 1:length(minw)){
		if (minw[i] < 2) minw[i] <- 2
		if (maxw[i] > (sizew-1)) maxw[i] <- sizew-1
	}
# now set everything outside these ranges to NaN
	for (i in 1:numlev){
		for (j in 1:minw[i]-1){
			section[i,j] <- NaN
		}
		for (j in (maxw[i]+1):sizew){
			section[i,j] <- NaN
		}
	}
	section2 <- section/backlength  #normalised section 
	vertsp <- seq(1,numlev)*0
	section3 <- matrix(0,numlev,sizew)
	lev2 <- seq(1,numlev)*0
	for (i in 1:numlev){
		vertsp[i] <- section2[i,spsectpix[i]] #find value of section at spinal mark
		lev2[i] <- (lev[i]-lev[1])/(lev[numlev]-lev[1]) #normalise levels
		section3[i,] <- 1*(section2[i,]-vertsp[i])-lev2[i] #create new sections
#    subtract off spinal level magnitude so all sections are zero at spinal
#    marker, subtract off normalised level for each
#    section, multiply by 2 (to give a visible result on paper and subtract the
#    zero position for that level. This spreads the sections equally over the 
#    normalised back length. (2 determined by trial and error - only for display, not
#    for quantitative results so OK.)
	}
	sections3 <- t(section3)
	x10 <- t(x9[1:19,]/backlength)
	par(mai=c(0.2,0.15,0.15,0.08))
	plot(x10,sections3,type="l",xlab="",ylab="",main="Transverse",xlim=c(-0.4,0.4),ylim=c(-1.2,0.05),axes=FALSE)
	# add paramedian lines
	paramin <- round(spsectpix-plotwidth)
	paramax <- round(spsectpix+plotwidth)
	#hold on
	for (i in 2:numlev){
		points(x10[paramin[i],1],sections3[paramin[i],i],col="red")
		points(x10[spsectpix[i],1],sections3[spsectpix[i],i],col="blue")
		points(x10[paramax[i],1],sections3[paramax[i],i],col="green")
#		points(x10[paramin[i],1]+0.019,sections3[paramin[i],i],pch=21,cex=3,col="red")
#		points(x10[spsectpix[i],1]+0.019,sections3[spsectpix[i],i],pch=21,cex=3,col="blue")
#		points(x10[paramax[i],1]+0.019,sections3[paramax[i],i],pch=21,cex=3,col="green")
	}
	# plot VP, SAC
	points(x10[vp9colnum,1],sections3[vp9colnum,1],pch=19,col="blue")
	points(x10[sac9colnum,1],sections3[sac9colnum,numlev],pch=23,col="magenta")
#	points(x10[vp9colnum,1]+0.019,sections3[vp9colnum,1],pch=21,bg="blue",col="blue",cex=3)
#	points(x10[sac9colnum,1],sections3[sac9colnum,numlev],pch=23,bg="magenta",col="magenta")
	#%calculate angle of skin across paramedian range for each level
	ang <- matrix(0,1,numlev);
	for (i in 2:numlev){
		ang[i] <- 57.2958*atan((section[i,paramax[i]]-section[i,paramin[i]])/(2*plotwidth*deltax))
		nantrig <- is.na(ang[i])
		if (nantrig){
			ang[i] <- 57.2958*atan((section[i,paramax[i]-5]-section[i,paramin[i]+5])/((2*plotwidth-10)*deltax));
		}
	}
	for (i in 2:numlev){
		angp <- round(ang[i])
		textn <- substitute(paste(angp,degree),list(angp=angp))
		text(0.35,-lev2[i]+0.01,textn,cex=1.0,pos=1,offset=0.05)
	}
	rotp <- round(rot)
	text1 <- substitute(paste("Rotation = ",rotp,degree),list(rotp=rotp))
	text3 <- sprintf("Backlength = %5.0f mm",backlength)
	text(0,-1.07,text1,cex=1.0)
	text(0,-1.13,text3,cex=1.0)
	text(0.3,0,"Skin angle",cex=1.0)
	
	text(-0.25,-1.22,"Code7701",cex=1.0)
	
	# calculate max absolute skin angle but with correct sign
	maxskinang <- round(max(abs(ang),na.rm=TRUE))
	maxskin <- round(max(ang,na.rm=TRUE))
	minskin <- round(min(ang,na.rm=TRUE))
	if (abs(minskin) > abs(maxskin)) maxskinang <- -maxskinang

#################### CORONAL PLOT WITH COBB ##################################
# calculate vertebral line from line of spinous processes and replot
	k <- c(0, 0.4, 1.8, 2.2, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.0, 0.0)
	deltab <- 1/(numlev-1)
	len <- matrix(0,1,numlev)
	p <- matrix(0,1,numlev)
	shiftx <- matrix(0,1,numlev)
	p[1] <- 0.03*backlength
	for (i in 2:numlev){
		len[i] <- deltab*(i-1)
		p[i] <- 0.03*(1+3*len[i])*backlength
	}
	shiftx <- p*sin(k*ang/57.2958)/backlength
	shiftx[numlev] <- 0 # just in case bottom value of ang is NaN
	xspinproc <- seq(1:numlev)*0
	xshift <- seq(1:numlev)*0
	yspinproc <- seq(1:numlev)*0
	for (i in 1:numlev){
		xspinproc[i] <- x1spine[lev[i]]
		xshift[i] <- xspinproc[i]+shiftx[i]
		yspinproc[i] <- yng[lev[i]]
	}
	xvert <- c(xshift[1],xshift,xshift[numlev])
	yvert <- c(yspinproc[1]+0.08,yspinproc,yspinproc[numlev]-0.08)
	non2 <- nlm(mypoly,pc,gradtol=1e-8,myx=yvert,myy=xvert)
	pv <- non2$estimate #coefficients of fifth order polynomial, ascending order 
	newspine <- polyval(pv,yng) #x values for the fitted 5th order poly for all yng
	vert <- seq(1:numlev)*0
	plot(xspinproc,yspinproc,type="l",main ="Coronal",col="blue",lwd=2,xlim=c(-0.4,0.4),ylim=c(-1.2,0.05),xlab="",ylab="",axes=FALSE)
#	title("Coronal")
	imb=round(sac9[1])
	text4 <- paste("Imbalance = ",imb," mm")
	text(0,-1.08,text4,cex=1.0)
	# plot from 0.06 to -1.06 just for the simulated spinal line
	yplus <- round((lev[19]-lev[1]+1)*0.06)
	toplev <- (lev[1]-yplus)
	toplev <- ifelse(toplev < 1,1,toplev)
	botlev <- (lev[19]+yplus)
	botlev <- ifelse((botlev>length(yng)),length(yng),botlev)
	lines(newspine[toplev:botlev],yng[toplev:botlev],lty=2,lwd=2)
#	lines(vert,yspinproc,lty=3,lwd=2)
	vert1 <- c(0,0)
	yvert <- c(0.06,-1.06)
	lines(vert1,yvert,col="blue")
	lines(xspinproc+paramed/backlength,yspinproc,col="green",lwd=2)
	lines(xspinproc-paramed/backlength,yspinproc,col="red",lwd=2)
	points(dvl8[1]/backlength,-dvl8[2]/backlength,pch=19,col="blue")
	points(dvr8[1]/backlength,-dvr8[2]/backlength,pch=19,col="blue")
	points(sac9[1]/backlength,-sac9[2]/backlength,pch=23,col="magenta")
	points(vp9[1]/backlength,-vp9[2]/backlength,pch=19,col="blue")
#	points(dvl8[1]/backlength+0.019,-dvl8[2]/backlength,pch=21,bg="blue",col="blue",cex=3)
#	points(dvr8[1]/backlength+0.019,-dvr8[2]/backlength,pch=21,bg="blue",col="blue",cex=3)
#	points(sac9[1]/backlength,-sac9[2]/backlength,pch=23,col="magenta")
#	points(vp9[1]/backlength+0.019,-vp9[2]/backlength,pch=21,bg="blue",col="blue",cex=3)
	# now calculate normals to inflection points in newspine and get Cobb
	# first find roots of second derivative equation to get inflection points
	coeffv <- c(2*pv[3],6*pv[4],12*pv[5],20*pv[6])
	rv <- polyroot(coeffv)#polyroot always gives complex results
	rv <- fliplr(sort(rv)) #gets roots into highest to lowest order, top to bottom of back
	xcalcv <- Re(polyval(pv,rv))
	# now calculate angle of perpendicular to tangent at these points
	slopev <- Re(pv[2]+2*pv[3]*rv+3*pv[4]*rv^2+4*pv[5]*rv^3+5*pv[6]*rv^4)
	# angle of normal is just negative of this
	cobbv <- 0.9*(-atan(slopev)*57.2957795)
	anglev <- -atan(slopev)
	# logic to eliminate inflection points outside range
	rv2 <- c(1,1,1)
	for (ii in 1:3){
		if (Re(rv[ii]) > 0.06) rv2[ii] <- 0
			if (Re(rv[ii]) < -1.06) rv2[ii] <- 0
#				if (abs(Im(rv[ii]))>0.0000001) rv2[ii] <- 0
	}
	# then check for complex conjugate pair and eliminate 1 if need be
	if ((Re(rv[1])-Re(rv[2])) < 0.00001) rv2[1] <- 0
	if ((Re(rv[2])-Re(rv[3])) < 0.00001) rv2[3] <- 0
	#
	nptv <- sum(rv2) #1 means no Cobb, 2 means 1 Cobb value, 3 means 2 Cobb values
	# get positions of inflection pts for possibly determining thoracic/thoracolumbar/lumbar later
	rvp <- round(Re(rv),2)  
	# do plotting only where rv2 is 1
	xgradplot <- c(-0.2,-0.1,0,0.1,0.2)
	for (ii in 1:3){
		if (rv2[ii] > 0){
			ygradv <- -xgradplot*slopev[ii]+Re(rv[ii])+slopev[ii]*xcalcv[ii]
			lines(xgradplot,ygradv,lty=3,lwd=2)
			cobbvp <- round(cobbv[ii],1)
			text1 <- substitute(paste(cobbvp,degree),list(cobbvp=cobbvp))
			text(0.3,Re(rv[ii]),text1,cex=1.0)
		}
	}
	# add lateral asymmetry text
	upcobbv <- round((cobbv[2]-cobbv[1]),6)
	textup <- "R"
	textlow <- "R"
	if (upcobbv > 0) textup <- "L"
	lowcobbv <- round((cobbv[3]-cobbv[2]),6)
	if (lowcobbv > 0) textlow <- "L"
	upcobbv2 <- round(abs(upcobbv))
	lowcobbv2 <- round(abs(lowcobbv))
	upcobbv2 <- round(abs(upcobbv))
	lowcobbv2 <- round(abs(lowcobbv))
	if (nptv == 3){
		txt5 <- substitute(paste("Lateral asymmetry = ",upcobbv2,degree,textup,
				" and ",lowcobbv2,degree,textlow),
				list(upcobbv2=upcobbv2, textup=textup, lowcobbv2=lowcobbv2,textlow=textlow))
		text(0,-1.14,txt5,cex=1.0)
	}
	if (nptv == 2){
		if (rv2[1]) txt5 <- substitute(paste("Lateral asymmetry = ",upcobbv2,degree,textup),
							list(upcobbv2=upcobbv2, textup=textup))
		if (rv2[3]) txt5 <- substitute(paste("Lateral asymmetry =  ",lowcobbv2,degree,textlow),
							list(lowcobbv2=lowcobbv2,textlow=textlow))
		text(0,-1.14,txt5,cex=1.0)
	} 
	
######### calculate bilateral asymmetry for plotting later #######
######### needed now for volumetric asymmetry ####################
	bam0 <- matrix(0,sizeh,sizew) #set up zeroed results matrix
	bam <- matrix(NaN,sizeh,sizew) #set up NaN results matrix
	# calculate difference between the two sides of the back
	# with the spinal curve as the axis of symmetry
	nw <- 3*plotwidth #normal
	if (nw > min(x1sppix)) nw <- min(x1sppix)-1
	if (nw > (sizew-max(x1sppix))) nw <- (sizew-max(x1sppix))-1
#	n2 <- 2.3*plotwidth #chris jones
	for (i in lev[1]:lev[numlev]){
		n <- x1sppix[i]
		test1 <- nh9d[i,(n+1):(n+nw)]
		test2 <- nh9d[i,(n-nw):(n-1)]
		test3 <- fliplr(test2)
		delt <- test1-test3
		bam0[i,(n+1):(n+nw)] <- delt
		delt2 <- fliplr(delt)
		bam0[i,(n-nw):(n-1)] <- delt2
	}
	# define locations of any NaNs
	tftemp <- is.na(bam0)
	# define area to be examined (exclude edges)
	const2 <- matrix(0,1,sizeh)
	# up to lev[4] is zero, below lev[19] is zero
	const2[lev[4]:lev[9]] <- 2.0 # main shoulder area
	const2[lev[13]:lev[numlev]] <- 1.5 #lower back region
	for (i in (lev[9]+1):(lev[13]-1)) { # transition shoulder area to lower back
		const2[i] <- 2.0-0.5*(i-lev[9])/(lev[13]-lev[9])
	}
	minw2 <- x1sppix-round(plotwidth*const2) #min x pixel of interest at each level
	maxw2 <- x1sppix+round(plotwidth*const2) #max x pixel of interest at each level
	# set values on the side that is higher and put opposite side to NaNs
	for (i in 1:sizeh){
		for (j in minw2[i]:(x1sppix[i]-1)) {
			if (!(tftemp[i,j])) {
				if (bam0[i,j] < 0.0) bam[i,j] <- -bam0[i,j]
			}
		}
			for (j in (x1sppix[i]+1):maxw2[i]){
				if (!(tftemp[i,j])) {
					if (bam0[i,j] > 0.0) bam[i,j] <- bam0[i,j]
				}
			}
	}
# calculate area under curves at each of the transverse levels from 2 to 19
	area <- seq(1:numlev)*0
	for (i in 2:19){
		area[i] <- sum(bam0[lev[i],minw2[lev[i]]:x1sppix[lev[i]]],na.rm=TRUE)*deltax
	}
# add volumetric data lines to coronal plot
	slice <- (lev[19]-lev[1]+1)*deltay/18
	volasym <- area * slice/(2*backlength*backlength) # 2 is fiddle factor
	volx <- cbind(xspinproc,xspinproc+volasym)
	voly <- cbind(yspinproc,yspinproc)
	for (ii in 4:19) {
		lines(volx[ii,],voly[ii,])
	}
# calculate overall volumetric parameters and add to plot
	volr <- sum(volasym[volasym>0],na.rm=TRUE)*50 #50 is fiddle factor
	voll <- abs(sum(volasym[volasym<0],na.rm=TRUE)*50)
	volr2 <- round(volr)
	voll2 <- round(voll)
		txt7 <- sprintf("Volumetric asymmetry = %4.0f %s   %4.0f %s",voll,"L",volr,"R")
		text(0,-1.19,txt7,cex=1.0)
		
#########  SAGITTAL PLOTS ####################
	off <- c(-0.24*backlength,0,0.24*backlength)
	st <- sin(-tilt*2*pi/360)
	ct <- cos(-tilt*2*pi/360)
	zoldvp <- seq(1,length(x1sppix))*0
	for (ii in 1:length(x1sppix)) {
		n2 <- x1sppix[ii]
		zoldvp[ii] <- nh9d[ii,n2]
	}
	yold <- yng*backlength
	znewc <- zoldvp*ct+yold*st
	ynewc <- -zoldvp*st+yold*ct
	plot(znewc[vp9rownum:sac9rownum],ynewc[vp9rownum:sac9rownum],type="l",lwd=2,col="blue",xlim=c(-0.24*backlength,0.36*backlength),ylim=c(-1.2*backlength,0.05*backlength),asp=1,xlab="",ylab="",axes=FALSE)
	title("Sagittal")
	# left paramedian 
	zoldparamin <- seq(1,length(x1sppix))*0
	for (ii in 1:length(x1sppix)) {
		n2 <- x1sppix[ii]-plotwidth
		zoldparamin[ii] <- nh9d[ii,n2]
	}
	znewpl <- zoldparamin*ct+yold*st
	ynewpl <- -zoldparamin*st+yold*ct
	lines(znewpl[vp9rownum:sac9rownum]+off[1],ynewpl[vp9rownum:sac9rownum],type="l",lwd=2,col="red")
	# right paramedian
	zoldparamax <- seq(1,length(x1sppix))*0
	for (ii in 1:length(x1sppix)) {
		n2 <- x1sppix[ii]+plotwidth
		zoldparamax[ii] <- nh9d[ii,n2]
	}
	znewpr <- zoldparamax*ct+yold*st
	ynewpr <- -zoldparamax*st+yold*ct
	lines(znewpr[vp9rownum:sac9rownum]+off[3],ynewpr[vp9rownum:sac9rownum],type="l",lwd=2,col="green")
	
	
# add lines between VP and SAC
	zvc <- zoldvp*0.0
	znewvc <- zvc*ct+yold*st
	ynewvc <- -zvc*st+yold*ct
	lines(znewvc[(vp9rownum-5):(sac9rownum+5)]+off[1],ynewvc[(vp9rownum-5):(sac9rownum+5)],col="red")
	lines(znewvc[(vp9rownum-5):(sac9rownum+5)]+off[2],ynewvc[(vp9rownum-5):(sac9rownum+5)],col="blue")
	lines(znewvc[(vp9rownum-5):(sac9rownum+5)]+off[3],ynewvc[(vp9rownum-5):(sac9rownum+5)],col="green")

# calculate max and min, relative to line through VP and SAC (perpendicular distance)
	#sagittalsect <- cbind(nh9d[vp9rownum:sac9rownum,paramin[1]], nh9d[vp9rownum:sac9rownum,vp9colnum], nh9d[vp9rownum:sac9rownum,paramax[1]])
	sagittalsect <- cbind(zoldparamin[vp9rownum:sac9rownum], zoldvp[vp9rownum:sac9rownum], zoldparamax[vp9rownum:sac9rownum])
	#only look for max in top two thirds
	#miss out top 10% of backlength in case of discontinuities over shoulders
	pos0 <- round(nrow(sagittalsect)/10)
	saglen <- round(2*nrow(sagittalsect)/3)
	mx <- apply(sagittalsect[pos0:saglen,],2,max,na.rm=TRUE)
	rmax <- apply(sagittalsect[pos0:saglen,],2,which.max)
	#only look below max kyphosis
	pos1 <- max(rmax)
	pos2 <- sac9rownum-vp9rownum+1
	mn <- apply(sagittalsect[pos1:pos2,],2,min,na.rm=TRUE)
	rmin <- apply(sagittalsect[pos1:pos2,],2,which.min)
	rmax <- rmax+vp9rownum+pos0-2 #positions of max values relative to beginning of array
	rmin <- rmin+vp9rownum+pos1-2 #positions of min values relative to beginning of array

# add max kyphosis and lordosis lines
	for (ii in 1:3){
		zh1 <- c(0, mx[ii])
		yh1 <- c(yold[rmax[ii]] ,yold[rmax[ii]])
		zmaxh2 <- zh1*ct+yh1*st+off[ii]
		ymaxh2 <- -zh1*st+yh1*ct
		lines(zmaxh2,ymaxh2)
		zh1 <- c(0, mn[ii])
		yh1 <- c(yold[rmin[ii]] ,yold[rmin[ii]])
		zminh2 <- zh1*ct+yh1*st+off[ii]
		yminh2 <- -zh1*st+yh1*ct
		lines(zminh2,yminh2)
		zpos <- zmaxh2[2]+14
		ypos <- ymaxh2[2]
		# stop max kyph value falling off edge of plot when it is > 60
#		if (max(mx) > 60) {
		if (max(mx)/backlength > 0.165) {
			zpos <- zmaxh2[2]-25
			ypos <- ymaxh2[2]+12
		}
	#	text(zmaxh2[2]+14,ymaxh2[2],formatC(mx[ii],digits=0,format="f"),cex=1.0)
	#	text(zmaxh2[2]-50,ymaxh2[2]+10,formatC(mx[ii],digits=0,format="f"),cex=1.0)
		text(zpos,ypos,formatC(mx[ii],digits=0,format="f"),cex=1.0)
		text(zminh2[2]-20,yminh2[2],formatC(mn[ii],digits=0,format="f"),cex=1.0)
	}
	tiltp <- round(tilt)
	text2 <- substitute(paste("Extension = ",tiltp,degree),list(tiltp=tiltp))
	if (tilt>0) text2 <- substitute(paste("Flexion = ",tiltp,degree),list(tiltp=tiltp))
	text((0.12*backlength),-(1.05*backlength),text2,cex=1.0)
	textkl <- "Max kyphosis and lordosis on plot in mm"
	text((0.12*backlength),-(1.10*backlength),textkl,cex=1.0)
	maxkyph <- round(max(mx))
	maxlord <- round(abs(min(mn)))
	
# calculate kyphosis and lordosis angles on a curve fitted to the central sagittal skin curve
# calculate fifth order polynomial coefficients for central sagittal section
	ykyph1 <- yold[vp9rownum:sac9rownum]/backlength # get rid of NaNs
	zkyph1 <- zoldvp[vp9rownum:sac9rownum]/backlength
	pc2 <- c(1,1,1,1,1,1) #starting guesses for coefficients
	non2 <- nlm(mypoly,pc2,gradtol=1e-8,myx=ykyph1,myy=zkyph1)
	ps2 <- non2$estimate #coefficients of fifth order polynomial, ascending order 
	zkyph2 <- polyval(ps2,ykyph1) #z1spine is z values for the fitted 5th order poly for all yng
# now calculate normals to inflection points in newspine and get kyphosis and lordosis angles
# first find roots of second derivative equation to get inflection points
	coeffvk <- c(2*ps2[3],6*ps2[4],12*ps2[5],20*ps2[6])
	rvk <- polyroot(coeffvk)#polyroot always gives complex results
	rvk <- fliplr(sort(rvk)) #gets roots into highest to lowest order, top to bottom of back
	xcalcvk <- Re(polyval(ps2,rvk))
# now calculate slope of tangent at the central point
	slopevk <- Re(ps2[2]+2*ps2[3]*rvk+3*ps2[4]*rvk^2+4*ps2[5]*rvk^3+5*ps2[6]*rvk^4)
# decide which root to take for calculation of middle slope
	rvkind <- c(1,1,1)
	# mark locations of inflection points between 0 and -1 with 1, 0 otherwise
	for (ii in 1:3){
		if(Re(rvk[ii]) > 0) rvkind[ii] <- 0
		if(Re(rvk[ii]) < -1) rvkind[ii] <- 0
	}
	# check for complex conjugate pairs
	for (ii in 1:3){
		if(abs(Im(rvk[ii])) > 0.000001) rvkind[ii] <- 0
	}
	slopevk2 <- slopevk*rvkind #only slopes left are for real inflection points on back
	locind <- which.max(slopevk2) #find position of positive slope (midback inflection)	
# position of point of inflection as proportion of backlength
# for sensible results this should be between -0.05 and -0.95
	inflect <- round(Re(rvk[locind]),2)
if ((inflect < -0.05)&(inflect > -0.95)) {
	kyphprint <- TRUE
} else {
	kyphprint <- FALSE
}

# angle of normal is just negative of tangent angle
	midangle <- -atan(slopevk[locind])*57.2957795 #only use central value
# get gradient of tangent at all locations
	myslope <- Re(ps2[2]+2*ps2[3]*ykyph1+3*ps2[4]*ykyph1^2+4*ps2[5]*ykyph1^3+5*ps2[6]*ykyph1^4)

# calculate average tangential angle 21 points around 0.05 down from VP and above sacrum	
	blpix <- sac9rownum-vp9rownum+1
	topkyphrow <- round(0.05*(blpix))
	if (topkyphrow < 11) topkyphrow <- 11 # make sure rows to be summed are real (row 1 or greater)
	topkyphslope <- sum(myslope[(topkyphrow-10):(topkyphrow+10)])/21
	topkyphangle <- -atan(topkyphslope)*57.2957795
	
# calculate average tangential angle over first 21 points down from VP	
	topkyphslope2 <- sum(myslope[1:21])/21
	topkyphangle2 <- -atan(topkyphslope2)*57.2957795
	
# calculate average tangential angle 21 points just above sacrum	
	botlordslope <- sum(myslope[(length(myslope)-20):(length(myslope))])/21
	botlordangle <- -atan(botlordslope)*57.2957795
	
	kyphangle <- round(midangle-topkyphangle,4)
	lordangle <- round(botlordangle-midangle,4)
	kyphangle2 <- round(midangle-topkyphangle2,4)
	
# add angles to plot
	kyphanglep <- round(abs(kyphangle))
	lordanglep <- round(lordangle)
	
	txtkl2 <- substitute(paste("Kyphosis = ",kyphanglep,degree," Lordosis = ",lordanglep,degree),
						list(kyphanglep=kyphanglep,lordanglep=lordanglep))
	if (kyphprint){
	text((0.12*backlength),-(1.15*backlength),txtkl2,cex=1.0)
#	txt6 <- sprintf("Location of kyph/lord inflection= %5.2f",inflect)
#	text((0.12*backlength),-(1.20*backlength),txt6,cex=1.0)
	}
	
#	add locations of slope measurement for kyph and lord angles
# point for upper end of kyphosis measurement
#points(znewc[vp9rownum+topkyphrow],ynewc[vp9rownum+topkyphrow],pch=19,col="black")
# point for lower end of kyphosis measurement/ upper end lordosis
inflectrownum <- vp9rownum+round(-inflect*(sac9rownum-vp9rownum+1))
#points(znewc[inflectrownum],ynewc[inflectrownum],pch=19,col="black")
# point for lower end of lordosis measurement 
#points(znewc[sac9rownum-10],ynewc[sac9rownum-10],pch=19,col="black")

if (kyphprint) {
#now add perpendiculars at these points
# central point
midx <- c(zoldvp[inflectrownum]-0.05*backlength,zoldvp[inflectrownum],zoldvp[inflectrownum]+0.05*backlength)
midy <- c(yold[inflectrownum]-0.05*backlength*sin(midangle/57.2957795),yold[inflectrownum],yold[inflectrownum]+0.05*backlength*sin(midangle/57.2957795))
midxnew <- midx*ct+midy*st
midynew <- -midx*st+midy*ct

lines(midxnew,midynew)
# upper point
midx2 <- c(zoldvp[vp9rownum+topkyphrow]-0.05*backlength,zoldvp[vp9rownum+topkyphrow],zoldvp[vp9rownum+topkyphrow]+0.05*backlength)
midy2 <- c(yold[vp9rownum+topkyphrow]-0.05*backlength*sin(topkyphangle/57.2957795),yold[vp9rownum+topkyphrow],yold[vp9rownum+topkyphrow]+0.05*backlength*sin(topkyphangle/57.2957795))
midxnew2 <- midx2*ct+midy2*st
midynew2 <- -midx2*st+midy2*ct
lines(midxnew2,midynew2)
# lower point
 midx3 <- c(zoldvp[sac9rownum-10]-0.05*backlength,zoldvp[sac9rownum-10],zoldvp[sac9rownum-10]+0.05*backlength)
 midy3 <- c(yold[sac9rownum-10]-0.05*backlength*sin(botlordangle/57.2957795),yold[sac9rownum-10],yold[sac9rownum-10]+0.05*backlength*sin(botlordangle/57.2957795))
 midxnew3 <- midx3*ct+midy3*st
 midynew3 <- -midx3*st+midy3*ct
 lines(midxnew3,midynew3)
}

if (!kyphprint) {
 txt6 <- sprintf("Kyphosis/lordosis angles - no results available")
 text((0.12*backlength),-(1.20*backlength),txt6,cex=1.0)
}

 dev.off()
 time32 <- Sys.time()
################## BILATERAL ASYMMETRY MAPS ######################
### bam0 (difference between the sides) was calculated in coronal plot section ####
#### define area to be examined (exclude edges)
	swidth <- 2.7
	if (round(nw/plotwidth,1) < swidth) swidth <- round(nw/plotwidth,1)
	vpwidth <- 0.7
	lbackwidth <- 2.0
	const2 <- matrix(0,1,sizeh)
	const2[lev[1]] <- 0.3 # VP level
	const2[lev[3]:lev[9]] <- swidth # main shoulder area
	const2[lev[12]:lev[numlev]] <- 2.0 #lower back region
	for (i in (lev[1]+1):(lev[3]-1)) { # transition VP to main shoulder area
		const2[i] <- vpwidth+(swidth-vpwidth)*(i-lev[1])/(lev[3]-lev[1])
	}
	for (i in (lev[9]+1):(lev[12]-1)) { # transition shoulder area to lower back
		const2[i] <- swidth-(swidth-lbackwidth)*(i-lev[9])/(lev[12]-lev[9])
	}
	minw2 <- x1sppix-round(plotwidth*const2) #min x pixel of interest at each level
	maxw2 <- x1sppix+round(plotwidth*const2) #max x pixel of interest at each level
	# set values on the side that is higher and put opposite side to NaNs
	for (i in 1:sizeh){
		for (j in minw2[i]:(x1sppix[i]-1)) {
			if (!(tftemp[i,j])) {
				if (bam0[i,j] < 0.0) bam[i,j] <- -bam0[i,j]
			}
		}
			for (j in (x1sppix[i]+1):maxw2[i]){
				if (!(tftemp[i,j])) {
					if (bam0[i,j] > 0.0) bam[i,j] <- bam0[i,j]
				}
			}
	}
	time33 <- Sys.time()

	# display bilateral asymmetry
	a<- bam
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
	zlim <- range(aplot,finite=FALSE,na.rm=TRUE)
	minw2[1:(lev[1]-2)] <- NaN
	minw2[(lev[numlev]+2):sizeh] <- NaN
	maxw2[1:(lev[1]-2)] <- NaN
	maxw2[(lev[numlev]+2):sizeh] <- NaN
	leftline <- fliplr(minw2)
	rightline <- fliplr(maxw2)
	#calculate suitable plot limits
	yup <- max(yplot)-lev[1]+fact+1
	ylow <- max(yplot)-lev[19]-fact-1
	yspread <- (yup-ylow)/2 #gives half the range to be plotted for y, x must be same
	# but spread equally on both sides of the sacrum
	xmid <-max(xplot)/2
	xlow <- xmid-yspread
	xup <- xmid+yspread
	#check that x is not outside the xplot range
	if (xlow < xplot[1]) xlow <- xplot[1]
	if (xup > xplot[ncols]) xup <- xplot[ncols]
	#first plot all differences
#	oldpar <- par(no.readonly=TRUE)
	pdffile <- "bamall.pdf"
	pdf(pdffile,width=7.0,height=7.5,bg="white",onefile=TRUE,family="Times")
	par(oma=c(0,0,0,1.5),mar=c(1,1,3,1))
	filled.contour(xplot,yplot,fliplr(t(aplot)),frame.plot=FALSE,ylim=c(ylow,yup),xlim=c(xlow,xup),asp=1,xlab="",ylab="",levels=seq(0,zlim[2],zlim[2]/64),col=mycol,plot.title=title(main="Bilateral Asymmetry Map\nall differences"),key.title=title(main="mm"),plot.axes={lines(leftline,1:sizeh,lwd=2);lines(rightline,1:sizeh,lwd=2)})
	#mtext(aText,line=0.5,outer=TRUE,font=2)
	dev.off()
	#add ;mtext(text,side=1,line=1) to within the plot.axes statement if you want text added at the bottom
	#now plot only differences>10mm, scaled to max of 30
	aplot2 <- ifelse(aplot<10,NaN,aplot) #less than 10 mm not to be plotted
	aplot2 <- ifelse(aplot>30,30,aplot2) #greater than 30 set to 30
	pdffile <- "bamhigh.pdf"
	pdf(pdffile,width=7.0,height=7.5,bg="white",onefile=TRUE,family="Times")
	par(oma=c(0,0,0,1.5),mar=c(1,1,3,1))
	filled.contour(xplot,yplot,fliplr(t(aplot2)),asp=1,frame.plot=FALSE,ylim=c(ylow,yup),xlim=c(xlow,xup),xlab="",ylab="",levels=seq(0,30,30/64),col=mycol,plot.title=title(main="Bilateral Asymmetry Map\ndifferences >10mm"),key.title=title(main="mm"),plot.axes={lines(leftline,1:sizeh,lwd=2);lines(rightline,1:sizeh,lwd=2)})
	#mtext(aText,line=0.5,outer=TRUE,font=2)
#	par(oldpar)
	dev.off()
	time34 <- Sys.time()

	#add a text file output for reading in SuperCard
	zz <- file("results.txt","wt")
	outputdata <- c(round(backlength),round(rot),round(tilt),imb,nptv,rv2,rvp,round(upcobbv),textup,
					round(lowcobbv),textlow,maxskin,minskin,maxkyph,maxlord,kyphanglep,lordanglep,
					inflect,scapdiffheight,maxleftoff,maxrightoff,voll2,volr2)
	write(outputdata,zz)
	close(zz)
	time35 <- Sys.time()
	deltat <- (time35-times)
	cat(sprintf('Processing time is %g\n',deltat))
	cat(sprintf('Run completed successfully\n'))
#END