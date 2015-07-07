
library(fields)
set.seed(3245)

n <- 800
x <- sort(runif(n, 0,1))
xgrid <- seq(0, 1,,1000)
f1 <- sin((2*pi*x))	
f2 <- sin((2*pi*x)^2)	
f3 <- sin((2*pi*x)^3)	
f4 <- ifelse(x < .5, sin((2*pi*(x+0.0026))^3), -sin((2*pi*(1-x-0.0026))^3))


fg1 <- sin((2*pi*xgrid))
fg2 <- sin((2*pi*xgrid)^2)
fg3 <- sin((2*pi*xgrid)^3)
fg4 <- ifelse(xgrid < .5, sin((2*pi*(xgrid+0.0026))^3), -sin((2*pi*(1-xgrid-0.0026))^3))

snr <- 5
sdf1 <- sd(f1)
sdf2 <- sd(f2)
sdf3 <- sd(f3)
sdf4 <- sd(f4)
y1 <- f1 + rnorm(length(f1), 0, sdf1/snr)
y2 <- f2 + rnorm(length(f2), 0, sdf2/snr)
y3 <- f3 + rnorm(length(f3), 0, sdf3/snr)
y4 <- f4 + rnorm(length(f4), 0, sdf4/snr)

bdx <- c(-rev(x[x < 0.2]), x, rev(2-x[x > 0.8]))
f1 <- c(-rev(f1[x < 0.2]), f1, -rev(f1[x > 0.8]))
f2 <- c(rev(f2[x < 0.2]), f2, rev(f2[x > 0.8]))
f3 <- c(rev(f3[x < 0.2]), f3, rev(f3[x > 0.8]))
f4 <- c(rev(f4[x < 0.2]), f4, rev(f4[x > 0.8]))
y1 <- c(-rev(y1[x < 0.2]), y1, -rev(y1[x > 0.8]))
y2 <- c(rev(y2[x < 0.2]), y2, rev(y2[x > 0.8]))
y3 <- c(rev(y3[x < 0.2]), y3, rev(y3[x > 0.8]))
y4 <- c(rev(y4[x < 0.2]), y4, rev(y4[x > 0.8]))
x <- bdx


test1 <- list(x=x, y=y1, f=f1)
test2 <- list(x=x, y=y2, f=f2)
test3 <- list(x=x, y=y3, f=f3)
test4 <- list(x=x, y=y4, f=f4)

pilot1 <- pilotQV(test1$x, test1$y)
pilot2 <- pilotQV(test2$x, test2$y)
pilot3 <- pilotQV(test3$x, test3$y)
pilot4 <- pilotQV(test4$x, test4$y)


lambda.grid1 <- pilot1$lambda.grid
lambda.grid2 <- pilot2$lambda.grid
lambda.grid3 <- pilot3$lambda.grid
lambda.grid4 <- pilot4$lambda.grid

nh <- length(lambda.grid1)

newx <- seq(0,1,,300)
rval <- 0.05
ptm <- proc.time()
out1 <- improve_local(pilot1, newx, rval, sparameter=0.1^7, boundary="none")
proc.time() - ptm

out2 <- improve_local(pilot2, newx, rval, sparameter=0.1^7, boundary="none")
out3 <- improve_local(pilot3, newx, rval, sparameter=0.1^7, boundary="none")
out4 <- improve_local(pilot4, newx, rval, sparameter=0.1^7, boundary="none")


newf1 <- out1$predy[, ncol(out1$predy)]
newsigma1 <- out1$predsigma[length(out1$predsigma)]
newf2 <- out2$predy[, ncol(out2$predy)]
newsigma2 <- out2$predsigma[length(out2$predsigma)]
newf3 <- out3$predy[, ncol(out3$predy)]
newsigma3 <- out3$predsigma[length(out3$predsigma)]
newf4 <- out4$predy[, ncol(out4$predy)]
newsigma4 <- out4$predsigma[length(out4$predsigma)]

elrisk_final1 <- localrisk(pilot1, rmethod="true", pilotf=newf1 , pilotsigma=newsigma1)
elrisk_final2 <- localrisk(pilot2, rmethod="true", pilotf=newf2 , pilotsigma=newsigma2)
elrisk_final3 <- localrisk(pilot3, rmethod="true", pilotf=newf3 , pilotsigma=newsigma3)
elrisk_final4 <- localrisk(pilot4, rmethod="true", pilotf=newf4 , pilotsigma=newsigma4)

elambda_final1 <- movinglambda(elrisk_final1, newx, rval, boundary="none")
elambda_final2 <- movinglambda(elrisk_final2, newx, rval, boundary="none")
elambda_final3 <- movinglambda(elrisk_final3, newx, rval, boundary="none")
elambda_final4 <- movinglambda(elrisk_final4, newx, rval, boundary="none")


save(test1, test2, test3, test4, pilot1, pilot2, pilot3, pilot4, 
out1, out2, out3, out4, lambda.grid1, lambda.grid2, lambda.grid3, lambda.grid4,
elrisk_final1, elrisk_final2, elrisk_final3, elrisk_final4,
elambda_final1, elambda_final2, elambda_final3, elambda_final4, file="repara.RData")

##################
########################################################################

elambda <- elambda_final2$lambda
tmpx <- newx[!is.na(elambda)]
tmpy <- log(elambda[!is.na(elambda)])
fite <- Tps(tmpx, tmpy, lambda=0.1^4)
plot(tmpx, tmpy)
lines(tmpx, fite$fitted.values, col=2, lwd=2)
########################################################################
##################

set.seed(3245)
seed.vals <- sort(unique(round(runif(100, 1000, 1000000))))[1:50]


for(i in 1:50){
    cat(i, "\n")
    set.seed(seed.vals[i])
    n <- 800
    x <- sort(runif(n, 0,1))
    xgrid <- seq(0, 1,,1000)
    f1 <- sin((2*pi*x))	
    f2 <- sin((2*pi*x)^2)	
    f3 <- sin((2*pi*x)^3)	
    f4 <- ifelse(x < .5, sin((2*pi*(x+0.0026))^3), -sin((2*pi*(1-x-0.0026))^3))
    
    
    fg1 <- sin((2*pi*xgrid))
    fg2 <- sin((2*pi*xgrid)^2)
    fg3 <- sin((2*pi*xgrid)^3)
    fg4 <- ifelse(xgrid < .5, sin((2*pi*(xgrid+0.0026))^3), -sin((2*pi*(1-xgrid-0.0026))^3))
    
    snr <- 5
    sdf1 <- sd(f1)
    sdf2 <- sd(f2)
    sdf3 <- sd(f3)
    sdf4 <- sd(f4)
    y1 <- f1 + rnorm(length(f1), 0, sdf1/snr)
    y2 <- f2 + rnorm(length(f2), 0, sdf2/snr)
    y3 <- f3 + rnorm(length(f3), 0, sdf3/snr)
    y4 <- f4 + rnorm(length(f4), 0, sdf4/snr)
    
    
    test1 <- list(x=x, y=y1, f=f1)
    test2 <- list(x=x, y=y2, f=f2)
    test3 <- list(x=x, y=y3, f=f3)
    test4 <- list(x=x, y=y4, f=f4)
    
    pilot1 <- pilotQV(test1$x, test1$y)
    pilot2 <- pilotQV(test2$x, test2$y)
    pilot3 <- pilotQV(test3$x, test3$y)
    pilot4 <- pilotQV(test4$x, test4$y)
    
    
    lambda.grid1 <- pilot1$lambda.grid
    lambda.grid2 <- pilot2$lambda.grid
    lambda.grid3 <- pilot3$lambda.grid
    lambda.grid4 <- pilot4$lambda.grid
    
    nh <- length(lambda.grid1)
    
    newx <- seq(0,1,,300)
    rval <- 0.05
    out1 <- improve_local(pilot1, newx, rval, sparameter=0.1^7)
    out2 <- improve_local(pilot2, newx, rval, sparameter=0.1^7)
    out3 <- improve_local(pilot3, newx, rval, sparameter=0.1^7)
    out4 <- improve_local(pilot4, newx, rval, sparameter=0.1^7)
    newf1 <- out1$predy[, ncol(out1$predy)]
    newsigma1 <- out1$predsigma[length(out1$predsigma)]
    newf2 <- out2$predy[, ncol(out2$predy)]
    newsigma2 <- out2$predsigma[length(out2$predsigma)]
    newf3 <- out3$predy[, ncol(out3$predy)]
    newsigma3 <- out3$predsigma[length(out3$predsigma)]
    newf4 <- out4$predy[, ncol(out4$predy)]
    newsigma4 <- out4$predsigma[length(out4$predsigma)]
    
    elrisk_final1 <- localrisk(pilot1, rmethod="true", pilotf=newf1 , pilotsigma=newsigma1)
    elrisk_final2 <- localrisk(pilot2, rmethod="true", pilotf=newf2 , pilotsigma=newsigma2)
    elrisk_final3 <- localrisk(pilot3, rmethod="true", pilotf=newf3 , pilotsigma=newsigma3)
    elrisk_final4 <- localrisk(pilot4, rmethod="true", pilotf=newf4 , pilotsigma=newsigma4)
    
    elambda_final1 <- movinglambda(elrisk_final1, newx, rval)
    elambda_final2 <- movinglambda(elrisk_final2, newx, rval)
    elambda_final3 <- movinglambda(elrisk_final3, newx, rval)
    elambda_final4 <- movinglambda(elrisk_final4, newx, rval)
    save(test1, test2, test3, test4, pilot1, pilot2, pilot3, pilot4, 
    out1, out2, out3, out4, lambda.grid1, lambda.grid2, lambda.grid3, lambda.grid4,
    elrisk_final1, elrisk_final2, elrisk_final3, elrisk_final4,
    elambda_final1, elambda_final2, elambda_final3, elambda_final4, file=paste("repara_",i, ".RData", sep=""))
}


######################################################
##########################################################################################
##########################################################################################
######################################################
# Simulation 20150527
#path <-"/Volumes/multiscale/"
path <-"~/DataDisk/MacHD2/Temp_Research/Changing point/Papers/R_code/repara/"
library(fields)

predlambda1 <- matrix(NA, 300, 50)
predlambda2 <- matrix(NA, 300, 50)
predlambda3 <- matrix(NA, 300, 50)
predlambda4 <- matrix(NA, 300, 50)

llambda1 <- list()
llambda2 <- list()
llambda3 <- list()
llambda4 <- list()
newx <- seq(0,1,,300)


for(i in 1:50){
  load(paste(path, "repara_test1_", i, ".RData", sep=""))
  fite1 <- sreg(newx, -log(elambda_final1$lambda), lambda=0.1^4)
  predlambda1[ ,i] <- predict(fite1, newx)
  llambda1[[i]] <- cbind(out1$inilrisk$xM, -log(out1$inilrisk$lambda.grid[apply(out1$inilrisk$risk, 1, which.min)]))
  
  load(paste(path, "repara_test2_", i, ".RData", sep=""))
  fite2 <- sreg(newx, -log(elambda_final2$lambda), lambda=0.1^4)
  predlambda2[ ,i] <- predict(fite2, newx)
  llambda2[[i]] <- cbind(out2$inilrisk$xM, -log(out2$inilrisk$lambda.grid[apply(out2$inilrisk$risk, 1, which.min)]))
  
  load(paste(path, "repara_test3_", i, ".RData", sep=""))
  fite3 <- sreg(newx, -log(elambda_final3$lambda), lambda=0.1^4)
  predlambda3[ ,i] <- predict(fite3, newx)
  llambda3[[i]] <- cbind(out3$inilrisk$xM, -log(out3$inilrisk$lambda.grid[apply(out3$inilrisk$risk, 1, which.min)]))
  
  load(paste(path, "repara_test4_", i, ".RData", sep=""))
  fite4 <- sreg(newx, -log(elambda_final4$lambda), lambda=0.1^4)
  predlambda4[ ,i] <- predict(fite4, newx)
  llambda4[[i]] <- cbind(out4$inilrisk$xM, -log(out4$inilrisk$lambda.grid[apply(out4$inilrisk$risk, 1, which.min)]))
}
matplot(newx, cbind(predlambda1, predlambda2, predlambda3, predlambda4), col=rep(1:4, each=50), type="l")

save(newx, predlambda1, predlambda2, predlambda3, predlambda4, 
     llambda1, llambda2, llambda3, llambda4,
     file=paste(path, "repara_all_local_n_moving.RData", sep=""))

lines(out1$inilrisk$xM, -log(out1$inilrisk$lambda.grid[apply(out1$inilrisk$risk, 1, which.min)]))


dev1 <- rep(NA, 50)
dev2 <- rep(NA, 50)
dev3 <- rep(NA, 50)
dev4 <- matrix(NA, 50, 2)
newx <- seq(0,1,,300)


for(i in 1:50){
    load(paste(path, "repara_1", i, ".RData", sep=""))
    fite1 <- sreg(newx, -log(elambda_final1$lambda), lambda=0.1^5)
    predlambda1 <- predict(fite1, newx)
    #dev1 <- predict(fite1, newx, derivative=1)
    ind <- which(newx> .18 & newx < .82)
    dev1[i] <- lm(predlambda1[ind] ~newx[ind])$coeff[2]
    
    load(paste(path, "repara_2", i, ".RData", sep=""))
    fite2 <- sreg(newx, -log(elambda_final2$lambda), lambda=0.1^5)
    predlambda2 <- predict(fite2, newx)
    #dev2 <- predict(fite2, newx, derivative=1)
    dev2[i] <- lm(predlambda2[ind] ~newx[ind])$coeff[2]
    
    load(paste(path, "repara_3", i, ".RData", sep=""))
    fite3 <- sreg(newx, -log(elambda_final3$lambda), lambda=0.1^5)
    predlambda3 <- predict(fite3, newx)
    #dev3 <- predict(fite3, newx, derivative=1)
    dev3[i] <- lm(predlambda3[ind] ~newx[ind])$coeff[2]
    
    load(paste(path, "repara_4", i, ".RData", sep=""))
    fite4 <- sreg(newx, -log(elambda_final4$lambda), lambda=0.1^5)
    predlambda4 <- predict(fite4, newx)
    #dev4 <- predict(fite3, newx, derivative=1)
    ind <- which(newx> .02 & newx < 0.4)
    #    dev4 <- c(0,0)
    dev4[i,1] <- lm(predlambda4[ind] ~newx[ind])$coeff[2]
    ind <- which(newx> 0.6 & newx < .98)
    dev4[i,2] <- lm(predlambda4[ind] ~newx[ind])$coeff[2]
}





load("/Volumes/MacHD2/Temp_Research/Changing point/Papers/R_code/repara.RData")


set.seed(3245)
n <- 800
x <- sort(runif(n, 0,1))
xgrid <- seq(0, 1,,1000)
f1 <- sin((2*pi*x))	
f2 <- sin((2*pi*x)^2)	
f3 <- sin((2*pi*x)^3)	
f4 <- ifelse(x < .5, sin((2*pi*(x+0.0026))^3), -sin((2*pi*(1-x-0.0026))^3))

fg1 <- sin((2*pi*xgrid))
fg2 <- sin((2*pi*xgrid)^2)
fg3 <- sin((2*pi*xgrid)^3)
fg4 <- ifelse(xgrid < .5, sin((2*pi*(xgrid+0.0026))^3), -sin((2*pi*(1-xgrid-0.0026))^3))

snr <- 5
sdf1 <- sd(f1)
sdf2 <- sd(f2)
sdf3 <- sd(f3)
sdf4 <- sd(f4)
y1 <- f1 + rnorm(length(f1), 0, sdf1/snr)
y2 <- f2 + rnorm(length(f2), 0, sdf2/snr)
y3 <- f3 + rnorm(length(f3), 0, sdf3/snr)
y4 <- f4 + rnorm(length(f4), 0, sdf4/snr)

test1 <- list(x=x, y=y1, f=f1)
test2 <- list(x=x, y=y2, f=f2)
test3 <- list(x=x, y=y3, f=f3)
test4 <- list(x=x, y=y4, f=f4)

newx <- seq(0,1,,300)
fite1 <- sreg(newx, -log(elambda_final1$lambda), lambda=0.1^5)
predlambda1 <- predict(fite1, newx)
dpredlambda1 <- predict(fite1, newx, derivative = 1)
fite2 <- sreg(newx, -log(elambda_final2$lambda), lambda=0.1^5)
predlambda2 <- predict(fite2, newx)
dpredlambda2 <- predict(fite2, newx, derivative = 1)
fite3 <- sreg(newx, -log(elambda_final3$lambda), lambda=0.1^5)
predlambda3 <- predict(fite3, newx)
dpredlambda3 <- predict(fite3, newx, derivative = 1)
fite4 <- sreg(newx, -log(elambda_final4$lambda), lambda=0.1^5)
predlambda4 <- predict(fite4, newx)
dpredlambda4 <- predict(fite4, newx, derivative = 1)





save(test1, test2, test3, test4, xgrid, fg1, fg2, fg3, fg4,
elambda_final1, elambda_final2, elambda_final3, elambda_final4,newx, 
predlambda1, predlambda2, predlambda3, predlambda4, dev1, dev2, dev3, dev4,
file="repara2.RData")


cairo_pdf(file="repara.pdf", width=13.26, height=9.06, onefile = TRUE);par(font.main=2)
#########################################################
par(mfrow=c(3,2), cex.main=2.5, cex.axis=1.5, cex.lab=1.5, mar=c(3,5.5,4,2))
plot(test1$x, test1$y, main="(a)", xlab="", ylab="y")
lines(xgrid, fg1, col=2 ,lwd=2)
plot(test2$x, test2$y, main="(b)", xlab="", ylab="y")
lines(xgrid, fg2, col=2 ,lwd=2)
plot(test3$x, test3$y, main="(c)", xlab="", ylab="y")
lines(xgrid, fg3, col=2 ,lwd=2)
plot(test4$x, test4$y, main="(d)", xlab="", ylab="y")
lines(xgrid, fg4, col=2 ,lwd=2)


matplot(newx, cbind(-log(elambda_final1$lambda), -log(elambda_final2$lambda), -log(elambda_final3$lambda) , -log(elambda_final4$lambda)), 
xlab="", ylab=expression(-log(hat(lambda)(x))) , type="n", main="(e)")
lines(newx, predlambda1, col=1, lwd=3)
lines(newx, predlambda2, col=2, lwd=3, lty=2)
lines(newx, predlambda3, col=3, lwd=3, lty=3)
lines(newx, predlambda4, col=4, lwd=3, lty=4)
boxplot(cbind(dev1, dev2, dev3), names=c("(a)", "(b)", "(c)"), main="(f)", ylab=expression("Slope of"~~-log(hat(lambda)(x))))

#plot(c(dev1, dev2, dev3, dev4), rep(1,5), type="n", axes=FALSE, main="(d)", xlab="", ylab="", xlim=c(-2,18))
#abline(v=c(dev1, dev2, dev3, dev4), col=c(1:4,4), lty=c(1:4,4), lwd=3)
#axis(1)
#box()
#matplot(newx, cbind(dev1, dev2, dev3), xlab="", ylab="First derivative" , type="n", main="(e)")
#lines(newx, dev1, col=1, lwd=2)
#lines(newx, dev2, col=2, lwd=2, lty=2)
#lines(newx, dev3, col=3, lwd=2, lty=3)

dev.off()
#############################################

cairo_pdf(file="repara2.pdf", width=13.26, height=9.06, onefile = TRUE);par(font.main=2)
#########################################################
par(mfrow=c(3,2), cex.main=2.5, cex.axis=1.5, cex.lab=1.5, mar=c(3,5.5,4,2))
plot(test1$x, test1$y, main="(a)", xlab="", ylab="y")
lines(xgrid, fg1, col=2 ,lwd=2)
plot(test2$x, test2$y, main="(b)", xlab="", ylab="y")
lines(xgrid, fg2, col=2 ,lwd=2)
plot(test3$x, test3$y, main="(c)", xlab="", ylab="y")
lines(xgrid, fg3, col=2 ,lwd=2)
plot(test4$x, test4$y, main="(d)", xlab="", ylab="y")
lines(xgrid, fg4, col=2 ,lwd=2)


matplot(newx, cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda)), 
xlab="", ylab=expression(log(hat(lambda)(x))) , type="n", main="(e)")
lines(newx, -predlambda1, col=1, lwd=3)
lines(newx, -predlambda2, col=2, lwd=3, lty=2)
lines(newx, -predlambda3, col=3, lwd=3, lty=3)
lines(newx, -predlambda4, col=4, lwd=3, lty=4)
boxplot(-cbind(dev1, dev2, dev3), names=c("(a)", "(b)", "(c)"), main="(f)", ylab=expression("Slope of"~~log(hat(lambda)(x))))

#plot(c(dev1, dev2, dev3, dev4), rep(1,5), type="n", axes=FALSE, main="(d)", xlab="", ylab="", xlim=c(-2,18))
#abline(v=c(dev1, dev2, dev3, dev4), col=c(1:4,4), lty=c(1:4,4), lwd=3)
#axis(1)
#box()
#matplot(newx, cbind(dev1, dev2, dev3), xlab="", ylab="First derivative" , type="n", main="(e)")
#lines(newx, dev1, col=1, lwd=2)
#lines(newx, dev2, col=2, lwd=2, lty=2)
#lines(newx, dev3, col=3, lwd=2, lty=3)

dev.off()




cairo_pdf(file="repara2.pdf", width=13.26, height=9.06, onefile = TRUE);par(font.main=2)
#########################################################
par(mfrow=c(3,2), cex.main=2.5, cex.axis=1.5, cex.lab=1.5, mar=c(3,5.5,4,2))
plot(test1$x, test1$y, main="(a)", xlab="", ylab="y")
lines(xgrid, fg1, col=2 ,lwd=2)
plot(test2$x, test2$y, main="(b)", xlab="", ylab="y")
lines(xgrid, fg2, col=2 ,lwd=2)
plot(test3$x, test3$y, main="(c)", xlab="", ylab="y")
lines(xgrid, fg3, col=2 ,lwd=2)
plot(test4$x, test4$y, main="(d)", xlab="", ylab="y")
lines(xgrid, fg4, col=2 ,lwd=2)


matplot(newx, cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda)), 
xlab="", ylab=expression(log(hat(lambda)(x))) , type="n", main="(e)")
lines(newx, -predlambda1, col=1, lwd=3)
lines(newx, -predlambda2, col=2, lwd=3, lty=2)
lines(newx, -predlambda3, col=3, lwd=3, lty=3)
lines(newx, -predlambda4, col=4, lwd=3, lty=4)
boxplot(-cbind(dev1, dev2, dev3), names=c("(a)", "(b)", "(c)"), main="(f)", ylab=expression("Slope of"~~log(hat(lambda)(x))))

#plot(c(dev1, dev2, dev3, dev4), rep(1,5), type="n", axes=FALSE, main="(d)", xlab="", ylab="", xlim=c(-2,18))
#abline(v=c(dev1, dev2, dev3, dev4), col=c(1:4,4), lty=c(1:4,4), lwd=3)
#axis(1)
#box()
#matplot(newx, cbind(dev1, dev2, dev3), xlab="", ylab="First derivative" , type="n", main="(e)")
#lines(newx, dev1, col=1, lwd=2)
#lines(newx, dev2, col=2, lwd=2, lty=2)
#lines(newx, dev3, col=3, lwd=2, lty=3)

dev.off()



cairo_pdf(file="repara3.pdf", width=13.26, height=9.06, onefile = TRUE);par(font.main=2)
#########################################################
par(mfrow=c(3,2), cex.main=2.5, cex.axis=1.5, cex.lab=1.5, mar=c(3,5.5,4,2))
plot(test1$x, test1$y, main="(a)", xlab="", ylab="y")
lines(xgrid, fg1, col=2 ,lwd=2)
plot(test2$x, test2$y, main="(b)", xlab="", ylab="y")
lines(xgrid, fg2, col=2 ,lwd=2)
plot(test3$x, test3$y, main="(c)", xlab="", ylab="y")
lines(xgrid, fg3, col=2 ,lwd=2)
plot(test4$x, test4$y, main="(d)", xlab="", ylab="y")
lines(xgrid, fg4, col=2 ,lwd=2)

aa1 <- cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda))
bb1 <- apply(aa1, 2, diff)/(newx[2] - newx[1])
dnewx <- newx[-length(newx)] + diff(newx)/2

ylim <- range(c(range(dpredlambda1), range(dpredlambda2), range(dpredlambda3), range(dpredlambda4)))
matplot(dnewx, bb1, xlab="", ylab=expression(log(hat(lambda)(x))) , type="n", main="(e)", ylim=ylim)
lines(newx, dpredlambda1, col=1, lwd=3)
lines(newx, dpredlambda2, col=2, lwd=3, lty=2)
lines(newx, dpredlambda3, col=3, lwd=3, lty=3)
lines(newx, dpredlambda4, col=4, lwd=3, lty=4)


matplot(newx, cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda)), 
xlab="", ylab=expression(log(hat(lambda)(x))) , type="n", main="(e)")
lines(newx, dpredlambda1, col=1, lwd=3)
lines(newx, dpredlambda2, col=2, lwd=3, lty=2)
lines(newx, dpredlambda3, col=3, lwd=3, lty=3)
lines(newx, dpredlambda4, col=4, lwd=3, lty=4)



#boxplot(-cbind(dev1, dev2, dev3), names=c("(a)", "(b)", "(c)"), main="(f)", ylab=expression("Slope of"~~log(hat(lambda)(x))))

#plot(c(dev1, dev2, dev3, dev4), rep(1,5), type="n", axes=FALSE, main="(d)", xlab="", ylab="", xlim=c(-2,18))
#abline(v=c(dev1, dev2, dev3, dev4), col=c(1:4,4), lty=c(1:4,4), lwd=3)
#axis(1)
#box()
#matplot(newx, cbind(dev1, dev2, dev3), xlab="", ylab="First derivative" , type="n", main="(e)")
#lines(newx, dev1, col=1, lwd=2)
#lines(newx, dev2, col=2, lwd=2, lty=2)
#lines(newx, dev3, col=3, lwd=2, lty=3)

dev.off()


###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

set.seed(3245)

n <- 800
x <- sort(runif(n, 0,1))
xgrid <- seq(0, 1,,1000)
f1 <- sin((2*pi*x))	
f2 <- sin((2*pi*x)^2)	
f3 <- sin((2*pi*x)^3)	
f4 <- ifelse(x < .5, sin((2*pi*(x+0.0026))^3), -sin((2*pi*(1-x-0.0026))^3))


fg1 <- sin((2*pi*xgrid))
fg2 <- sin((2*pi*xgrid)^2)
fg3 <- sin((2*pi*xgrid)^3)
fg4 <- ifelse(xgrid < .5, sin((2*pi*(xgrid+0.0026))^3), -sin((2*pi*(1-xgrid-0.0026))^3))

snr <- 5
sdf1 <- sd(f1)
sdf2 <- sd(f2)
sdf3 <- sd(f3)
sdf4 <- sd(f4)
y1 <- f1 + rnorm(length(f1), 0, sdf1/snr)
y2 <- f2 + rnorm(length(f2), 0, sdf2/snr)
y3 <- f3 + rnorm(length(f3), 0, sdf3/snr)
y4 <- f4 + rnorm(length(f4), 0, sdf4/snr)

bdx <- c(-rev(x[x < 0.2]), x, rev(2-x[x > 0.8]))
f1 <- c(-rev(f1[x < 0.2]), f1, -rev(f1[x > 0.8]))
f2 <- c(rev(f2[x < 0.2]), f2, rev(f2[x > 0.8]))
f3 <- c(rev(f3[x < 0.2]), f3, rev(f3[x > 0.8]))
f4 <- c(rev(f4[x < 0.2]), f4, rev(f4[x > 0.8]))
y1 <- c(-rev(y1[x < 0.2]), y1, -rev(y1[x > 0.8]))
y2 <- c(rev(y2[x < 0.2]), y2, rev(y2[x > 0.8]))
y3 <- c(rev(y3[x < 0.2]), y3, rev(y3[x > 0.8]))
y4 <- c(rev(y4[x < 0.2]), y4, rev(y4[x > 0.8]))
x <- bdx


test1 <- list(x=x, y=y1, f=f1)
test2 <- list(x=x, y=y2, f=f2)
test3 <- list(x=x, y=y3, f=f3)
test4 <- list(x=x, y=y4, f=f4)

pilot1 <- pilotQV(test1$x, test1$y)
pilot2 <- pilotQV(test2$x, test2$y)
pilot3 <- pilotQV(test3$x, test3$y)
pilot4 <- pilotQV(test4$x, test4$y)


lambda.grid1 <- pilot1$lambda.grid
lambda.grid2 <- pilot2$lambda.grid
lambda.grid3 <- pilot3$lambda.grid
lambda.grid4 <- pilot4$lambda.grid

nh <- length(lambda.grid1)

newx <- seq(0,1,,300)
rval <- 0.05
out1 <- improve_local(pilot1, newx, rval, sparameter=0.1^7, boundary="none")
out2 <- improve_local(pilot2, newx, rval, sparameter=0.1^7, boundary="none")
out3 <- improve_local(pilot3, newx, rval, sparameter=0.1^7, boundary="none")
out4 <- improve_local(pilot4, newx, rval, sparameter=0.1^7, boundary="none")


newf1 <- out1$predy[, ncol(out1$predy)]
newsigma1 <- out1$predsigma[length(out1$predsigma)]
newf2 <- out2$predy[, ncol(out2$predy)]
newsigma2 <- out2$predsigma[length(out2$predsigma)]
newf3 <- out3$predy[, ncol(out3$predy)]
newsigma3 <- out3$predsigma[length(out3$predsigma)]
newf4 <- out4$predy[, ncol(out4$predy)]
newsigma4 <- out4$predsigma[length(out4$predsigma)]

elrisk_final1 <- localrisk(pilot1, rmethod="true", pilotf=newf1 , pilotsigma=newsigma1)
elrisk_final2 <- localrisk(pilot2, rmethod="true", pilotf=newf2 , pilotsigma=newsigma2)
elrisk_final3 <- localrisk(pilot3, rmethod="true", pilotf=newf3 , pilotsigma=newsigma3)
elrisk_final4 <- localrisk(pilot4, rmethod="true", pilotf=newf4 , pilotsigma=newsigma4)

elambda_final1 <- movinglambda(elrisk_final1, newx, rval, boundary="none")
elambda_final2 <- movinglambda(elrisk_final2, newx, rval, boundary="none")
elambda_final3 <- movinglambda(elrisk_final3, newx, rval, boundary="none")
elambda_final4 <- movinglambda(elrisk_final4, newx, rval, boundary="none")


save(test1, test2, test3, test4, pilot1, pilot2, pilot3, pilot4, 
out1, out2, out3, out4, lambda.grid1, lambda.grid2, lambda.grid3, lambda.grid4,
elrisk_final1, elrisk_final2, elrisk_final3, elrisk_final4,
elambda_final1, elambda_final2, elambda_final3, elambda_final4, file="repara.RData")


newx <- seq(0,1,,300)
fite1 <- sreg(newx, -log(elambda_final1$lambda), lambda=0.1^4)
predlambda1 <- predict(fite1, newx)
dpredlambda1 <- predict(fite1, newx, derivative = 1)
fite2 <- sreg(newx, -log(elambda_final2$lambda), lambda=0.1^4)
predlambda2 <- predict(fite2, newx)
dpredlambda2 <- predict(fite2, newx, derivative = 1)
fite3 <- sreg(newx, -log(elambda_final3$lambda), lambda=0.1^4)
predlambda3 <- predict(fite3, newx)
dpredlambda3 <- predict(fite3, newx, derivative = 1)
fite4 <- sreg(newx, -log(elambda_final4$lambda), lambda=0.1^5)
predlambda4 <- predict(fite4, newx)
dpredlambda4 <- predict(fite4, newx, derivative = 1)


cairo_pdf(file="repara3.pdf", width=13.26, height=9.06, onefile = TRUE);par(font.main=2)
par(mfrow=c(3,2), cex.main=2.5, cex.axis=1.5, cex.lab=1.5, mar=c(3,5.5,4,2))
plot(test1$x[test1$x > 0 & test1$x < 1], test1$y[test1$x > 0 & test1$x < 1], main="(a)", xlab="", ylab="y")
lines(xgrid, fg1, col=2 ,lwd=2)
plot(test2$x[test1$x > 0 & test1$x < 1], test2$y[test1$x > 0 & test1$x < 1], main="(b)", xlab="", ylab="y")
lines(xgrid, fg2, col=2 ,lwd=2)
plot(test3$x[test1$x > 0 & test1$x < 1], test3$y[test1$x > 0 & test1$x < 1], main="(c)", xlab="", ylab="y")
lines(xgrid, fg3, col=2 ,lwd=2)
plot(test4$x[test1$x > 0 & test1$x < 1], test4$y[test1$x > 0 & test1$x < 1], main="(d)", xlab="", ylab="y")
lines(xgrid, fg4, col=2 ,lwd=2)

matplot(newx, cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda)), 
xlab="", ylab=expression(log(hat(lambda)(x))) , type="n", main="(e)")
lines(newx, -predlambda1, col=1, lwd=3)
lines(newx, -predlambda2, col=2, lwd=3, lty=2)
lines(newx, -predlambda3, col=3, lwd=3, lty=3)
lines(newx, -predlambda4, col=4, lwd=3, lty=4)
legend("bottomleft", legend=c("(a)", "(b)", "(c)", "(d)"), col=1:4, lty=1:4, lwd=2)

aa1 <- cbind(log(elambda_final1$lambda), log(elambda_final2$lambda), log(elambda_final3$lambda) , log(elambda_final4$lambda))
bb1 <- apply(aa1, 2, diff)/(newx[2] - newx[1])
dnewx <- newx[-length(newx)] + diff(newx)/2
ylim <- -range(c(range(dpredlambda1), range(dpredlambda2), range(dpredlambda3), range(dpredlambda4)))
matplot(dnewx, bb1, xlab="", ylab=expression(log^"'"~(hat(lambda)(x))) , type="n", main="(f)", ylim=sort(ylim))
lines(newx, -dpredlambda1, col=1, lwd=3)
lines(newx, -dpredlambda2, col=2, lwd=3, lty=2)
lines(newx, -dpredlambda3, col=3, lwd=3, lty=3)
lines(newx, -dpredlambda4, col=4, lwd=3, lty=4)
legend("topleft", legend=c("(a)", "(b)", "(c)", "(d)"), col=1:4, lty=1:4, lwd=2)
dev.off()