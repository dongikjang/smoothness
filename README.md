=========
#R-script  for Nonparametric Estimation of Local Smoothness for Spatially Inhomogeneous Function


=========
### Preliminary Steps for Windows user
  * This R-script has been thoroughly tested under Mac OS X 10.10 and Ubuntu 14.04. If you have R for Windows version, you need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) and compile two Fortran codes.
  * On the install Rtool, You should check the "Edit the System Path" and add the directory which has `R.exe` file something like `C:\Program Files\R\R-3.2.x\bin\i386;` for 32-bit OS (if you use 64-bit os then the last folder may be `x64`). 
  * After install `Rtools`, please check the status of system path on the Command Prompt window (`CMD`). On the `CMD`, type `R` then `R` will be excuted with non-gui mode on the `CMD`.
 
=========
### Compile source codes 
 * For 64-bit OS
 
	``` 
R --arch=x86_64 CMD SHLIB -dynamiclib -O3 -fopenmp localrisk.f90 localfit.f90  -o localmethod64.dll 
	```
 * For 32-bit OS	
	```
R --arch=i386 CMD SHLIB -dynamiclib -O3 localrisk.f90 localfit.f90 -o localmethod32.dll
	```

 * On the R	
	```	
	path <- "where/the/compiled/dll/file"
  if(.Platform$r_arch == "x86_64"){
        dyn.load(paste(path, "localmethod64.dll", sep=""))
  }else{
        dyn.load(paste(path, "localmethod32.dll", sep=""))
  }
	```
	
=========
### An example
	```
source("cp_source.R")	
library(fields)

# generate a test set
set.seed(3245)
n <- 800
x <- sort(runif(n, 0,1))
xgrid <- seq(0, 1,,1000)
f4 <- ifelse(x < .5, sin((2*pi*(x+0.0026))^3), -sin((2*pi*(1-x-0.0026))^3))
fg4 <- ifelse(xgrid < .5, sin((2*pi*(xgrid+0.0026))^3), -sin((2*pi*(1-xgrid-0.0026))^3))

snr <- 5
sdf4 <- sd(f4)
y4 <- f4 + rnorm(length(f4), 0, sdf4/snr)

bdx <- c(-rev(x[x < 0.2]), x, rev(2-x[x > 0.8]))
f4 <- c(rev(f4[x < 0.2]), f4, rev(f4[x > 0.8]))
y4 <- c(rev(y4[x < 0.2]), y4, rev(y4[x > 0.8]))
x <- bdx
test4 <- list(x=x, y=y4, f=f4)


# pilot estimation 
pilot4 <- pilotQV(test4$x, test4$y)
lambda.grid4 <- pilot4$lambda.grid
nh <- length(lambda.grid1)

# update the pilot estimation
newx <- seq(0,1,,300)
rval <- 0.05
ptm <- proc.time()
out4 <- improve_local(pilot4, newx, rval, sparameter=0.1^7, boundary="none")
proc.time() - ptm
  ```
	
=========
### Created by Dongik Jang 06/20/2015 