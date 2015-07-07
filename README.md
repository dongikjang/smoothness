=========
#R-script  for Nonparametric Estimation of Local Smoothness for Spatially Inhomogeneous Function


=========
### Preliminary Steps for Windows user
  * This R-script has been thoroughly tested under Mac OS X 10.10 and Ubuntu 14.04. If you have R for Windows version, you need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) and compile two Fortran codes.
  * On the install Rtool, You should check the "Edit the System Path" and add the directory which has `R.exe` file something like `C:\Program Files\R\R-3.2.x\bin\i386;` for 32-bit OS (if you use 64-bit os then the last folder may be `x64`). 
  * After install `Rtools`, please check the status of system path on the Command Prompt window (`CMD`). On the `CMD`, type `R` then `R` will be excuted with non-gui mode on the `CMD`.
  
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
### Created by Dongik Jang 06/20/2015 