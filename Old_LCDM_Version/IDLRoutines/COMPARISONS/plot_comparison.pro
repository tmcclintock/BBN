;Tom McClintock
;Read in the header, then read in the dynamics, then plot the dynamics

fdata="dynamics.dat"
nsteps=1023
print,"Toms nsteps = ",nsteps
DynTom=read_bbn_dynamics(fdata,nsteps)

fdata="ABBNdynamics.dat"
nsteps=1591
print,"ABBN nsteps = ",nsteps
DynABBN=read_bbn_dynamics(fdata,nsteps)


set_plot,'ps';change device to postscript

;set the output file name to Abundances.eps, make it and encapsulated,
;8-bit depth color postscript
;set the x and y size to 5 inches
device,filename="DensityCompare.eps",/encapsulated,/color,bits=8,$
  xsize=5.0,ysize=5.0,/inches

loadct,13,/silent ;load the rainbow color table silently

;plot t vs Yp
plot,/XLOG,DynTom(0,*),DynTom(2,*),xrange=[0.01,10000],$
  ;yrange=[2.7e-5,1e-4],$
  title="Density Comparison",xtitle="Time (sec)",$
  ytitle="Baryon Density * T^3 (g*GK^3/cm^3)",$
  xstyle=1,ystyle=1,font=1,color=4

oplot,DynABBN(0,*),DynABBN(2,*),linestyle=4                   ;


device,/close ;close the postscript filt
set_plot,'x'  ;rever to x-window plotting
end
