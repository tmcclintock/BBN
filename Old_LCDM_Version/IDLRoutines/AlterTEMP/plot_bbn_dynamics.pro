;Tom McClintock
;Read in the header, then read in the dynamics, then plot the dynamics

;fdata="../../BinaryDataFiles/dynamics.dat"
;hdata="../../BinaryDataFiles/header.dat"
fdata="ABBNdynamics.dat"

;header=dblarr(3)
;header=read_bbn_header(hdata);,header)

;nnuc=header(0)
;nreacs=header(1)
;nsteps=header(2)
nsteps=1591;read off from within AlterBBN
print,"nsteps = ",nsteps

Dyn=read_bbn_dynamics(fdata,nsteps)

set_plot,'ps';change device to postscript

;set the output file name to Abundances.eps, make it and encapsulated,
;8-bit depth color postscript
;set the x and y size to 5 inches
device,filename="ABBNTemp.eps",/encapsulated,/color,bits=8,$
  xsize=5.0,ysize=5.0,/inches

loadct,13,/silent ;load the rainbow color table silently

;plot t vs Yp
plot,/XLOG,/YLOG,Dyn(0,*),Dyn(1,*),yrange=[0.1,100],xrange=[0.01,10000],$
title="ABBN Temperature",xtitle="Time (sec)",ytitle="Temperature (GK)",$
  xstyle=1,ystyle=1,font=1,color=4

device,/close ;close the postscript filt
set_plot,'x'  ;rever to x-window plotting
end
