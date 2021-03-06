fdata="../../BinaryDataFiles/abundances.dat"
hdata="../../BinaryDataFiles/header.dat"

;header=dblarr(3)
header=read_bbn_header(hdata);,header)

nnuc=header(0)
nreacs=header(1)
nsteps=header(2)
print,"nnucs = ",nnuc
print,"nsteps = ",nsteps

Y=read_bbn_abundances(fdata,nsteps,nnuc)

set_plot,'ps';change device to postscript

;set the output file name to Abundances.eps, make it and encapsulated,
;8-bit depth color postscript
device,filename="Abundances.eps",/encapsulated,/color,bits=8

;set the x and y size to 5 inches
device,xsize=5.0,ysize=5.0,/inches

;plot t vs Yp
plot,/XLOG,/YLOG,Y(0,*),Y(2,*),$;xrange=[7,120],
yrange=[1e-20,1e2],xtitle="time (sec)",ytitle="abundances",$
  xstyle=1,ystyle=1,font=1

loadct,13,/silent ;load the rainbow color table silently
for i=3, nnuc+2-1 do begin
    oplot,Y(0,*),Y(i,*),color=4         ;over plot neutrons
endfor
loadct,0,/silent  ;load the black and white color table silently

device,/close ;close the postscript filt
set_plot,'x'  ;rever to x-window plotting
end
