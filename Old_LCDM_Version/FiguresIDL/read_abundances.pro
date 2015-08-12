function read_dynamics, time,Temp,h,phie

fdata = '../BinaryDataFiles/dynamics.dat' ;name of the abundance data file
ndata = '../BinaryDataFiles/header.dat'   ;name of the file that has the number of lines
openr,1,ndata

nnuc=0L
nreac=0L
nsteps=0L
readu,1,nnuc                  ;read in the number of lines
readu,1,nread
readu,1,nsteps
close, 1
openr,1,fdata
time=dblarr(nsteps)
Temp=dblarr(nsteps)
h   =dblarr(nsteps)
phie=dblarr(nsteps)

t_temp=0.0D
Temp_temp=0.0D
h_temp=0.0D
phie_temp=0.0D

;print,nsteps
I = 0L
for I=0,10 do begin
    readu,1,t_temp
    readu,1,Temp_temp
    readu,1,h_temp
    readu,1,phie_temp
    time(i)=t_temp
    Temp(i)=Temp_temp
    h(i)=h_temp
    phie(i)=phie_temp
endfor

close,1
;print, nsteps
return, 10

end

pro main
n=read_dynamics(a,b,c,d)
;print, n
i=0L
for i=0,n-1 do begin
    print, a(i)
    print, b(i)
;    print, c(i)
;    print, d(i)
endfor

end
