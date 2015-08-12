function read_bbn_abundances, fdata, nsteps, nnuc

openr,1,fdata
print,"Number of steps = ",nsteps
print,"Number of Nuclides = ",nnuc

;abundance array
Y=dblarr(2+nnuc,nsteps);time+temp+nnucs

;temporary variable
tin=0.0D
T9in=0.0D
Yin=0.0D
for j=0,nsteps-1 do begin
    readu,1,tin
    Y(0,j)=tin;the time
    readu,1,T9in
    Y(1,j)=T9in;the temperature
    for i=2,nnuc+2-1 do begin
        readu,1,Yin;the abundances
        Y(i,j)=Yin
        ;print,"Y_",i," = ",Y(i,j)
    endfor                      ;nnucs
endfor                          ;nsteps

close,1
return,Y
end
