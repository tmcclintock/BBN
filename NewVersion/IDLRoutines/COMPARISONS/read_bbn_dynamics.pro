function read_bbn_dynamics, fdata, nsteps

openr,1,fdata
print,"Number of steps = ",nsteps

;Dynamics are in this arrays
;They appear as: Time, temperature, h, phie, and H
Dyn=dblarr(5,nsteps)

;temporary variable
tin=0.0D
T9in=0.0D
hin=0.0D
phiin=0.0D
ERin=0.0D
for j=0,nsteps-1 do begin

    readu,1,tin
    Dyn(0,j)=tin;the time

    readu,1,T9in
    Dyn(1,j)=T9in;the temperature
    ;print,T9in

    readu,1,hin
    Dyn(2,j)=hin;baryon density * T^3
    print,hin

    readu,1,phiin
    Dyn(3,j)=phiin;electron chemical potential

    readu,1,ERin
    Dyn(4,j)=ERin;expansion rate
    ;print,ERin
endfor                          ;nsteps

close,1
return,Dyn
end
