function read_bbn_header,hdata;,header

openr,1,hdata
nnucs=0L
nreacs=0L
nsteps=0L

readu,1,nnucs
readu,1,nreacs
readu,1,nsteps

header=dblarr(3)
header(0)=nnucs
header(1)=nreacs
header(2)=nsteps

close,1
return,header
end
