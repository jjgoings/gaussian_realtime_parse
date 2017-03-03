%cpu=2
%mem=2GB
%chk=water_RHF.chk
#p hf/sto-3g complex nosymm 

water ground state

0 1
O                   0.000000   -0.075792    0.000000
H                   0.866812    0.601436    0.000000
H                  -0.866812    0.601436    0.000000


--Link1--
%cpu=2
%mem=2GB
%oldchk=water_RHF.chk
%chk=water_delta_RHF.chk
%subst l512 .
%subst l118 .
#p hf geom=check chkbasis Ehrenfest(Electroniconly) complex pop=always guess=read nosymm 
iop(5/177=200,5/134=500,5/138=1,5/140=-1,5/141=1)
extralinks(308)

Water with z-delta kick

0 1

0

1 0 0.0 0.0 0.02 0.0 0.0 0.0 0.0 0.0 0.0 0.00001

