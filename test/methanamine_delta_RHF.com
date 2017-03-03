%cpu=2
%mem=2GB
%chk=methanamine.chk
#p hf/sto-3g complex nosymm 

CH2NH2+ ground state

1 1
 C,0,0.0000003838,-0.6887902204,0.
 H,0,-0.954473116,-1.2422749536,0.
 H,0,0.9545045493,-1.2422327762,0.
 N,0,0.0000003838,0.6034511167,0.
 H,0,0.8814678965,1.150635869,0.
 H,0,-0.8815092313,1.1506178278,0.


--Link1--
%cpu=2
%mem=2GB
%oldchk=methanamine.chk
%chk=methanamine_delta_RHF.chk
%subst l512 .
%subst l118 .
#p hf geom=check chkbasis Ehrenfest(Electroniconly) complex pop=always guess=read nosymm 
iop(5/177=200,5/134=500,5/138=1,5/140=-1,5/141=1)
extralinks(308)

CH2NH2+ activated for torsional for initial density for Ehrenfest

1 1

0

1 0 0.0000 0.0 0.02 0.0 0.0 0.0 0.0 0.0 0.0 0.00001

