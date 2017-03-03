%cpu=2
%mem=2GB
%oldchk=test_rabi.chk
%subst l512 ..
%subst l118 ..
#p hf/sto-3g Ehrenfest(Electroniconly,restart) complex pop=always nosymm 
iop(5/177=10000,5/134=500,5/138=1,5/140=0,5/141=1,4/13=-2)
extralinks(308)

H2+ Rabi oscillation

1 2
H
H 1 B1

B1 1.0603

0

1 1 0.0 0.0 0.05 0.0 0.0 0.0 0.4746 1.57079632679 0.0 1000000.0


