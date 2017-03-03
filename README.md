## Gaussian RealTime Parse 
Gaussian RealTime Parse is a suite of Python scripts to ease the post-processing of Gaussian Real-Time TDDFT and Ehrenfest dynamics output.

You'll probably get a better idea of what you can do just opening up some of the scripts, but here is a basic workflow.

Say you have a RT-TDDFT output file called `rt.log`.

If you want to extract the z-dipole moment and plot it as a function of time. You can do so like this:

```
from realtime import RealTime
import matplotlib.pylot as plt

rt = RealTime('rt.log')
plt.plot(rt.time,rt.electricDipole.z)
plt.savefig('dipole.pdf')

```

This would plot the time-oscillating z-electric dipole moment (`electricDipole.z`) as a function of time (`time`) in au. Then it saves it to a PDF called `dipole.pdf`.

Pretty simple.

You can even plot spectra if you perturbed your system properly. Say you gave a Cadmium atom an x-type electric delta perturbation, and saved the output in `cd.log`.

Then this script should do the trick.

```
from spectra import Spectra

cd = Spectra(x='cd.log')
cd.plot(save='cd.pdf')

```

Now you're done! You have the fourier transformed x-dipole moment (the "spectra") saved to file `cd.pdf`.

# Debug tools for Gaussian L512

It is possible to run some sanity checks on the output of the latest version of Gaussian's realtime code.
We'll use the example `test_rabi.log` found in the `test` subdirectory. This is basically the Rabi test found 
in X. Li's 2005 paper "A time-dependent Hartree–Fock approach for studying the electronic optical response of molecules in intense fields."

In the sample script `debug.py`, you'll see we first load the output RT job.

```
rabi = RealTime('./test/test_rabi')
```

Sanity checks are simple, just do

```
rabi.test()
```
 This should test the applied field, what energy conservation (if any) exists,
that the iOps are called correctly, etc. The expected output is:

```
Energy conserved to:  4.82e-01  au
Max energy at time:   49.65  au
Min energy at time:   16.545  au
Time step             [OK]:  0.005  au
Number total steps    [OK]:  10000  steps
Field on:             [OK]
External field:       Linear
Ex field matches:     True
Ey field matches:     True
Ez field matches:     True
Bx field matches:     True
By field matches:     True
Bz field matches:     True
Field using:  Electric Dipole
Orthonormality        [OK]: Lowdin
Gauge                 [OK]: Length
```

So energy varies 0.482 Hartree and the time step matches the input. We know a field is on, the type of field is linear, and, importantly, the parsed applied field matches what we computed it should be internally. The test function will actually generate the field that was specified by the iOPs to compare. We also dump out what types of multipole moments went into our field definition, and check that we are doing the gauge we wanted.

For the specific Rabi case (this is not general) we can also dump out important HOMO and LUMO occupation information.
```
rabi.test_rabi()
```
This generates a CSV with the energy, appleid field, and HOMO/LUMO occupations as a function of time. You can plot the output then in Excel, etc. The CSV is called `rabi-analysis.csv`.

It looks like

```
#  Time (au), Ener. (au), Ez  (au), HOMO Occ, LUMO Occ
0.00000000,-0.58269663,1.00190378,1.00000000,0.00000000
0.00000000,-0.58269663,1.00190378,1.00000000,0.00000000
0.00500000,-0.58269667,1.00190378,1.00000000,0.00000000
0.01000000,-0.58269681,1.00190378,1.00000000,0.00000000
0.01500000,-0.58269703,1.00190378,1.00000000,0.00000000
0.02000000,-0.58269735,1.00190378,1.00000000,0.00000000
0.02500000,-0.58269775,1.00190378,1.00000000,0.00000000
0.03000000,-0.58269824,1.00190378,1.00000000,0.00000000
0.03500000,-0.58269883,1.00190378,1.00000000,0.00000000

...

```
 and so on.
