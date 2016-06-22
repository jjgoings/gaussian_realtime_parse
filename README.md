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


