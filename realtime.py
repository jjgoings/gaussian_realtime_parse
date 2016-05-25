from __future__ import division
import numpy as np
import sys
from properties import *
from parse_file import *

class RealTime(object):
    """A RealTime object contains important parsed data from a Gaussian RealTime
    log file.

    Attributes:
        name:            A string containing primary filename
        logfile:         A string representing the Gaussian realtime log file 
        electricDipole:  Object containing x, y, z electric dipole moments (au)
        magneticDipole:  Object containing x, y, z magnetic dipole moments (au)
        electricField:   Object containing x, y, z electric field strengths (au)
        magneticField:   Object containing x, y, z magnetic field strengths (au)
        time:            Array containing time (au)
        energy:          Array containing energy (au)
        frequency:       Array containing frequencies from *time* (au)
        fourier:         Array containing fourier transformed signal (au)
    """

    def __init__(self, name):
        """Return a RealTime object whose logfile is *logfile*.""" 
        # Initialize data
        self.name           = name
        self.logfile        = name+'.log'
        self.electricDipole = ElectricDipole()
        self.magneticDipole = MagneticDipole()
        self.electricField  = ElectricField()
        self.magneticField  = MagneticField()
        self.time           = None
        self.energy         = None
        self.frequency      = None
        self.fourier        = None
        self.mmut_restart   = 10000000000 # e.g. never restart

        # Call parser 
        parse_file(self)

        # Make all arrays consistent length
        clean_data(self)

    def fft(self,dipole_direction='x',spectra='abs',damp_const=150):
        from scipy.fftpack import fft, fftfreq 
        # Choose which signal to FFT
        if spectra.lower() == 'abs':  
            if dipole_direction.lower() == 'x':
                dipole = self.electricDipole.x
                kick_strength = self.electricField.x[0]
            elif dipole_direction.lower() == 'y':
                dipole = self.electricDipole.y
                kick_strength = self.electricField.y[0]
            elif dipole_direction.lower() == 'z':
                dipole = self.electricDipole.z
                kick_strength = self.electricField.z[0]
            else:
                print "Not a valid direction for the dipole! Try: x,y,z "
        elif spectra.lower() == 'ecd':
            if dipole_direction.lower() == 'x':
                dipole = self.magneticDipole.x
                kick_strength = self.electricField.x[0]
            elif dipole_direction.lower() == 'y':
                dipole = self.magneticDipole.y
                kick_strength = self.electricField.y[0]
            elif dipole_direction.lower() == 'z':
                dipole = self.magneticDipole.z
                kick_strength = self.electricField.z[0]
            else:
                print "Not a valid direction for the dipole! Try: x,y,z "
        else: 
            print "Not a valid spectra choice"

        if np.isclose(kick_strength,0.0):
            print "Kick = 0, you are trying to FFT the wrong file"
            print "Try to change dipole direction!"
            sys.exit(0)

        dipole = dipole - dipole[0]
        damp = np.exp(-(self.time-self.time[0])/float(damp_const))
        dipole = dipole * damp
        

        fw = fft(dipole)
        fw_re = np.real(fw)
        fw_im = np.imag(fw)
       
        n = len(fw_re)
        timestep = self.time[1] - self.time[0]
        self.frequency = fftfreq(n,d=timestep)*2.0*np.pi
        
        if spectra.lower() == 'abs':
            self.fourier = \
                -(4.0*self.frequency*np.pi*fw_im)/(3.0*137*kick_strength)
        elif spectra.lower() == 'ecd':
            self.fourier = \
                (17.32*fw_re)/(np.pi*kick_strength)
        

    
 
if __name__ == '__main__':
    x = RealTime('he')
    x.fft(dipole_direction='x',damp_const=800)
    import matplotlib.pyplot as plt
    #plt.plot(x.time,x.electricDipole.x)
    #plt.plot(x.time,x.energy)
    plt.plot(x.frequency*27.2114,x.fourier)
    plt.xlim(60,80)
    #plt.savefig('he.pdf')
    plt.show()
    
            

