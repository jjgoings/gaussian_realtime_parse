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
        propertyarrays:  List containing names of properties stored as arrays.
        truncate:        Method to truncate propertyarrays to a given length
        mmut_restart:    Integer containing how often MMUT restarts
    """

    def __init__(self, name):
        """Return a RealTime object whose logfile is *logfile*.""" 
        # Initialize data
        self.name           = name
        self.logfile        = name+'.log'
        self.envelope       = {}
        self.electricDipole = ElectricDipole()
        self.magneticDipole = MagneticDipole()
        self.electricField  = ElectricField()
        self.magneticField  = MagneticField()
        self.time           = None
        self.energy         = None
        self.frequency      = None
        self.fourier        = None
        # TODO May want to look at a better way of defining which attributes are
        # arrays instead of just hard-coding them in.
        self.propertyarrays = ['electricDipole',
                               'magneticDipole',
                               'electricField',
                               'magneticField',
                               'time',
                               'energy']
        self.truncate       = truncate
        self.min_length     = None
        self.mmut_restart   = 10000000000 # e.g. never restart

        # Call parser 
        parse_file(self)

        # Make all arrays consistent length
        clean_data(self)

    def fourier_tx(self,dipole_direction='x',spectra='abs',damp_const=150,
                    zero_pad=None,auto=False):
        """Return a set of frequencies and fourier transforms of a time
        dependent signal, e.g. return fourier transform of the x component of
        the time varying electric dipole"""
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

        if auto:
            dt = self.time[1] - self.time[0]
            damp_const = self.time[-1]/10.0
            #print "Damp const = ", damp_const
            print "Line width (eV) = ", \
                "{0:.3f}".format((2.0/damp_const)*27.2114)
             
            dipole = dipole - dipole[0]
            damp = np.exp(-(self.time-self.time[0])/float(damp_const))
            dipole = dipole * damp

            resolution = 0.025 #eV
            zero_pad   = int(np.floor((2.0*np.pi*27.2114)/(resolution*dt))\
                - len(self.time))
            if(zero_pad < 0.0):
                zero_pad = 0.0
            print "Number zeros = ", zero_pad

            zero = np.linspace(0,0,zero_pad)
            dipole = np.hstack((dipole,zero))

        else:
            dipole = dipole - dipole[0]
            damp = np.exp(-(self.time-self.time[0])/float(damp_const))
            dipole = dipole * damp

            if zero_pad:
                zero = np.linspace(0,0,zero_pad)
                dipole = np.hstack((dipole,zero))
    
        fw = fft(dipole)
        fw_re = np.real(fw)
        fw_im = np.imag(fw)
       
        n = len(fw_re)
        m = n / 2
        timestep = self.time[1] - self.time[0]
        self.frequency = fftfreq(n,d=timestep)*2.0*np.pi

        if spectra.lower() == 'abs':
            self.fourier = \
                -(4.0*self.frequency*np.pi*fw_im)/(3.0*137*kick_strength)
        elif spectra.lower() == 'ecd':
            self.fourier = \
                (17.32*fw_re)/(np.pi*kick_strength)

        # Grab positive values only
        self.frequency = self.frequency[1:m]
        self.fourier   = self.fourier[1:m]

    def test(self):
        self.check_energy()
        self.check_field()
        pass

    def check_energy(self):
        dE = abs(max(self.energy) - min(self.energy)) 
        t_maxE = self.time[np.argmax(self.energy)]
        t_minE = self.time[np.argmin(self.energy)]
        print "Energy conserved to: ", "{0:.2e}".format(dE), " au"
        print "Max energy at time: ", t_maxE, " au"
        print "Min energy at time: ", t_minE, " au"

    def check_field(self):
        if self.envelope['Field']:
            print "Checking the external field for: ", self.envelope['Envelope']
            print "Ex field matches: ", np.allclose(self.electricField.x,
                self.expected_field('Ex'))
            print "Ey field matches: ", np.allclose(self.electricField.y,
                self.expected_field('Ey'))
            print "Ez field matches: ", np.allclose(self.electricField.z,
                self.expected_field('Ez'))
            print "Bx field matches: ", np.allclose(self.magneticField.x,
                self.expected_field('Bx'))
            print "By field matches: ", np.allclose(self.magneticField.y,
                self.expected_field('By'))
            print "Bz field matches: ", np.allclose(self.magneticField.z,
                self.expected_field('Bz'))
        else:
            print "No external field applied"

    def expected_field(self,component):
        Time  = self.time
        TOn   = self.envelope['TOn']
        TOff  = self.envelope['TOff']
        Omega = self.envelope['Frequency']
        Phase = self.envelope['Phase']
        OmegT = Omega*(Time - TOn) + Phase
        if self.envelope['Envelope'] == 'Constant':
            field = np.zeros_like(self.time)
            # Step function, depending on how TOn and TOff are defined
            idx = np.where((self.time >= TOn) & (self.time < TOff))
            field[idx] = self.envelope[component]*np.cos(OmegT)
        else:
            print "Not a valid field!"
            sys.exit(0) 
        return field

             

      
        

    
 
if __name__ == '__main__':
    x = RealTime('h2')
    import matplotlib.pyplot as plt 
    plt.plot(x.time,x.electricField.z)
    plt.show()
    #x.test()
    
            

