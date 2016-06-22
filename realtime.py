from __future__ import division
import numpy as np
from scipy.linalg import solve
import sys
import time
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
        iops:            Dict containing IOps for 512 
        envelope:        Dict containing field parameters printed in logfile
        time:            Array containing time (au)
        energy:          Array containing energy (au)
        frequency:       Array containing frequencies from *time* (au)
        fourier:         Array containing fourier transformed signal (au)
        propertyarrays:  List containing names of properties stored as arrays.
        truncate:        Method to truncate propertyarrays to a given length
        mmut_restart:    Integer containing how often MMUT restarts
        au2fs:           Scalar constant to convert au to femtoseconds 
    """

    def __init__(self, name):
        """Return a RealTime object whose logfile is *logfile*.""" 
        # Initialize data
        self.name           = name
        self.logfile        = name+'.log'
        self.envelope       = {}
        self.iops           = {'132':['0'],
                               '134':['0'], 
                               '177':['0'], 
                               '136':['0'], 
                               '137':['0'], 
                               '138':['0'], 
                               '139':['0'], 
                               '140':['0'], 
                               '141':['0'], 
                               '142':['0'], 
                               '143':['0'], 
                               '144':['0']}
        self.electricDipole = ElectricDipole()
        self.magneticDipole = MagneticDipole()
        self.electricField  = ElectricField()
        self.magneticField  = MagneticField()
        self.orthonorm      = None
        self.step_size      = None
        self.total_steps    = None
        self.time           = None
        self.energy         = None
        self.frequency      = None
        self.fourier        = None
        self.au2fs          = 0.0241888425
        # TODO May want to look at a better way of defining which attributes are
        # arrays instead of just hard-coding them in.
        self.propertyarrays = ['electricDipole',
                               'magneticDipole',
                               'electricField',
                               'magneticField',
                               'time',
                               #FIXME valid for H2+ Rabi ONLY
                               'HOMO',
                               'LUMO',
                               'energy']
        self.truncate       = truncate
        self.min_length     = None
        self.mmut_restart   = 10000000000 # e.g. never restart
        #FIXME: ONLY FOR H2+ RABI
        self.HOMO           = None
        self.LUMO           = None

        # Call parser 
        parse_file(self)
        decode_iops(self)

        # Make all arrays consistent length
        clean_data(self)

    def pade_tx(self,dipole_direction='x',spectra='abs',damp_const=5500,
        num_pts=10000):
        # num_pts: number of points to sample for pade transformation

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
 

        # skip is integer to skip every n-th value
        # skip = 1 would not skip any values, but skip = 10 would only
        # consider every tenth value
        skip = 1 
        dipole = dipole - dipole[0]
        dipole = dipole[::skip]
        damp = np.exp(-(self.time-self.time[0])/float(damp_const))
        damp = damp[::skip]
        dipole = dipole * damp

        timestep = skip*(self.time[2] - self.time[1])
        M = len(dipole)
        N = int(np.floor(M / 2))

        if N > num_pts:
            N = num_pts

        # G and d are (N-1) x (N-1)
        # d[k] = -dipole[N+k] for k in range(1,N)
        d = -dipole[N+1:2*N] 

        # Old code, which works with regular Ax=b linear solver. 
        # G[k,m] = dipole[N - m + k] for m,k in range(1,N)
        #G = dipole[N + np.arange(1,N)[:,None] - np.arange(1,N)]
        #b = solve(G,d,check_finite=False)

        # Toeplitz linear solver using Levinson recursion
        # Should be O(n^2), and seems to work well, but if you get strange
        # results you may want to switch to regular linear solver which is much
        # more stable.
        try:
            from scipy.linalg import toeplitz, solve_toeplitz
        except ImportError:
            print "You'll need SciPy version >= 0.17.0"
            
        # Instead, form G = (c,r) as toeplitz
        #c = dipole[N:2*N-1]
        #r = np.hstack((dipole[1],dipole[N-1:1:-1]))
        b = solve_toeplitz((dipole[N:2*N-1],\
            np.hstack((dipole[1],dipole[N-1:1:-1]))),d,check_finite=False)
      
        # Now make b Nx1 where b0 = 1 
        b = np.hstack((1,b)) 

        # b[m]*dipole[k-m] for k in range(0,N), for m in range(k) 
        a = np.dot(np.tril(toeplitz(dipole[0:N])),b)

        p = np.poly1d(a)
        q = np.poly1d(b)

        # If you want energies greater than 2*27.2114 eV, you'll need to change
        # the default frequency range to something greater.
        self.frequency = np.arange(0,2,0.0001)
        W = np.exp(-1j*self.frequency*timestep)

        fw_re = np.real(p(W)/q(W))
        fw_im = np.imag(p(W)/q(W))

        if np.any(np.isinf(self.frequency)) or np.any(np.isnan(self.frequency)):
            print "Check your dT: frequency contains NaNs and/or Infs!"
            sys.exit(0)

        if spectra.lower() == 'abs':
            self.fourier = \
                (4.0*self.frequency*np.pi*fw_im)/(3.0*137*kick_strength)
        elif spectra.lower() == 'ecd':
            self.fourier = \
                (17.32*fw_re)/(np.pi*kick_strength)

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
            dt = self.time[2] - self.time[1]
            damp_const = self.time[-1]/10.0
            line_width = (2.0/damp_const)*27.2114
            #print "Damp const = ", damp_const
            if line_width > 2.0:
                print "Large line width: ", "{0:.3f}".format(line_width)," eV"
                print "Spectra not meaningful. Exiting..."
                sys.exit(0)
            else:
                print "Line width (eV) = ", "{0:.3f}".format(line_width)
             
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
        m = int(n / 2)
        timestep = self.time[2] - self.time[1]
        self.frequency = fftfreq(n,d=timestep)*2.0*np.pi
        if np.any(np.isinf(self.frequency)) or np.any(np.isnan(self.frequency)):
            print "Check your dT: frequency contains NaNs and/or Infs!"
            sys.exit(0)

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
        self.check_iops()
        pass

    def check_energy(self):
        dE = abs(max(self.energy) - min(self.energy)) 
        t_maxE = self.time[np.argmax(self.energy)]
        t_minE = self.time[np.argmin(self.energy)]
        print "Energy conserved to: ", "{0:.2e}".format(dE), " au"
        print "Max energy at time:  ", t_maxE, " au"
        print "Min energy at time:  ", t_minE, " au"

    def check_field(self,tol=1e-6):
        if self.envelope['Field']:
            print "External field:      ", self.envelope['Envelope']
            print "Ex field matches:    ", np.allclose(self.electricField.x,
                self.expected_field('Ex'),atol=tol)
            print "Ey field matches:    ", np.allclose(self.electricField.y,
                self.expected_field('Ey'),atol=tol)
            print "Ez field matches:    ", np.allclose(self.electricField.z,
                self.expected_field('Ez'),atol=tol)
           # print "Bx field matches: ", np.allclose(self.magneticField.x,
           #     self.expected_field('Bx'),atol=tol)
           # print "By field matches: ", np.allclose(self.magneticField.y,
           #     self.expected_field('By'),atol=tol)
           # print "Bz field matches: ", np.allclose(self.magneticField.z,
           #     self.expected_field('Bz'),atol=tol)
        else:
            print "No external field applied"

    def check_iops(self):
        """ Check internal consistency of some set iops and values printed out
        to the logfile, as well as some derived quantities"""
        # Check the step size
        if self.step_size == (self.time[2] - self.time[1]):
            if ((self.step_size == 0.05) \
                and (int(self.iops['134'][0]) == 0)) or\
               (self.step_size == float(self.iops['134'][0])*0.00001):
                print "Time step             [OK]: ", self.step_size, " au"
        else:
            print "Inconsistent time step: "
            print "  IOps:                  ", self.iops['134'][1]
            print "  logfile header showing ", self.step_size
            print "  logfile showing        ", self.time[2] - self.time[1]
        # Check the total propagation steps
        if ((self.total_steps == 50) \
           and (int(self.iops['177'][0]) == 0)) or\
          (self.total_steps == int(self.iops['177'][0])):
               print "Total propagation     [OK]: ", self.total_steps, " steps"
        else:
            print "Inconsistent propagation time: "
            print "  IOps:                  ", self.iops['177'][1]
            print "  logfile header showing ", self.total_steps
        # Check if external field is indeed On or OFF
        if ((self.envelope['Field'] == False) and\
           (int(self.iops['138'][0]) == 0)):
            print "Field off:            [OK]"
        elif (self.envelope and int(self.iops['138'][0]) != 0):
            print "Field on:             [OK]"
            self.check_field()
        else:
            print "Inconsistency in field:"
            print "IOps:                     ", self.iops['138'] 
        
        # Check Orthonormalization
        if ((self.orthonorm == self.iops['136'][1])):
            print "Orthonormality        [OK]:", self.orthonorm
        else:
            print "Inconsistency in orthonormality"
            print "IOps:                      ", self.iops['136'][1]
            print "logfile showing:           ", self.iops['136'][1]

    def expected_field(self,component):
        Time  = self.time
        TOn   = self.envelope['TOn']
        TOff  = self.envelope['TOff']
        Omega = self.envelope['Frequency']
        Phase = self.envelope['Phase']
        OmegT = Omega*(Time - TOn) + Phase
        field = np.zeros_like(self.time)
        if self.envelope['Envelope'] == 'Constant':
            # Step function, depending on how TOn and TOff are defined
            idx = np.where((Time >= TOn) & (Time < TOff))
            # in GDV OmegT begins at TOn as well
            field[idx] = self.envelope[component]*np.cos(OmegT[idx])
        elif self.envelope['Envelope'] == 'Linear':
            TMax = (2.0*np.pi)/Omega
            # Linearly ramp off to zero 
            idx = np.where((Time >= TOn) & (Time <= TOff) & \
                (Time > TOff-TMax))
            field[idx] = self.envelope[component]*\
                         ((TOff-Time[idx])/TMax)*np.cos(OmegT[idx])
            # Constant envelope 
            idx = np.where((Time >= TOn) & (Time <= TOff) & \
                (Time > TOn+TMax) & (Time <= TOff-TMax))
            field[idx] = self.envelope[component]*np.cos(OmegT[idx])
            # Linearly ramp up to maximum in first cycle
            idx = np.where((Time >= TOn) & (Time <= TOff) & \
                (Time <= TOn+TMax))
            field[idx] = self.envelope[component]*\
                         ((Time[idx]-TOn)/TMax)*np.cos(OmegT[idx])
        elif self.envelope['Envelope'] == 'Gaussian':
            idx = np.where((Time >= TOn) & (Time < TOff))
            #FIXME: Sigma is hard-coded for testing...need to print it in the 
            # output and then search for it during parsing.
            Sigma = 0.01
            TCntr = np.sqrt(np.log(1000.0))/Sigma
            field[idx] = self.envelope[component]*\
                         np.cos(OmegT[idx])*\
                         np.exp(-(Sigma*(Time[idx]-TCntr))**2)
        else:
            print "Not a valid field!"
            sys.exit(0) 
        return field

 
if __name__ == '__main__':
    a = RealTime('test')
    import matplotlib.pyplot as plt 
    plt.plot(a.time,a.electricDipole.z)
    plt.savefig('dipole.pdf')
    #plt.show()
    
            

