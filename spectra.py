import sys
import numpy as np
from realtime import *

class Spectra(object):
    """Return a spectra object that can plot the total absorption spectrum or
    circular dichroism spectra for a given system. Can accept x,y,and z RT-TDDFT
    log files"""
    def __init__(self,x=None,y=None,z=None,s='abs',d=150,zp=None,auto=False):

        self.spectra_type = s
        self.damp_const   = d
        self.zero_pad     = zp

        # Load all the RealTime objects
        self.directions = []
        if x:
            self.x = RealTime(x)
            self.directions.append('x')
        if y:
            self.y = RealTime(y)
            self.directions.append('y')
        if z:
            self.z = RealTime(z)
            self.directions.append('z')

        # Enforce consistent data lengths
        self.align_data()
 
        # Do the isotropic fourier transform 
        for q in self.directions:
            self.__dict__[q].fourier_tx(q,self.spectra_type,self.damp_const,
                self.zero_pad,auto=auto)


    def plot(self,xlim=[0,15],ylim=None):
        toEV = 27.2114 
        import matplotlib.pyplot as plt
        ax = plt.subplot(111)
       
        S = np.zeros_like(self.__dict__[self.directions[0]].fourier)
        for q in self.directions:
            frequency = self.__dict__[q].frequency
            S += self.__dict__[q].fourier

        ax.plot(frequency*toEV,S)
        ax.set_xlim(xlim)
        if not ylim:
            if self.spectra_type == 'abs':
                ax.set_ylim([0,4])
            elif self.spectra_type == 'ecd':
                ax.set_ylim([-400,400])
        plt.show()

    def align_data(self):
        lengths = []
        if self.x:
            lengths.append(self.x.min_length) 
        if self.y:
            lengths.append(self.y.min_length) 
        if self.z:
            lengths.append(self.z.min_length) 
        min_length = min(lengths)
        if self.x:
            self.x.truncate(self.x,min_length)
        if self.y:
            self.y.truncate(self.y,min_length)
        if self.z:
            self.z.truncate(self.z,min_length)

         


if __name__ == '__main__':
    spectra = Spectra(x='AuH-x2c_xx',
                      y='AuH-x2c_yy',
                      z='AuH-x2c_zz',
                      s='abs',auto=True)
    spectra.plot()

