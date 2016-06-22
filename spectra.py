import sys
import time
import numpy as np
from realtime import *

class Spectra(object):
    """Return a spectra object that can plot the total absorption spectrum or
    circular dichroism spectra for a given system. Can accept x,y,and z RT-TDDFT
    log files"""
    def __init__(self,x=None,y=None,z=None,s='abs',d=150,zp=None,auto=False,
        num_pts=10000):
        self.spectra = None
        self.frequency = None

        self.num_pts = num_pts

        self.x = None
        self.y = None
        self.z = None

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
            t0 = time.time()
            self.__dict__[q].pade_tx(q,self.spectra_type,self.damp_const,
                self.num_pts)
            t1 = time.time()
            print "Pade done in: ", t1-t0
            #self.__dict__[q].fourier_tx(q,self.spectra_type,self.damp_const,
            #    self.zero_pad,auto=auto)

        self.spectra = np.zeros_like(self.__dict__[self.directions[0]].fourier)
        for q in self.directions:
            self.frequency = self.__dict__[q].frequency
            self.spectra += self.__dict__[q].fourier


    def plot(self,xlim=[0,15],ylim=None,save=None,show=True,
        xlabel='Energy / eV',ylabel='S($\omega$) / arb. units',legend=None,
        grid=True,no_xticks=False,no_yticks=False,color='blue'):
        # save is a filename 
        toEV = 27.2114 
        import matplotlib.pyplot as plt
        ax = plt.subplot(111)
       
        ax.plot(self.frequency*toEV,self.spectra,label=legend,color=color)
        if legend:
            if isinstance(legend,str):
                plt.legend()
            else:
                raise ValueError, "'legend' needs to be a string"

        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if grid:
            plt.grid()
        if no_yticks:
            ax.set_yticklabels([]) 
        if no_xticks:
            ax.set_xticklabels([]) 
        ax.set_xlim(xlim)
        if not ylim:
            if self.spectra_type == 'abs':
                ax.set_ylim([0,4])
            elif self.spectra_type == 'ecd':
                ax.set_ylim([-400,400])
        else:
          ax.set_ylim(ylim)
        if save:
            if isinstance(save,str):
                plt.tight_layout()
                plt.savefig(save)
            else:
                raise ValueError, "'save' needs to be your filename!"
        if show:
            plt.show()
        plt.close()

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

    def peaks(self,number=3,thresh=0.01):
        """ Return the peaks from the Fourier transform"""
        from scipy.signal import argrelextrema as pks
        # find all peak indices [idx], and remove those below thresh [jdx]
        idx = pks(self.spectra,np.greater,order=5)
        jdx = np.where((self.spectra[idx] >= thresh))
        kdx = idx[0][jdx[0]] # indices of peaks matching criteria
        if number > len(kdx):
            number = len(kdx)
        print "First "+str(number)+" peaks (eV) found: "
        for i in xrange(number): 
            print "{0:.2f}".format(self.frequency[kdx][i]*27.2114),\
                  "{0:.2f}".format(self.spectra[kdx][i])
        
            

         


if __name__ == '__main__':
    N = 20000
    damp = 2000
    zinc     = Spectra(x='zn',num_pts=N,d=damp)
    cadmium  = Spectra(x='cd',num_pts=N,d=damp)
    mercury  = Spectra(x='hg',num_pts=N,d=damp)
    #spectra = Spectra(x='TlH-x2c-DFT_xx',
    #                  y='TlH-x2c-DFT_yy',
    #                  z='TlH-x2c-DFT_zz',
    #                  auto=True)
    #spectra = Spectra(x='AuH-x2c-DFT_xx',
    #                  y='AuH-x2c-DFT_yy',
    #                  z='AuH-x2c-DFT_zz',
    #                  auto=True)
    #spectra = Spectra(x='AuH-x2c_xx',
    #                  y='AuH-x2c_yy',
    #                  z='AuH-x2c_zz',
    #                  auto=True)
    #spectra = Spectra(x='cd')
    #spectra.z.test()
    zinc.peaks(9)
    zinc.plot(xlim=[0,20],ylim=[-0.1,20],save='zn.pdf',show=False,
        no_yticks=True,color='blue',legend='Zn')
    cadmium.peaks(9)
    cadmium.plot(xlim=[0,20],ylim=[-0.1,20],save='cd.pdf',show=False,
        no_yticks=True,color='green',legend='Cd')
    mercury.peaks(9)
    mercury.plot(xlim=[0,20],ylim=[-0.1,20],save='hg.pdf',show=False,
        no_yticks=True,color='red',legend='Hg')
    

