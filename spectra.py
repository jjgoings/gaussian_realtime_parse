import sys
import time
import numpy as np
from realtime import *

class Spectra(object):
    """Return a spectra object that can plot the total absorption spectrum or
    circular dichroism spectra for a given system. Can accept x,y,and z RT-TDDFT
    log files

    Attributes:
        spectra:        array. Contains the signal in the frequency domain.
        frequency:      array. Contains the frequencies considered.
        num_pts:        integer. Number of points in signal to consider for Pade
        x, y, z:        RealTime objects for x,y, and z pulses/runs
        spectra_type:   string. Type of spectra to plot, either 'abs' or 'ecd'
        damp_const:     float. Damp signal to give FWHM of (2/damp_const) in au.
        zero_pad:       integer. Pad signal with integer of zeros for aesthetics
        directions:     list of RT directions to consider for generating spectra
        tranformation:  string. Either 'fourier' or 'pade' for transformation.
    """
    def __init__(self,x=None,y=None,z=None,s='abs',d=150,zp=None,auto=False,
        num_pts=10000,trans='fourier'):
        self.spectra = None
        self.frequency = None

        self.num_pts = num_pts

        self.x = None
        self.y = None
        self.z = None

        self.spectra_type   = s
        self.damp_const     = d
        self.zero_pad       = zp
        self.transformation = trans

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
 
        # Do the isotropic transform to frequency domain 
        for q in self.directions:
            if self.transformation == 'pade':
                t0 = time.time()
                self.__dict__[q].pade_tx(q,self.spectra_type,self.damp_const,
                    self.num_pts)
                t1 = time.time()
                print "Pade done in: ", t1-t0
            elif self.transformation == 'fourier':
                self.__dict__[q].fourier_tx(q,self.spectra_type,\
                    self.damp_const,self.zero_pad,auto=auto)

        self.spectra = np.zeros_like(self.__dict__[self.directions[0]].fourier)
        for q in self.directions:
            self.frequency = self.__dict__[q].frequency
            self.spectra += self.__dict__[q].fourier


    def plot(self,xlim=[0,15],ylim=None,save=None,show=True,
        xlabel='Energy / eV',ylabel='S($\omega$) / arb. units',legend=None,
        grid=True,no_xticks=False,no_yticks=False,color='blue'):
        """ Plots the spectra you have obtained
            Variables:
                xlim:      list that defines range of x-axis, e.g. [xmin,xmax]
                ylim:      list that defines range of y-axis, e.g. [ymin,ymax]
                save:      string. filename of your saved output. you can change
                            the extension, e.g. save='output.pdf' will create 
                            pdf. default is png.
                show:      Boolean. Show plot interactively?
                xlabel:    string label for x-axis
                ylabel:    string label for y-axis
                legend:    string label for spectral line
                grid:      boolean. True means grid=on.
                no_xticks: boolean. True means turn xticks off
                no_yticks: boolean. True means turn yticks off
                color:     string. color of your spectral line
        """

        toEV = 27.2114 
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print "You need matplotlib to plot spectra"

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
        # Some defaults
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
        """ Return the peaks from the Fourier transform
            Variables:
            number:     integer. number of peaks to print.
            thresh:     float. Threshhold intensity for printing.

            Returns: Energy (eV), Intensity (depends on type of spectra)
        """
    
        from scipy.signal import argrelextrema as pks
        # find all peak indices [idx], and remove those below thresh [jdx]
        idx = pks(self.spectra,np.greater,order=5)
        jdx = np.where((np.abs(self.spectra[idx]) >= thresh))
        kdx = idx[0][jdx[0]] # indices of peaks matching criteria
        if number > len(kdx):
            number = len(kdx)
        print "First "+str(number)+" peaks (eV) found: "
        for i in xrange(number): 
            print "{0:.2f}".format(self.frequency[kdx][i]*27.2114),\
                  "{0:.2f}".format(self.spectra[kdx][i])
        
            

         


if __name__ == '__main__':
    cadmium = Spectra(x='cd')
    cadmium.peaks(9)
    cadmium.plot(xlim=[0,20],ylim=[-0.1,20],save='cd.pdf',show=False,
        no_yticks=True,color='green',legend='Cd')
    

