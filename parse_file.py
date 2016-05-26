import numpy as np

def parse_file(self):
    """Extract important attributes from the Gaussian realtime logfile."""
    filename = self.logfile
    lines = [line.rstrip('\n') for line in open(filename)] 

    muX = []
    muY = []
    muZ = []
    mX  = []
    mY  = []
    mZ  = []
    eX  = []
    eY  = []
    eZ  = []
    bX  = []
    bY  = []
    bZ  = []
    t   = []
    en  = []
   
    for idx, line in enumerate(lines):
        if line[1:26] == 'External field Parameters':
            self.envelope['Field']     = True
            self.envelope['Envelope']  = lines[idx+1].split()[2]   # string
            self.envelope['Gauge']     = lines[idx+2].split()[2]   # string 
            self.envelope['Ex']        = float(lines[idx+3].split()[2])  # au 
            self.envelope['Ey']        = float(lines[idx+4].split()[2])  # au
            self.envelope['Ez']        = float(lines[idx+5].split()[2])  # au
            self.envelope['Bx']        = float(lines[idx+6].split()[2])  # au
            self.envelope['By']        = float(lines[idx+7].split()[2])  # au
            self.envelope['Bz']        = float(lines[idx+8].split()[2])  # au
            self.envelope['Frequency'] = float(lines[idx+9].split()[2])  # au
            self.envelope['Phase']     = float(lines[idx+11].split()[2]) # au
            self.envelope['TOn']       = float(lines[idx+12].split()[2]) # au
            self.envelope['TOff']      = float(lines[idx+13].split()[2]) # au
            self.envelope['Terms']     = lines[idx+14].split()[3:] # mult str
        elif line[1:27] == 'No external field applied.':
            self.envelope['Field']     = False
        elif line[1:7] == 'Time =':
            time = line.split()
            t.append(float(time[2]))
        elif line[1:22] == 'Dipole Moment (Debye)': 
            dipole = lines[idx+1].split()
            muX.append(float(dipole[1])*0.393456) 
            muY.append(float(dipole[3])*0.393456) 
            muZ.append(float(dipole[5])*0.393456) 
        elif line[1:31] == 'Magnetic Dipole Moment (a.u.):': 
            dipole = lines[idx+1].split()
            mX.append(float(dipole[1])) 
            mY.append(float(dipole[3])) 
            mZ.append(float(dipole[5])) 
        elif line[1:9] == 'Energy =':
            energy = line.split()
            en.append(float(energy[2]))
        elif line[1:38] == 'Current electromagnetic field (a.u.):':
            efield = lines[idx+1].split()
            bfield = lines[idx+2].split()
            eX.append(float(efield[1])) 
            eY.append(float(efield[3])) 
            eZ.append(float(efield[5])) 
            bX.append(float(bfield[1])) 
            bY.append(float(bfield[3])) 
            bZ.append(float(bfield[5])) 
        elif line[1:27] ==  '        Restart MMUT every':
           self.mmut_restart = line.split()[3]

    # Save to object, if it exists
    if(muX and muY and muZ):
        self.electricDipole.x = np.asarray(muX)
        self.electricDipole.y = np.asarray(muY)
        self.electricDipole.z = np.asarray(muZ)
    if(mX and mY and mZ):
        self.magneticDipole.x = np.asarray(mX)
        self.magneticDipole.y = np.asarray(mY)
        self.magneticDipole.z = np.asarray(mZ)
    if(eX and eY and eZ):
        self.electricField.x  = np.asarray(eX)
        self.electricField.y  = np.asarray(eY)
        self.electricField.z  = np.asarray(eZ)
    if(bX and bY and bZ):
        self.magneticField.x  = np.asarray(bX)
        self.magneticField.y  = np.asarray(bY)
        self.magneticField.z  = np.asarray(bZ)
    if(t):
        self.time             = np.asarray(t)
    if(en):
        self.energy           = np.asarray(en)

def clean_data(self):
    """Make all the data arrays the same length, in case the log file
    did not finish a full time step (e.g. you killed the job early or are
    monitoring a job in progess. Furthermore, delete redundant time steps
    corresponding to when MMUT restarts"""
    
    def get_length(data):
        """Get length of array. If array is 'None', make it seem impossibly
        large"""
        if data.size:
            return len(data) 
        else:
            return 1e100
    # if doMMUT == True, we will delete duplicate data from MMUT restart
    doMMUT = True 
    lengths = []
    for x in self.propertyarrays:
       try:
           # If it is an array, remove MMUT steps, and grab its length 
           #FIXME Not sure if MMUT steps are actually double printed in latest
           if (doMMUT):
               self.__dict__[x] = np.delete(self.__dict__[x],
                   list(range(int(self.mmut_restart)-1, 
                   self.__dict__[x].shape[0], 
                   int(self.mmut_restart))), 
                   axis=0)
           lengths.append(get_length(self.__dict__[x]))
       except AttributeError:
           try:
               # Dipoles, fields, etc., are objects and we want their x/y/z
               for q in ['_x','_y','_z']:
                   #FIXME Again, not sure about MMUT duplicates
                   if (doMMUT):
                       self.__dict__[x].__dict__[q] = \
                           np.delete(self.__dict__[x].__dict__[q],
                           list(range(int(self.mmut_restart)-1, 
                           self.__dict__[x].__dict__[q].shape[0], 
                           int(self.mmut_restart))), 
                           axis=0)
                   lengths.append(get_length(self.__dict__[x].__dict__[q]))
           except:
               #print "Unknown data type: "+str(x)+str(q)
               pass

    self.min_length = min(lengths)
    # truncate all the arrays so they are the same length 
    truncate(self,self.min_length)

def truncate(self,length):
    """ Truncates the property arrays to a given *length* (integer) """
    for x in self.propertyarrays:
       try:
           # If it is an array, truncate its length 
           self.__dict__[x] = self.__dict__[x][:length]
       except TypeError:
           try:
               # Dipoles, fields, etc., are objects and we want their x/y/z
               for q in ['_x','_y','_z']:
                   self.__dict__[x].__dict__[q] = \
                       self.__dict__[x].__dict__[q][:length]
           except:
               #print "Unknown data type: "+str(x)+str(q)
               pass





