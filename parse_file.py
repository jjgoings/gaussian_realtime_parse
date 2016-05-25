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
        if line[1:7] == 'Time =':
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
    #TODO: replace attributes with self.__dict__ but remove 'logfile'
    attributes = ['electricDipole',
                  'magneticDipole',
                  'electricField',
                  'magneticField',
                  'time',
                  'energy']
    lengths = []
    for x in attributes:
        try:
            # If it is an array, remove MMUT steps, and grab its length 
           #FIXME Not sure if MMUT steps are actually double printed in latest
           # self.__dict__[x] = np.delete(self.__dict__[x],
           #     list(range(int(self.mmut_restart)-1, 
           #     self.__dict__[x].shape[0], 
           #     int(self.mmut_restart))), 
           #     axis=0)
            lengths.append(get_length(self.__dict__[x]))
        except AttributeError:
            try:
                # Dipoles, fields, etc., are objects and we want their x/y/z
                for q in ['_x','_y','_z']:
               #FIXME Again, not sure about MMUT duplicates
               #     self.__dict__[x].__dict__[q] = \
               #         np.delete(self.__dict__[x].__dict__[q],
               #         list(range(int(self.mmut_restart)-1, 
               #         self.__dict__[x].__dict__[q].shape[0], 
               #         int(self.mmut_restart))), 
               #         axis=0)
                    lengths.append(get_length(self.__dict__[x].__dict__[q]))
            except:
                print "Unknown data type: "+str(x)+str(q)

    min_length = min(lengths)

    for x in attributes:
        try:
            # If it is an array, truncate its length 
            self.__dict__[x] = self.__dict__[x][:min_length]
        except TypeError:
            try:
                # Dipoles, fields, etc., are objects and we want their x/y/z
                for q in ['_x','_y','_z']:
                    self.__dict__[x].__dict__[q] = \
                        self.__dict__[x].__dict__[q][:min_length]
            except:
                print "Unknown data type: "+str(x)+str(q)



