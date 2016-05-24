class ElectricDipole(object):
    """An ElectricDipole object contains the x, y, and z electric dipole 
    components

    Attributes:
        x: x component of the electric dipole
        y: y component of the electric dipole
        z: z component of the electric dipole
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self,value):
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self,value):
        self._y = value
 
    @property
    def z(self):
        return self._z

    @z.setter
    def z(self,value):
        self._z = value

class MagneticDipole(object):
    """A MagneticDipole object contains the x, y, and z magnetic dipole 
    components

    Attributes:
        x: x component of the magnetic dipole
        y: y component of the magnetic dipole
        z: z component of the magnetic dipole
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self,value):
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self,value):
        self._y = value
 
    @property
    def z(self):
        return self._z

    @z.setter
    def z(self,value):
        self._z = value

class ElectricField(object):
    """An ElectricField object contains the x, y, and z electric field 
    components

    Attributes:
        x: x component of the electric field
        y: y component of the electric field
        z: z component of the electric field
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self,value):
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self,value):
        self._y = value
 
    @property
    def z(self):
        return self._z

    @z.setter
    def z(self,value):
        self._z = value

class MagneticField(object):
    """A MagneticField object contains the x, y, and z magnetic field 
    components

    Attributes:
        x: x component of the magnetic field
        y: y component of the magnetic field
        z: z component of the magnetic field
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self,value):
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self,value):
        self._y = value
 
    @property
    def z(self):
        return self._z

    @z.setter
    def z(self,value):
        self._z = value


