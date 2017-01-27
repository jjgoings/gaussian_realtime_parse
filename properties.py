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

class ElectricQuadrupole(object):
    """An ElectricQuadrupole object contains the xx, xy, xz, yy, yz and zz 
    unique electric quadrupole components

    Attributes:
        xx: xx component of the electric quadrupole
        xy: xy component of the electric quadrupole
        xz: xz component of the electric quadrupole
        yy: yy component of the electric quadrupole
        yz: yz component of the electric quadrupole
        zz: zz component of the electric quadrupole
    """
    def __init__(self):
        self._xx = None
        self._xy = None
        self._xz = None
        self._yy = None
        self._yz = None
        self._zz = None

    @property
    def xx(self):
        return self._xx

    @xx.setter
    def xx(self,value):
        self._xx = value

    @property
    def xy(self):
        return self._xy

    @xy.setter
    def xy(self,value):
        self._xy = value


    @property
    def xz(self):
        return self._xz

    @xz.setter
    def xz(self,value):
        self._xz = value


    @property
    def yy(self):
        return self._yy

    @yy.setter
    def yy(self,value):
        self._yy = value

    @property
    def yz(self):
        return self._yz

    @yz.setter
    def yz(self,value):
        self._yz = value

    @property
    def zz(self):
        return self._zz

    @zz.setter
    def zz(self,value):
        self._zz = value

