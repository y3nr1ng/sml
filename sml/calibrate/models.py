from abc import abstractmethod
from pprint import pprint

class BaseModel(object):
    __slots__ = ()

    def __init__(self, **kwargs):
        if not all(attr in kwargs for attr in self.__slots__):
            raise ValueError("incomplete parameter list")
        for key, value in kwargs.items():
            setattr(self, key, value)

    @abstractmethod
    def _equation(self, x):
        raise NotImplementedError

    def __call__(self, x):
        return _equation(x)

    def __str__(self):
        dump = '('
        dump += ', '.join(self.__slots__)
        dump += ')=('
        dump += ', '.join(str(getattr(self, attr)) for attr in self.__slots__)
        dump += ')'
        return dump

class PolynomialModel(BaseModel):
    """
    Describe the relationship between the axial position of a molecule in its
    imaged widths along two perpendicular axes by a pair of second degree
    polynomials.
    """
    __slots__ = ('w0', 'A', 'B', 'zc')

    def _equation(self, z):
        w0, A, B, zc = self.w0, self.A, self.B, self.zc
        return A * (z-zc)**2 + B

class HuangModel(BaseModel):
    """
    The relationship between the axial position of a molecule and its imaged
    widths along perpendicular axes is given by Huang et al.
    """
    __slots__ = ('w0', 'A', 'B', 'zc', 's')

    def _equation(self, z):
        #TODO suggest using (w0/2) to follow the original model
        w0, A, B, zc, s = self.w0, self.A, self.B, self.zc, self.s
        return w0 * np.sqrt( 1 + ((x-zc)/s)**2 * (1 + A*((x-zc)/s) + B*((x-zc)/s)**2) )
