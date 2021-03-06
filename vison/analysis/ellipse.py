

"""

Auxiliary module with functions to generate generalized ellipse masks.

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import collections
import numpy as np
from scipy.special import gamma
import unittest
from functools import reduce
# END IMPORT


def dist_superellipse(n, center, q=1., pos_ang=0., c=0.):
    """Form an array in which the value of each element is equal to the
    semi-major axis of the superellipse of specified center, axial ratio,
    position  angle, and c parameter which passes through that element.
    Useful for super-elliptical aperture photometry.

    Inspired on dist_ellipse.pro from AstroLib (IDL).

    Note: this program doesn't take into account the change in the order of axes
    from IDL to Python. That means, that in 'n' and in 'center', the order of the
    coordinates must be reversed with respect to the case for dist_ellipse.pro, in
    order to get expected results. Nonetheless, the polar angle means the
    counter-clock wise angle with respect to the 'y' axis.

    :param n: shape of array (N1,N2), it can be an integer (squared shape NxN)
    :param center: center of superellipse radii: (c1,c2)
    :param q: axis ratio r2/r1
    :param pos_ang: position angle of isophotes, in degrees, CCW from axis 1
    :param c: boxyness (c>0) /diskyness (c<0)


    """

    # CHECK INPUTS

    criterio1 = isinstance(n, int) or \
        (isinstance(n, collections.Sequence) and len(n) == 2)

    assert criterio1, 'n must be an scalar or a 2 elements tuple'

    if not isinstance(n, collections.Sequence):
        n = (n, n)
    else:
        n = tuple([int(x) for x in n])

    criterio2 = isinstance(q, (int, float))

    if not criterio2:
        print('q must be an integer or float')
        criterio3 = False
    else:
        criterio3 = q >= 0 and q <= 1

    assert criterio3, 'q = %f is out of range (0,1)' % q

    if criterio2 and criterio3:
        q = float(q)

    criterio4 = isinstance(pos_ang, (int, float))
    if not criterio4:
        print('pos_ang must be an integer or float')
        criterio5 = 'False'
    else:
        criterio5 = -180 <= pos_ang <= 180

    assert criterio5, 'pos_ang out of range (-180,180)'

    if criterio4 and criterio5:
        pos_ang = float(pos_ang)

    criterio6 = isinstance(c, (int, float))
    if not criterio6:
        print('c must be an integer or float')
    else:
        c = float(c)

    criterio7 = isinstance(center, collections.Sequence) and len(center) == 2
    if not criterio7:
        print('center must be a 2 element tuple')
    else:
        center = tuple([float(x) for x in center])

    # END CHECK OF INPUTS

    criterios = reduce(np.logical_and, np.array([criterio1, criterio2, criterio3,
                                                 criterio4, criterio5, criterio6, criterio7]))

    if not criterios:
        return None

    radeg = 180. / np.pi

    ang = pos_ang / radeg
    cosang = np.cos(ang)
    sinang = np.sin(ang)

    if len(n) == 2:
        nx = n[1]
        ny = n[0]
    else:
        nx = ny = n

    x = np.arange(nx, dtype='Float32') - center[1]
    y = np.arange(ny, dtype='Float32') - center[0]
    im = np.zeros(shape=(ny, nx), dtype='Float32')
    xcosang = x * cosang
    xsinang = x * sinang

    for i in range(ny):
        xtemp = xcosang + y[i] * sinang
        ytemp = -xsinang + y[i] * cosang
        im[i, :] = ((np.abs(xtemp / q))**(c + 2.) +
                    (np.abs(ytemp))**(c + 2.))**(1. / (c + 2.))

    return im


def area_superellip(r, q, c=0):
    """Returns area of superellipse, given the semi-major axis length"""

    a_dummie = 4.**(1. - (c + 2)**(-1.))
    b_dummie = r**2. * q * ((np.pi)**(0.5))
    c_dummie = gamma(1. + (c + 2.)**(-1.)) / gamma(0.5 + (c + 2.)**(-1.))
    area = a_dummie * b_dummie * c_dummie
    return area


def effective_radius(area, q=1., c=0.):
    """Returns semi-major axis length of superellipse, given the area"""

    a_dummie = 4.**(1. - (c + 2.)**(-1.))
    c_dummie = gamma(1. + (c + 2.)**(-1.)) / gamma(0.5 + (c + 2.)**(-1.))

    b_dummie = area / (a_dummie * c_dummie)
    r = (b_dummie / (q * (np.pi ** 0.5)))**(0.5)

    return r


class TestEllipse(unittest.TestCase):
    """
    Unit tests for the ellipse module.
    """

    def setUp(self):

        self.tolerance = 1.e-7
        self.n = 5
        self.center = (2., 2.)
        self.q = 1.
        self.pos_ang = 0.
        self.c = 0.
        self.actual_ellipse = np.array([[2.82842708,
                                         2.23606801,
                                         2.,
                                         2.23606801,
                                         2.82842708],
                                        [2.23606801,
                                         1.41421354,
                                         1.,
                                         1.41421354,
                                         2.23606801],
                                        [2.,
                                         1.,
                                         0.,
                                         1.,
                                         2.],
                                        [2.23606801,
                                         1.41421354,
                                         1.,
                                         1.41421354,
                                         2.23606801],
                                        [2.82842708,
                                         2.23606801,
                                         2.,
                                         2.23606801,
                                         2.82842708]],
                                       dtype='float32')

    def test_dse_squared(self):

        n = self.n
        ret = dist_superellipse((n, n), self.center,
                                self.q, self.pos_ang, self.c)
        ans = ((ret - self.actual_ellipse)**2.).sum()
        self.assertAlmostEqual(0., ans, msg='expected=%f, got=%f' % (0., ans),
                               delta=self.tolerance)

    def test_area_superellip(self):
        r = 1.
        ans = area_superellip(r, q=1., c=0.)
        self.assertAlmostEqual(np.pi, ans, msg='expected=%f, got=%f' % (np.pi, ans),
                               delta=self.tolerance)

    def test_effective_radius(self):

        area = np.pi
        ans = effective_radius(area, q=1., c=0.)
        self.assertAlmostEqual(1., ans, msg='expected=%f, got=%f' % (1., ans),
                               delta=self.tolerance)


if __name__ == '__main__':
    # testing section
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEllipse)
    unittest.TextTestRunner(verbosity=3).run(suite)
