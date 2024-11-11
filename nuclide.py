from fractions import Fraction
from os.path import realpath
from sqlite3 import connect

if __name__ == '__main__' or __name__ == 'nuclide':
    from physics_value import PhysicsValue
    from constant import c, u
else:
    from .physics_value import PhysicsValue
    from .constant import c, u


__all__ = ['DataManager', 'Nuclide']
PY_PATH = realpath(__file__)
DATABASE_PATH = PY_PATH.split('.')[0] + '.db'


ELEMENTS = ['n', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
            'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
            'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
INVERT_ELEMENTS: dict[str, int] = {}
for i, n in enumerate(ELEMENTS):
    INVERT_ELEMENTS[n] = i

ISOMERS = {0: '', 1: 'm', 2: 'n', 3: 'p', 4: 'q', 5: 'r', 8: 'i', 9: 'j'}
ISO_NUM = {'': 0, 'm': 1, 'n': 2, 'p': 3, 'q': 4, 'r': 5, 'i': 8, 'j': 9}


PROTON_ENERGY = PhysicsValue(7288.971064, 'keV')
NEUTRON_ENERGY = PhysicsValue(8071.3181, 'keV')


class UnknownElement(ValueError):
    """
    raised when meet an illegal element symbol or illegal atomic number.
    """


class UnknownNuclide(ValueError):
    """
    raised when a nuclide is not found in NuBase 2020.
    """


class Nuclide:
    """
    The Nuclide object stores basic data for a nuclide, with all non-dimensionless values presented
    in the 'PhysicsValue' type accompanied by appropriate units.

    @property z: atomic number;
    @property a: mass number;
    @property m: isomer id, for ground state m = 0;
    @property name: the nuclide symbol, in format like 'Ag-110m';
    @property mass: the rest mass of the nuclide, in mass dimension;
    @property mass_excess: mass excess in energy dimension;
    @property binding_energy and @property rest_energy: same as which in physics;
    @property half_life: half life in time dimension, positive with two special exceptions:
                         -1s for stable particle and 0 for unstable particle;
    @property spin: spin in type 'Fraction';
    @property parity: integer, 1 for positive and -1 for negative.
    """
    def __init__(self, z, a, m, me, exc, t, spin, parity):
        self.z: int = z
        self.a: int = a
        self.m: int = m
        self.mass_excess: PhysicsValue = PhysicsValue(me, 'keV')
        self.half_life: PhysicsValue = PhysicsValue(t, 's')
        self.spin: Fraction = Fraction(spin)
        self.parity: int = parity

    @property
    def name(self) -> str:
        if self.z == 0:
            return 'n'
        else:
            return '%s-%d%s' % (ELEMENTS[self.z], self.a, ISOMERS[self.m])

    @property
    def mass(self) -> PhysicsValue:
        return self.a * u + self.mass_excess * c ** -2

    @property
    def rest_energy(self) -> PhysicsValue:
        return self.a * u * c ** 2 + self.mass_excess

    @property
    def binding_energy(self) -> PhysicsValue:
        return self.z * PROTON_ENERGY + (self.a - self.z) * NEUTRON_ENERGY - self.mass_excess

    def __str__(self):
        return 'Nuclide<%s>' % self.name

    __repr__ = __str__


class DataManager:
    def __init__(self):
        self.conn = connect(DATABASE_PATH)
        self.cursor = self.conn.cursor()

    def __del__(self):
        self.cursor.close()
        self.conn.close()

    def find_nuclide(self, z: int | str, a: int, m: int | str = 0) -> Nuclide:
        """
        search a nuclide by its atomic number 'z', mass number 'a' and isomer id 'm'.

        :param z: atomic number, also accept element symbol;
        :param a: mass number;
        :param m: isomer id, m = 0 means ground state, and m = 1, 2, ... for isomers, default to 0, also accept symbol
                  m, n, p, q, r, i, j;
        :return: the nuclide object, see help(physics_value.nuclide.Nuclide) for more information.
        :raise UnknownElement: if atomic number z or isomer number m is illegal;
        :raise UnknownNuclide: if the nuclide is not in NuBase 2020.
        """
        if isinstance(z, str):
            try:
                z = INVERT_ELEMENTS[z]
            except KeyError:
                raise UnknownElement('Unknown element symbol %s.' % z)
        else:
            if z < 0 or z > 118:
                raise UnknownElement('Illegal atomic number %d.' % z)

        if isinstance(m, str):
            try:
                m = ISO_NUM[m]
            except KeyError:
                raise UnknownElement('Unknown isomers symbol %s.' % m)
        else:
            if m < 0 or 5 < m < 8 or m > 9:
                raise UnknownElement('Illegal isomers number %d.' % m)

        self.cursor.execute('select * from nuclide where z = %d and a = %d and m = %d' % (z, a, m))
        dt = self.cursor.fetchone()
        if dt is None:
            raise UnknownNuclide('Can\'t find nuclide %s-%d%s in NuBase 2020.' % (ELEMENTS[z], a, ISOMERS[m]))
        else:
            _, nz, na, nm, nd, nexc, nt, ns, np = dt
            return Nuclide(nz, na, nm, nd, nexc, nt, ns, np)

    def find_stable_isotope(self, z: int | str) -> list[int]:
        if isinstance(z, str):
            try:
                z = INVERT_ELEMENTS[z]
            except KeyError:
                raise UnknownElement('Unknown element symbol %s.' % z)
        else:
            if z < 0 or z > 118:
                raise UnknownElement('Illegal atomic number %d.' % z)
        self.cursor.execute('select * from nuclide where z = %d and half_life = -1' % z)
        dt = self.cursor.fetchall()
        return [n[2] for n in dt]


if __name__ == '__main__':
    dm = DataManager()
    nuc = dm.find_nuclide('U', 235)
