from fractions import Fraction
from math import e as exp_base, pi

if __name__ == '__main__' or __name__ == 'constant':
    from physics_value import PhysicsValue
else:
    from .physics_value import PhysicsValue


one2nd = Fraction(1, 2)
one3rd = Fraction(1, 3)
one4th = Fraction(1, 4)
one5th = Fraction(1, 5)
two3rd = Fraction(2, 3)
two5th = Fraction(2, 5)
three4th = Fraction(3, 4)
three5th = Fraction(3, 5)
four5th = Fraction(4, 5)

exp_base = PhysicsValue(exp_base)
pi = PhysicsValue(pi)

c = PhysicsValue(2.99792458, 8, [0, 0, 0, 0, 1, 0, -1])
h = PhysicsValue(6.62607015, -34, [0, 0, 0, 1, 2, 0, -1])
hbar = PhysicsValue(1.054571817, -34, [0, 0, 0, 1, 2, 0, -1])
e = PhysicsValue(1.602176634, -19, [1, 0, 0, 0, 0, 0, 1])
kB = PhysicsValue(1.3806505, -23, [0, 0, -1, 1, 2, 0, -2])
Vm = PhysicsValue(2.2413996, -2, [0, 0, 0, 0, 3, -1, 0])
NA = PhysicsValue(6.02214076, 23, [0, 0, 0, 0, 0, -1, 0])
R = PhysicsValue(8.314462618, 0, [0, 0, -1, 1, 2, -1, -2])
sigma = PhysicsValue(5.670374419, -8, [0, 0, -4, 1, 0, 0, -3])
G = PhysicsValue(6.67430, -11, [0, 0, 0, -1, 3, 0, -2])
mu0 = PhysicsValue(1.25663706212, -6, [-2, 0, 0, 1, 1, 0, -2])
epsl0 = PhysicsValue(8.8541878128, -12, [2, 0, 0, -1, -3, 0, 4])
mp = PhysicsValue(1.672621777, -27, [0, 0, 0, 1, 0, 0, 0])
mn = PhysicsValue(1.674927351, -27, [0, 0, 0, 1, 0, 0, 0])
u = PhysicsValue(1.66053886, -27, [0, 0, 0, 1, 0, 0, 0])
me = PhysicsValue(9.10938291, -31, [0, 0, 0, 1, 0, 0, 0])
