"""
A library for physical values and their units. Here PhysicsValue object will record and calculate numbers with unit,
and physical constants can be imported too.

Author: Xmijybke;
Python version: 3.12 or newer publication;
"""


from fractions import Fraction
import math
from typing import Any, Self, overload


class UnitIdentifyError(ValueError):
    """
    raise when program find a strange unit
    """


class UnmatchUnitError(ValueError):
    """
    raise when two value in different unit are summing or comparing
    """


class ExponentWithUnitError(ValueError):
    """
    raise when a PhysicsValue with not 1 unit participate in power and logarithm
    """


class InsignificantResultError(FloatingPointError):
    """
    raise when a PhysicsValue more than 10^15 participate in trigonometric function
    """


ENGINEER_SYMBOL = {'q': -30, 'r': -27, 'y': -24, 'z': -21, 'a': -18, 'f': -15, 'p': -12, 'n': -9, 'u': -6, 'μ': -6,
                   'm': -3, 'c': -2, 'd': -1, '': 0, 'k': 3, 'M': 6, 'G': 9, 'T': 12, 'P': 15, 'E': 18, 'Z': 21,
                   'Y': 24, 'R': 27, 'Q': 30}
INVERT_ENGINEER_SYMBOL = {-30: 'q', -27: 'r', -24: 'y', -21: 'z', -18: 'a', -15: 'f', -12: 'p', -9: 'n', -6: 'μ',
                          -3: 'm', 0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T', 15: 'P', 18: 'E', 21: 'Z', 24: 'Y',
                          27: 'R', 30: 'Q'}
SI_BASIC_UNITS = ['A', 'cd', 'K', 'kg', 'm', 'mol', 's']
CONSTANT = [0, 0, 0, 0, 0, 0, 0]

SI_EXTEND_UNITS = {(0, 0, 0, 1, 0, 0, 0): 'kg',
                   (0, 0, 0, 0, 1, 0, 0): 'm',
                   (0, 0, 0, 0, 0, 0, 1): 's',
                   (1, 0, 0, 0, 0, 0, 0): 'A',
                   (0, 0, 1, 0, 0, 0, 0): 'K',
                   (0, 0, 0, 0, 0, 1, 0): 'mol',
                   (0, 1, 0, 0, 0, 0, 0): 'cd',
                   (0, 0, 0, 1, 1, 0, -2): 'N',
                   (0, 0, 0, 1, -1, 0, -2): 'Pa',
                   (0, 0, 0, 1, 2, 0, -2): 'J',
                   (0, 0, 0, 1, 2, 0, -3): 'W',
                   (1, 0, 0, 0, 0, 0, 1): 'C',
                   (-1, 0, 0, 1, 2, 0, -3): 'V',
                   (-2, 0, 0, -1, -2, 0, 4): 'F',
                   (-2, 0, 0, 1, 2, 0, -3): 'Ω',
                   (2, 0, 0, -1, -2, 0, 3): 'S',
                   (-2, 0, 0, 1, 2, 0, -2): 'H',
                   (-1, 0, 0, 1, 0, 0, -2): 'T',
                   (-1, 0, 0, 1, 2, 0, -2): 'Wb',
                   (0, 0, 0, 0, 0, 0, -1): 'Hz',
                   (0, 0, 0, 0, 2, 0, -2): 'Gy'}

UNITS = {'g': (1, -3, [0, 0, 0, 1, 0, 0, 0]),
         'm': (1, 0, [0, 0, 0, 0, 1, 0, 0]),
         's': (1, 0, [0, 0, 0, 0, 0, 0, 1]),
         'A': (1, 0, [1, 0, 0, 0, 0, 0, 0]),
         'K': (1, 0, [0, 0, 1, 0, 0, 0, 0]),
         'mol': (1, 0, [0, 0, 0, 0, 0, 1, 0]),
         'cd': (1, 0, [0, 1, 0, 0, 0, 0, 0]),
         '': (1, 0, [0, 0, 0, 0, 0, 0, 0]),
         'b': (1, -28, [0, 0, 0, 0, 2, 0, 0]),
         'min': (6, 1, [0, 0, 0, 0, 0, 0, 1]),
         'h': (3.6, 3, [0, 0, 0, 0, 0, 0, 1]),
         'd': (8.64, 4, [0, 0, 0, 0, 0, 0, 1]),
         'a': (3.1558150, 7, [0, 0, 0, 0, 0, 0, 1]),
         'yr': (3.1558150, 7, [0, 0, 0, 0, 0, 0, 1]),
         'N': (1, 0, [0, 0, 0, 1, 1, 0, -2]),
         'J': (1, 0, [0, 0, 0, 1, 2, 0, -2]),
         'W': (1, 0, [0, 0, 0, 1, 2, 0, -3]),
         'Pa': (1, 0, [0, 0, 0, 1, -1, 0, -2]),
         'Bar': (1, 5, [0, 0, 0, 1, -1, 0, -2]),
         'atm': (1.01325, 5, [0, 0, 0, 1, -1, 0, -2]),
         'mHg': (1.3332, -1, [0, 0, 0, 1, -1, 0, -2]),
         'C': (1, 0, [1, 0, 0, 0, 0, 0, 1]),
         'V': (1, 0, [-1, 0, 0, 1, 2, 0, -3]),
         'F': (1, 0, [-2, 0, 0, -1, -2, 0, 4]),
         'Ohm': (1, 0, [-2, 0, 0, 1, 2, 0, -3]),
         'Ω': (1, 0, [-2, 0, 0, 1, 2, 0, -3]),
         'S': (1, 0, [2, 0, 0, -1, -2, 0, 3]),
         'H': (1, 0, [-2, 0, 0, 1, 2, 0, -2]),
         'T': (1, 0, [-1, 0, 0, 1, 0, 0, -2]),
         'Wb': (1, 0, [-1, 0, 0, 1, 2, 0, -2]),
         'eV': (1.602176634, -19, [0, 0, 0, 1, 2, 0, -2]),
         'Hz': (1, 0, [0, 0, 0, 0, 0, 0, -1]),
         'Bq': (1, 0, [0, 0, 0, 0, 0, 0, -1]),
         'Ci': (3.7, 10, [0, 0, 0, 0, 0, 0, -1]),
         'Gy': (1, 0, [0, 0, 0, 0, 2, 0, -2]),
         'Sv': (1, 0, [0, 0, 0, 0, 2, 0, -2]),
         'rem': (1, -2, [0, 0, 0, 0, 2, 0, -2]),
         'R': (2.58, -4, [1, 0, 0, -1, 0, 0, 1]),
         'rad': (1, -2, [0, 0, 0, 0, 2, 0, -2]),
         'cal': (4.1859, 0, [0, 0, 0, 1, 2, 0, -2]),
         'Gs': (1, -4, [-1, 0, 0, 1, 0, 0, -2])}


def scientific_notation(_x: float) -> tuple[float, int]:
    """
    convert a number(n) to scientific_notation(a, p)\n
    let x = a * 10 ^ p (1 <= |a| < 10, p ∈ Z)\n
    >>> scientific_notation(2.5)
    (2.5, 0)
    >>> scientific_notation(0.25)
    (2.5, -1)
    >>> scientific_notation(25)
    (2.5, 1)

    :return: a, p
    """
    if _x >= 1:
        _p = int(math.log10(_x))
        _a = _x / 10 ** _p
        return _a, _p
    elif _x > 0:
        _p = math.ceil(-math.log10(_x))
        _a = _x * 10 ** _p
        return _a, -_p
    elif _x == 0:
        return 0, 0
    elif _x > -1:
        _p = math.ceil(-math.log10(abs(_x)))
        _a = _x * 10 ** _p
        return _a, -_p
    else:
        _p = int(math.log10(abs(_x)))
        _a = _x / 10 ** _p
        return _a, _p


def _unit_power(unit: str, power: int | Fraction = 1) -> tuple[float, int, list[int | Fraction]]:
    """
    raise a unit to nth power
    """
    if unit in UNITS:
        _a, _p, _u = UNITS[unit]
        return _a ** power, _p * power, [i * power for i in _u]
    else:
        es, unit = unit[0], unit[1:]
        try:
            _e = ENGINEER_SYMBOL[es]
            _a, _p, _u = UNITS[unit]
        except KeyError:
            raise UnitIdentifyError('Unknown unit %s%s.' % (es, unit))
        else:
            return _a ** power, (_p + _e) * power, [i * power for i in _u]


class PhysicsValue(object):
    """
    this class storage physical values and their units

    Excepted __init__() arguments
    -----------------------------
    PhysicsValue(x: str) return x;
    PhysicsValue(v: float, **units: int | Fraction) return v in units;
    PhysicsValue(v: float, unit: str) return v in unit;
    PhysicsValue(a: float, p: int, **units: int | Fraction) return a * 10^p in units;
    PhysicsValue(a: float, p: int, unit: str) return a * 10^p in unit.

    Examples
    --------
    >>> PhysicsValue('1.4kJ')
    1.400 * 10^3 J

    >>> PhysicsValue(487, 'cal')
    2.039 * 10^3 J

    """
    __slots__ = ['a', 'p', 'u']

    @overload
    def __init__(self, value: str): ...

    @overload
    def __init__(self, value: int | float, **kwargs: int): ...

    @overload
    def __init__(self, value: int | float, unit: str): ...

    @overload
    def __init__(self, _a: int | float, _p: int, **kwargs: int): ...

    @overload
    def __init__(self, _a: int | float, _p: int, unit: str): ...

    @overload
    def __init__(self, _a: int | float, _p: int, _u: list[int | Fraction]): ...

    def __init__(self, *args: Any, **kwargs: int):
        self.a: float = 0
        self.p: int = 0
        self.u: list[int | Fraction] = [0, 0, 0, 0, 0, 0, 0]

        if len(args) == 1 and isinstance(args[0], str):
            value = args[0]
            for i in range(len(value) - 1, -1, -1):
                if value[i].isdigit():
                    _n, _u = value[:i + 1], value[i + 1:]
                    try:
                        self.a, self.p = scientific_notation(float(_n))
                    except ValueError:
                        continue
                    else:
                        self.__set_unit(_u)
                    break
            else:
                self.a, self.p = 1, 0
                self.__set_unit(value)

        elif len(args) == 1 and isinstance(args[0], (int, float)):
            self.a, self.p = scientific_notation(args[0])
            for unit in kwargs:
                self.__set_unit(unit, kwargs[unit])

        elif len(args) == 2 and isinstance(args[0], (int, float)) and isinstance(args[1], str):
            self.a, self.p = scientific_notation(args[0])
            self.__set_unit(args[1])

        elif len(args) == 2 and isinstance(args[0], (int, float)) and isinstance(args[1], int):
            _a, _p = args
            if 1 <= abs(_a) < 10:
                self.a, self.p = _a, _p
            else:
                _a, dp = scientific_notation(_a)
                self.a, self.p = _a, _p + dp
            for unit in kwargs:
                self.__set_unit(unit, kwargs[unit])

        elif len(args) == 3 and isinstance(args[0], (int, float)) and isinstance(args[1], int) and isinstance(args[2], str):
            _a, _p, _u = args
            if 1 <= abs(_a) < 10:
                self.a, self.p = _a, _p
            else:
                _a, dp = scientific_notation(_a)
                self.a, self.p = _a, _p + dp
            self.__set_unit(_u)

        elif len(args) == 3 and isinstance(args[0], (int, float)) and isinstance(args[1], int) and isinstance(args[2], list):
            _a, _p, _u = args
            if len(_u) == 7:
                if 1 <= abs(_a) < 10:
                    self.a, self.p, self.u = _a, _p, _u
                else:
                    _a, dp = scientific_notation(_a)
                    self.a, self.p, self.u = _a, _p + dp, _u
            else:
                raise UnitIdentifyError('Meow!')

    def __set_unit(self, unit: str, power: int = 1) -> None:
        if unit not in UNITS:
            eng_smb = unit[0]
            if eng_smb in ENGINEER_SYMBOL:
                self.p += power * ENGINEER_SYMBOL[eng_smb]
            else:
                raise UnitIdentifyError(f'Can\'t indentify unit {unit}')
            unit = unit[1:]
        if unit in UNITS:
            unit_a, unit_p, unit_u = UNITS[unit]
            self.a *= unit_a ** power
            self.p += unit_p * power
            self.u = [unit_u[i] * power + self.u[i] for i in range(7)]
            if abs(self.a) < 1 or abs(self.a) > 10:
                self.a, dp = scientific_notation(self.a)
                self.p += dp
        else:
            raise UnitIdentifyError(f'Can\'t indentify unit {unit}')

    @property
    def unit_str(self) -> str:
        """
        return the value's unit in str, default to SI basic and extend units.

        Example
        -------
        >>> value = PhysicsValue(1.3, J=1, s=-1)
        >>> value.unit_str
        'W'

        :return: the value's unit.
        """
        unit = tuple(self.u)
        if unit in SI_EXTEND_UNITS:
            return SI_EXTEND_UNITS[unit]
        else:
            unit_list = []
            for i in range(7):
                up = self.u[i]
                if up:
                    if isinstance(up, Fraction) and up.denominator != 1:
                        unit_list.append('%s^%s' % (SI_BASIC_UNITS[i], str(up)))
                    elif up != 1:
                        unit_list.append('%s^%d' % (SI_BASIC_UNITS[i], up))
                    else:
                        unit_list.append(SI_BASIC_UNITS[i])
            if unit_list:
                return '*'.join(unit_list)
            else:
                return '1'

    def fmt(self, unit: str = None, *, check: bool = False, dgt: int = 4, ues: str | bool = True,
            usn: bool = False, **units: int | Fraction) -> str:
        """
        convert this value to string in specific format

        :param check: check this value in 'unit' or not. If check but value not in the unit, raise UnmatchUnitError,
                      else convert it anyway;

        :param dgt: number of significant digits, default to 4;

        :param ues: Using Engineer Symbol, accept y(1e-24), z(1e-21), a(1e-18), ..., Y(1e24) and boolean.
                    If set as 'True', a suitable symbol will be used to make 1 <= n < 1000;

        :param usn: Using Scientific Notation, default to True, if disable and the number is no more than
                    1e308, the result number will be output directly;

        :param unit: assign the unit as a sigle SI or extend unit.

        :param units: assign a complex unit by keyword argument.

        Example
        -------
        >>> length = PhysicsValue('148m')
        >>> area = length ** 2
        >>> print(area)
        2.190 * 10^4 m^2

        >>> print(area.fmt(dgt=6, ues=False, usn=False))
        21904.0m^2

        >>> print(area.fmt(dgt=3, ues=True, usn=False))
        21.9k m^2

        >>> print(area.fmt(dgt=3, ues='M', usn=False))
        0.0219M m^2

        >>> print(area.fmt(km=2, dgt=3, ues=False, usn=True))
        2.19 * 10^-2 km^2

        >>> print(area.fmt('s', check=True))
        Traceback (most recent call last):
        UnmatchUnitError: This value has unit in m^2 but output as s.

        >>> # NOTICE: if scientific notation is necessary to show precision, kwarg 'usn=False' will never work.
        >>> print(area.fmt(dgt=4, ues=False, usn=False))
        2.190 * 10^4 m^2

        (Not 21900, it shows the wrong precision)
        """
        if not isinstance(dgt, int) or dgt <= 0:
            raise ValueError('Argument dgt (significant digit) expected positive integer, got %s.' % str(dgt))

        _a, _p, _u, output_units = self.a, self.p, [i for i in self.u], []
        if unit is not None:
            output_units.append(unit)
            a2, p2, u2 = _unit_power(unit, 1)
            _a, _p, _u = _a / a2, _p - p2, [_u[i] - u2[i] for i in range(7)]
        else:
            for vu, vp in units.items():
                output_units.append('%s^%s' % (vu, str(vp)))
                a2, p2, u2 = _unit_power(vu, vp)
                _a, _p, _u = _a / a2, _p - p2, [_u[i] - u2[i] for i in range(7)]


        if abs(_a) < 1:
            _a, _p = _a * 10, _p - 1

        if _u != CONSTANT:
            if check:
                raise UnmatchUnitError(f'This value has unit in {self.unit_str} but output as {unit}.')
            else:
                if tuple(_u) in SI_EXTEND_UNITS:
                    output_units.append(SI_EXTEND_UNITS[tuple(_u)])
                else:
                    for i in range(7):
                        if _u[i] == 1:
                            output_units.append(SI_BASIC_UNITS[i])
                        elif _u[i] != 0:
                            output_units.append('%s^%s' % (SI_BASIC_UNITS[i], str(_u[i])))
        output_string = '*'.join(output_units)

        if ues:
            if tuple(self.u) not in SI_EXTEND_UNITS:
                output_string = ' ' + output_string
            if ues is True:
                if _p <= -24:
                    ues = 'y'
                elif _p <= 23:
                    ues = INVERT_ENGINEER_SYMBOL[_p - _p % 3]
                else:
                    ues = 'Y'
            output_string = ues + output_string
            if ues in ENGINEER_SYMBOL:
                _p -= ENGINEER_SYMBOL[ues]
            else:
                raise UnitIdentifyError('Unknown engineer symbol %s.' % ues)

        if dgt > _p and _p < 308 and usn is False:
            _a *= 10 ** _p
        else:
            output_string = ' * 10^%d ' % _p + output_string
            _p = 0

        # noinspection PyStringFormat
        output_string = f'%.{dgt - _p - 1}f' % _a + output_string
        return output_string

    def in_unit(self, unit: str = None, **kwargs: int) -> float:
        if unit is None:
            a, p, u = 1, 0, CONSTANT[:]
            output_name = []
            for unit, power in kwargs.items():
                a0, p0, u0 = _unit_power(unit, power)
                a *= a0
                p += p0
                u = [u[i] + u0[i] for i in range(7)]
                output_name.append('%s^%d' % (unit, power))
            if self.u == u:
                return self.a / a * 10 ** (self.p - p)
            else:
                raise UnmatchUnitError(f'This value has unit in {self.unit_str} but output as {'*'.join(output_name)}.')

        else:
            a, p, u = _unit_power(unit, 1)
            if self.u == u:
                return self.a / a * 10 ** (self.p - p)
            else:
                raise UnmatchUnitError(f'This value has unit in {self.unit_str} but output as {unit}.')

    def __str__(self):
        if self.u == CONSTANT:
            if self.p:
                return '%.3f * 10^%d' % (self.a, self.p)
            else:
                return '%.3f' % self.a
        else:
            if self.p:
                return '%.3f * 10^%d %s' % (self.a, self.p, self.unit_str)
            else:
                return '%.3f %s' % (self.a, self.unit_str)

    def __repr__(self):
        unit_lst = []
        for i in range(7):
            if self.u[i]:
                if isinstance(self.u[i], Fraction) and self.u[i].denominator != 1:
                    unit_lst.append('%s=%s' % (SI_BASIC_UNITS[i], repr(self.u[i])))
                else:
                    unit_lst.append('%s=%d' % (SI_BASIC_UNITS[i], self.u[i]))
        return 'PhysicsValue(%s, %d, %s)' % (str(self.a), self.p, ', '.join(unit_lst))

    __repr__ = __str__

    def __float__(self):
        return float(self.a * 10 ** self.p)

    def __abs__(self):
        return PhysicsValue(abs(self.a), self.p, self.u)

    def __mul__(self, other: int | float | Self) -> Self:
        if isinstance(other, (int, float)):
            return PhysicsValue(self.a * other, self.p, self.u)

        elif isinstance(other, PhysicsValue):
            return PhysicsValue(self.a * other.a, self.p + other.p, [self.u[i] + other.u[i] for i in range(7)])

        else:
            return NotImplemented

    def __rmul__(self, other: int | float) -> Self:
        if isinstance(other, (int, float)):
            return PhysicsValue(other * self.a, self.p, self.u)
        else:
            return NotImplemented

    def __truediv__(self, other: int | float | Self) -> Self:
        if isinstance(other, (int, float)):
            return PhysicsValue(self.a / other, self.p, self.u)

        elif isinstance(other, PhysicsValue):
            return PhysicsValue(self.a / other.a, self.p - other.p, [self.u[i] - other.u[i] for i in range(7)])

        else:
            return NotImplemented

    def __rdiv__(self, other: int | float) -> Self:
        if isinstance(other, (int, float)):
            return PhysicsValue(other / self.a, -self.p, [-self.u[i] for i in range(7)])
        else:
            return NotImplemented

    __rtruediv__ = __rdiv__

    def __pow__(self, power: int | Fraction, modulo=None) -> Self:
        if isinstance(power, Fraction) and power.denominator != 1:
            rest = power - int(power)
            return PhysicsValue(self.a ** power * 10 ** (self.p * rest), self.p * int(power),
                                [self.u[i] * power for i in range(7)])

        else:
            return PhysicsValue(self.a ** power, self.p * power, [self.u[i] * power for i in range(7)])

    def __add__(self, other: int | float | Self) -> Self:
        if isinstance(other, (int, float)):
            if other == 0:
                return self
            elif self.u == CONSTANT:
                _a, _p = scientific_notation(other)
                if _p <= self.p:
                    return PhysicsValue(self.a + other * 10 ** -self.p, self.p, self.u)
                else:
                    return PhysicsValue(self.a * 10 ** (self.p - _p) + _a, _p, self.u)
            else:
                raise UnmatchUnitError(f'Trying to add two PhysicsValue with unit {self.unit_str} and 1.')

        elif isinstance(other, PhysicsValue):
            if self.u == other.u:
                if self.p < other.p:
                    return PhysicsValue(self.a * 10 ** (self.p - other.p) + other.a, other.p, self.u)
                else:
                    return PhysicsValue(self.a + other.a * 10 ** (other.p - self.p), self.p, self.u)
            else:
                raise UnmatchUnitError(f'Trying to add two PhysicsValue with unit {self.unit_str} and {other.unit_str}.')

        else:
            return NotImplemented

    def __radd__(self, other: int | float) -> Self:
        if isinstance(other, (int, float)):
            if other == 0:
                return self
            elif self.u == CONSTANT:
                _a, _p = scientific_notation(other)
                if _p <= self.p:
                    return PhysicsValue(self.a + other * 10 ** -self.p, self.p, self.u)
                else:
                    return PhysicsValue(self.a * 10 ** (self.p - _p) + _a, _p, self.u)
            else:
                raise UnmatchUnitError(f'Trying to add two PhysicsValue with unit 1 and {self.unit_str}.')

        else:
            return NotImplemented

    def __sub__(self, other: int | float | Self) -> Self:
        if isinstance(other, (int, float)):
            if other == 0:
                return self
            elif self.u == CONSTANT:
                _a, _p = scientific_notation(other)
                if _p <= self.p:
                    return PhysicsValue(self.a - other * 10 ** -self.p, self.p, self.u)
                else:
                    return PhysicsValue(self.a * 10 ** (self.p - _p) - _a, _p, self.u)
            else:
                raise UnmatchUnitError(f'Trying to subtract two PhysicsValue with unit {self.unit_str} and 1.')

        elif isinstance(other, PhysicsValue):
            if self.u == other.u:
                if self.p < other.p:
                    return PhysicsValue(self.a * 10 ** (self.p - other.p) - other.a, other.p, self.u)
                else:
                    return PhysicsValue(self.a - other.a * 10 ** (other.p - self.p), self.p, self.u)
            else:
                raise UnmatchUnitError(
                    f'Trying to subtract two PhysicsValue with unit {self.unit_str} and {other.unit_str}.')

        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return self
            elif self.u == CONSTANT:
                _a, _p = scientific_notation(other)
                if _p <= self.p:
                    return PhysicsValue(other * 10 ** -self.p - self.a, self.p, self.u)
                else:
                    return PhysicsValue(_a - self.a * 10 ** (self.p - _p), _p, self.u)
            else:
                raise UnmatchUnitError(f'Trying to subtract two PhysicsValue with unit 1 and {self.unit_str}.')

        else:
            return NotImplemented

    def __neg__(self):
        return PhysicsValue(-self.a, self.p, self.u)

    def __eq__(self, other: int | float | Self) -> bool:
        if isinstance(other, (int, float)):
            if other == 0:
                return self.a == 0
            elif self.u == CONSTANT:
                return self.a * 10 ** self.p == other
            else:
                return False

        elif isinstance(other, PhysicsValue):
            return (self.u == other.u) and (self.p == other.p) and (self.a == other.a)

        else:
            return NotImplemented

    def __lt__(self, other: int | float | Self) -> bool:
        if isinstance(other, (int, float)):
            if other == 0:
                return self.a < 0
            elif self.u == CONSTANT:
                return self.a * 10 ** self.p < other
            else:
                raise UnmatchUnitError(f'Trying to compare two PhysicsValue with unit {self.unit_str} and 1.')

        elif isinstance(other, PhysicsValue):
            if self.u == other.u:
                if self.p < other.p:
                    return True
                elif self.p > other.p:
                    return False
                else:
                    return self.a < other.a
            else:
                raise UnmatchUnitError(
                    f'Trying to compare two PhysicsValue with unit {self.unit_str} and {other.unit_str}.')

        else:
            return NotImplemented

    def __gt__(self, other) -> bool:
        if isinstance(other, (int, float)):
            if other == 0:
                return self.a > 0
            elif self.u == CONSTANT:
                return self.a * 10 ** self.p > other
            else:
                raise UnmatchUnitError(f'Trying to compare two PhysicsValue with unit {self.unit_str} and 1.')

        elif isinstance(other, PhysicsValue):
            if self.u == other.u:
                if self.p < other.p:
                    return False
                elif self.p > other.p:
                    return True
                else:
                    return self.a > other.a
            else:
                raise UnmatchUnitError(
                    f'Trying to compare two PhysicsValue with unit {self.unit_str} and {other.unit_str}.')

        else:
            return NotImplemented


def lg(value: int | float | PhysicsValue,
       unit: PhysicsValue | str | tuple[str, int] | list[str | tuple[str, int]]) -> PhysicsValue:
    """
    lg(x)
    """

    if isinstance(value, (int, float)):
        return PhysicsValue(math.log10(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            return PhysicsValue(value.p + math.log10(value.a))
        else:
            raise ExponentWithUnitError(f'PhysicsValue with unit {value.unit_str} participate in logarithm.')

    else:
        return NotImplemented


def ln(value: int | float | PhysicsValue):
    """
    ln(x) = log_e(x)
    """
    if isinstance(value, (int, float)):
        return math.log(value)

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            return value.p * math.log(10) + math.log(value.a)
        else:
            raise ExponentWithUnitError(f'PhysicsValue with unit {value.unit_str} participate in logarithm.')

    else:
        return NotImplemented


def exp(value: int | float | PhysicsValue) -> PhysicsValue:
    """
    exp(x) = e^x
    """
    if isinstance(value, (int, float)):
        multi = 1
        while True:
            try:
                result = PhysicsValue(math.exp(float(value) / multi))
            except OverflowError:
                multi *= 2
            else:
                return result ** multi

    elif isinstance(value, PhysicsValue):
        if value.a == 0:
            return 1
        if value.u == CONSTANT:
            multi = 1
            while True:
                try:
                    result = PhysicsValue(math.exp(float(value) / multi))
                except OverflowError:
                    multi *= 2
                else:
                    return result ** multi
        else:
            raise ExponentWithUnitError(f'PhysicsValue with unit {value.unit_str} act as exponent.')

    else:
        return NotImplemented


def sqrt(value: int | float | PhysicsValue) -> PhysicsValue:
    """
    sqrt(x) = x ** 0.5
    """
    if isinstance(value, (int, float)):
        return PhysicsValue(math.sqrt(value))

    elif isinstance(value, PhysicsValue):
        if value.p % 2:
            new_a = value.a ** 0.5 * 10 ** 0.5
        else:
            new_a = value.a ** 0.5
        new_lst = [Fraction(value.u[i], 2) if value.u[i] % 2 else value.u[i] // 2 for i in range(7)]
        return PhysicsValue(new_a, value.p // 2, new_lst)


def cbrt(value: int | float | PhysicsValue) -> PhysicsValue:
    """
    cbrt(x) = x ** (1 / 3)
    """
    if isinstance(value, (int, float)):
        return PhysicsValue(math.cbrt(value))

    elif isinstance(value, PhysicsValue):
        new_lst = [Fraction(value.u[i], 3) if value.u[i] % 3 else value.u[i] // 3 for i in range(7)]
        return PhysicsValue(math.cbrt(value.a * 10 ** (value.p % 3)), value.p // 3, new_lst)


def sin(value: int | float | PhysicsValue) -> PhysicsValue:
    """sin(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(math.sin(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return value
            else:
                return PhysicsValue(math.sin(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


def cos(value: int | float | PhysicsValue) -> PhysicsValue:
    """cos(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(math.cos(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return PhysicsValue(1)
            else:
                return PhysicsValue(math.cos(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


def tan(value: int | float | PhysicsValue) -> PhysicsValue:
    """tan(x) = sin(x) / cos(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(math.tan(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return value
            else:
                return PhysicsValue(math.tan(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


def cot(value: int | float | PhysicsValue) -> PhysicsValue:
    """cot(x) = cos(x) / sin(x) = 1 / tan(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(1 / math.tan(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return 1 / value
            else:
                return PhysicsValue(1 / math.tan(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


def sec(value: int | float | PhysicsValue) -> PhysicsValue:
    """sec(x) = 1 / cos(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(1 / math.cos(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return PhysicsValue(1)
            else:
                return PhysicsValue(1 / math.cos(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


def csc(value: int | float | PhysicsValue) -> PhysicsValue:
    """csc(x) = 1 / sin(x)"""
    if isinstance(value, (int, float)):
        return PhysicsValue(1 / math.sin(value))

    elif isinstance(value, PhysicsValue):
        if value.u == CONSTANT:
            if value.p >= 15:
                raise InsignificantResultError(
                    'Due to floating-point number\'s limited precision, the result is meaningless.')
            elif value.p <= -308:
                return 1 / value
            else:
                return PhysicsValue(1 / math.sin(float(value)))
        else:
            raise ExponentWithUnitError(
                f'PhysicsValue with unit {value.unit_str} participate in trigonometric function.')


if __name__ == '__main__':
    from doctest import testmod
    testmod()
