### 概览

physics_value（以下简称phv）自定义包共包括三个子模块，physics_value实现了自定义类PhysicsValue，提供物理量（数值+单位）的存储与计算，同时可用作量纲分析；constant提供若干物理学常数的数值；nuclide提供核素数据库NuBase2020以供调用。

### 环境配置

Python 3.12+，无其他要求。

### 用法示例

##### 1. 声明一个含有单位的物理量

若单位为单个SI单位或其导出单位，可以直接使用PhysicsValue(数值, 单位)的格式输入；若单位为多个SI单位的复合单位，可以使用关键字参数进行设置：

```python
from physics_value import PhysicsValue

Q = PhysicsValue(1.4, 'kJ')
g = PhysicsValue(9.802, m=1, s=-2)
```

每个单位的前面都可以添加一个工程符号，程序收录了从q（1e-30）到Q（1e30）共23个工程符号，以及基础物理学、核物理中常见的40个单位。

若数值已经用科学计数法表示，则可以分别设置系数和指数：

```python
ME = PhysicsValue(5.965, 24, 'kg')
```

若该物理量是从其他程序输出中读取的或网页上复制的，且单位为单个SI单位或其导出单位，可以不区分数值与单位，作为单个字符串输入：

```python
Q = PhysicsValue('1.4kJ')
```

##### 2. 基础数学运算

phv中的物理量可直接进行加、减、乘、除、乘方五种基本运算与比较运算，进行加、减、比较时程序会检查量纲是否相同（与0相加减、比较时不检查量纲），若不同则拒绝计算并报错，乘除运算无要求，乘方运算要求幂次为整数（`int`）或分数（`fractions.Fraction`），为便于调用，constant中提供了所有分母不超过5的真分数。

```python
x = PhysicsValue(2, 'm')
y = PhysicsValue(3, 'm')
t = PhysicsValue(1, 's')

print(2 * x + 2 * y)
# 1.000 * 10^1 m

print(x * y)
# 6.000 m^2

print(x > y)
# False

print(x + t)
# UnmatchUnitError: ...

print(x > 0)
# True

print(x / t)
# 2.000 m*s^-1

print(x ** 2)
# 4.000 m^2
```

此外，phv提供了开平方（`sqrt`）与开立方（`cbrt`）两种常见的开方运算，其他开方运算需要转换成幂运算进行计算。

```python
from physics_value import sqrt
from physics_value.constant import one4th

a = PhysicsValue(16, m=4)
print(sqrt(a))
# 4.000 m^2

print(a ** one4th)
# 2.000 m
```

##### 3. 高级数学运算

phv提供了幂函数（`exp`）、自然对数（`ln`）与以10为底的对数（`lg`）、六种三角函数（`sin`、`cos`、`tan`、`cot`、`sec`、`csc`），进行这些运算时phv会检查该物理量是否具有常数量纲，若有则输出正常结果，否则拒绝运算并报错。特别地，由于phv使用科学计数法存储数值，允许计算结果超过1e308（无上限）

```python
from physics_value import exp

a = PhysicsValue(15)
print(exp(a))
# 3.269 * 10^6

x = PhysicsValue(2, 'm')
print(exp(x))
# ExponentWithUnitError: ...

print(exp(114514))
# 6.285 * 10^49732
```

phv还提供了类型转换函数`float`，可以将物理量转换为在国际单位制下的数值。对于其他高级数学运算需自行检验是否为常数量纲，再使用标准库`math`中的函数进行运算。

##### 4. 输出

物理量的默认输出格式（`str` `repr`）使用具有4为有效数字的科学计数法表示数值，若量纲对应某个已有的SI导出单位则直接使用该单位，否则使用SI基础单位的幂次表示（不会尝试表示为某两个SI导出单位的积）。类方法`fmt`允许设置输出的有效数字数量（参数`dgt`），是否使用科学计数法（参数`usn`）和工程符号（参数`ues`），是否进行量纲检查（参数`check`），并可以指定输出的单位。类方法`in_unit`可以输出物理量在某个指定单位下的数字，可传入其他程序中进行计算。

```python
length = PhysicsValue('148m')
area = length ** 2
print(area)
# 2.190 * 10^4 m^2

print(area.fmt(dgt=6, ues=False, usn=False))
# 21904.0m^2

print(area.fmt(dgt=3, ues=True, usn=False))
# 21.9k m^2

print(area.fmt(dgt=3, ues='M', usn=False))
# 0.0219M m^2

print(area.fmt(km=2, dgt=3, ues=False, usn=True))
# 2.19 * 10^-2 km^2

print(area.fmt('s', check=True))
# UnmatchUnitError: ...

# NOTICE: if scientific notation is necessary to show precision, kwarg 'usn=False' will never work.
print(area.fmt(dgt=4, ues=False, usn=False))
# 2.190 * 10^4 m^2

print(area.in_unit(m=2), type(area.in_unit(m=2)))
# 21904.0, <class 'float'>
```

##### 5. 物理学常数调用

constant中存储了16个物理学常数，直接调用即可。

nuclide存储了NuBase收录的5843中核素的质量、半衰期、自旋宇称，对于激发态还记录了激发态能量。若有需要可以调用管理器`DataManager`的实例，利用方法`find_nuclide`根据原子序数（或元素符号）、质量数、激发态编号（或符号）查找核素

```python
from physics_value.constant import kB
from physics_value.nuclide import DataManager

print(kB)
# 1.381 * 10^-23 K^-1*kg*m^2*s^-2 (Boltzmann常数)

dm = DataManager()
Cs137 = dm.find_nuclide('Cs', 137)
print(Cs137, Cs137.mass)
# Nuclide<Cs-137> 2.273 * 10^-25 kg

Ag110m = dm.find_nuclide('Ag', 110, 'm')
print(Ag110m.half_life)
# 6.600 * 10^-7 s
```

