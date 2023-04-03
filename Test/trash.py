'''
import sympy

x = sympy.Symbol("x")
cos = sympy.cos
sympy_input = cos(x) + x ** 2


import fenics
mesh = fenics.UnitIntervalMesh(10)
cos = fenics.cos
x = fenics.SpatialCoordinate(mesh)[0]
ufl_output = eval(str(sympy_input))
print(ufl_output)



'''
import sympy
from sympy.vector import Laplacian
from fenics import *

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)
u = TrialFunction(V)
v = TestFunction(V)


#u = SpatialCoordinate(mesh)
#v = SpatialCoordinate(mesh)

st = "inner(grad(u), grad(v))"

fe_ex = eval(st)

print("")
print(type(fe_ex))
print(fe_ex)


inner_args = fe_ex.ufl_operands
print(inner_args)


