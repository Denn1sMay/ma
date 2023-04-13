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
#from fenics import *
from dolfinx import *
from ufl import *
from mpi4py import MPI


domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.quadrilateral)
#mesh = UnitSquareMesh(8, 8)
V = fem.FunctionSpace(domain, ("CG", 1))
u = TrialFunction(V)
v = TestFunction(V)


#u = SpatialCoordinate(mesh)
#v = SpatialCoordinate(mesh)

st = "inner(grad(u), grad(v)) * dx"

fe_ex = eval(st)

print("")
print(type(fe_ex))
print("")

print(fe_ex)

pure_ufl = inner(grad(u), grad(v)) * dx
print(pure_ufl)
#inner_args = fe_ex.ufl_operands
#print(inner_args)


