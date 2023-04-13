import sympy
from sympy.vector import CoordSys3D, Laplacian
import main
from dolfinx import *
from mpi4py import MPI
from ufl import *


N = CoordSys3D('N')  # Erstellen eines Koordinatensystems

func_space = sympy.Symbol("V")

# Erstellen von Symbolen und Funktionen
u = sympy.Function('u')(func_space)  # Trial
v = sympy.Function('v')(func_space)
x = sympy.Symbol("x") #sympy.Function('x')(N.x, N.y, N.z)  # input
# Berechnung der Poisson-Gleichung
lap_u = Laplacian(u)  # Laplace-Operator auf u

poisson_equation = sympy.Eq(-lap_u, x)

# ------------------------
domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.quadrilateral)
#mesh = UnitSquareMesh(8, 8)
V = fem.FunctionSpace(domain, ("CG", 1))
trial = TrialFunction(V)
test = TestFunction(V)


result = main.get_weak_form(poisson_equation, u, v, trial, test)