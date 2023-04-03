import sympy
from sympy.vector import CoordSys3D, Laplacian
import main
from fenics import *

N = CoordSys3D('N')  # Erstellen eines Koordinatensystems

V = sympy.Symbol("V")

# Erstellen von Symbolen und Funktionen
u = sympy.Function('u')(V)  # Trial
test_function = sympy.Function("test_function")(V)
f = sympy.Function('f')(N.x, N.y, N.z)  # input
# Berechnung der Poisson-Gleichung
lap_u = Laplacian(u)  # Laplace-Operator auf u

lap_f = Laplacian(f)

poisson_equation = sympy.Eq(-lap_u, f)


result = main.solve(poisson_equation, test_function)

