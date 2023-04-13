import sympy
from sympy.vector import Laplacian
from dolfinx import *

V = sympy.Symbol("V")

u = sympy.Function('u')(V)  # Trial
test_function = sympy.Function("test_function")(V)
# Berechnung der Poisson-Gleichung
lap_u = Laplacian(u)  # Laplace-Operator auf u


omega = sympy.Symbol("omega")
term = sympy.Integral(3 + sympy.Pow(lap_u, 3) + 5, omega)

lap = term.atoms(sympy.Pow)

sympy.pprint(lap)
print(len(lap))
#sympy.pprint(term)

