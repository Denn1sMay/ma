import sympy
from sympy.vector import CoordSys3D, Laplacian
import main
from fenics import *


print("")
f = sympy.parse_expr("lap(u(x,y,z))", evaluate=True)
sympy.pprint(f)
print(f.args)

N = CoordSys3D('N')  # Erstellen eines Koordinatensystems


V = sympy.Symbol("V")

# Erstellen von Symbolen und Funktionen
u = sympy.Function('u')(V)  # Trial
f = sympy.Function('f')(N.x, N.y, N.z)  # input
test_function = sympy.Function("test_function")(V)
# Berechnung der Poisson-Gleichung
lap_u = Laplacian(u)  # Laplace-Operator auf u

lap_f = Laplacian(f)

poisson_equation = sympy.Eq((lap_u*2)/3 + 3, f)


result = main.solve(poisson_equation, test_function)

string_res = str(result.lhs)
st = "2*inner_func(nabla_grad(u(V)), nabla_grad(test_function(V)))/3 + Integral(3*test_function(V))"
print(st)

fe = eval(st)
print(fe)

