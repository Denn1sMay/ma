import sympy


y = sympy.Symbol("y")
x = sympy.Symbol("x")
term = (y+2)*x

print(sympy.expand(term))