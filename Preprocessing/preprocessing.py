import sympy
from sympy.vector import CoordSys3D, Laplacian, Del
import Util.utils as utils

def find_denendent_variables(sympy_equation: sympy.Eq):
    utils.print_space("Get dependent variables")
    print("Laplactian Term:")
    a = sympy_equation.atoms(Laplacian)
    print(a)
    #TODO equation can contain multiple laplace operators
    print(type(a))
    laplacian_function = list(a)[0]
    laplacian_args = laplacian_function.args
    print(type(laplacian_args))
    first_arg = laplacian_args[0]
    print("First argument: ")
    print(first_arg)

    u_args = first_arg.args
    print("u arguments: ")
    print(u_args)
    print(type(u_args))
    return u_args



def get_test_function_v(sympy_equation: sympy.Eq):
    utils.print_space("Creating Test Function with correct args")
    #lap_term = sympy_equation.lhs.coeff(Laplacian())
    
    u_args = find_denendent_variables(sympy_equation)
    
    test_function_v = sympy.Function("test_function")(u_args)

    return test_function_v, u_args

def multiply_with_test_function(sympy_equation: sympy.Eq, test_function: sympy.Function):
    utils.print_space("Multiply with Test Function")
    #Test Funktion muss nach der Variablen differenzierbar sein
    m = sympy.Eq(sympy_equation.lhs * test_function, sympy_equation.rhs * test_function)
    print("Multiplied with test function - result:")
    sympy.pprint(m)
    return m


def integrate_summands(functions, omega: sympy.Symbol):
    ###### Linke Seite
    # Distributivgesetz - ausmultiplizieren
    splitted = sympy.expand(functions)
    splitted_args = sympy.Add.make_args(splitted)

    # Integral über jeden Summand bilden und als Linke Seite der Gleichung festlegen
    expr = 0
    for term in splitted_args:
        expr = expr + sympy.Integral(term, omega)

    return expr


def integrate_over_domain(sympy_equation: sympy.Eq):
    utils.print_space("Integrate over Domain")
    omega = sympy.Symbol('omega')

    lhs = integrate_summands(sympy_equation.lhs, omega)
    rhs = integrate_summands(sympy_equation.rhs, omega)
    integrals = sympy.Eq(lhs, rhs)

    print("Integrated over domain:")
    sympy.pprint(integrals)
    return integrals, omega
