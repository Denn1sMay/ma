import sympy
from sympy.vector import CoordSys3D, Laplacian, Del
import Util.utils as utils


def use_sympy_laplace_operator(sympy_term):
    undefined_laplacian_function = sympy.Function("Laplacian")
    lap_atoms = sympy_term.atoms(undefined_laplacian_function)
    new_term = sympy_term
    for lap_atom in lap_atoms:
        laplacian_args = lap_atom.args
        print(type(laplacian_args[0]))
        lap_operator = Laplacian(laplacian_args[0])
        new_term = new_term.subs(lap_atom, lap_operator)
    return new_term


def parse_string_equation(string_equation: str):
    utils.print_space("Parsing String to sympy Equation")
    print(string_equation)
    equation_sides = string_equation.split("=")
    if len(equation_sides) > 2:
        raise ValueError("Equation contained more than one '='")

    if len(equation_sides) == 2:
        lhs_parsed = sympy.parse_expr(equation_sides[0], evaluate=False)
        rhs_parsed = sympy.parse_expr(equation_sides[1], evaluate=False)
        print("lhs:")
        sympy.pprint(lhs_parsed)
        print("rhs:")
        sympy.pprint(rhs_parsed)
        lhs_with_operators = use_sympy_laplace_operator(lhs_parsed)
        rhs_with_operators = use_sympy_laplace_operator(rhs_parsed)
    else:
        lhs_parsed = sympy.parse_expr(equation_sides[0], evaluate=False)
        print("lhs:")
        sympy.pprint(lhs_parsed)
        rhs_with_operators = sympy.parse_expr("0", evaluate=False)
        lhs_with_operators = use_sympy_laplace_operator(lhs_parsed)

    parsed_equation = sympy.Eq(lhs_with_operators, rhs_with_operators)
    print("Result with sympy operators:")
    sympy.pprint(parsed_equation)
    return parsed_equation


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

def multiply_with_test_function(sympy_equation: sympy.Eq, test_function: sympy.Symbol):
    utils.print_space("Multiply with Test Function")
    m = sympy.Eq(sympy_equation.lhs * test_function, sympy_equation.rhs * test_function)
    print("Multiplied with test function - result:")
    sympy.pprint(m)
    return m


def integrate_summands(functions, omega: sympy.Symbol):
    # Distributivgesetz - ausmultiplizieren
    splitted = sympy.expand(functions)
    splitted_args = sympy.Add.make_args(splitted)

    # Integral über jeden Summand bilden und Integrale addieren
    expr = 0
    for term in splitted_args:
        expr = expr + sympy.Integral(term, omega)

    return expr


def sort_equation(lhs, rhs, trial_function_u: sympy.Symbol):
    lhs_atoms = lhs.atoms(sympy.Integral)
    rhs_atoms = rhs.atoms(sympy.Integral)

    lhs_sorted = 0
    rhs_sorted = 0

    for lhs_atom in lhs_atoms:
        sympy.pprint(lhs_atom)
        if lhs_atom.has(trial_function_u):
            lhs_sorted = lhs_sorted + lhs_atom
        else:
            rhs_sorted = rhs_sorted + lhs_atom

    for rhs_atom in rhs_atoms:
        if rhs_atom.has(trial_function_u):
            lhs_sorted = lhs_sorted + rhs_atom
        else:
            rhs_sorted = rhs_sorted + rhs_atom
            
    return sympy.Eq(lhs_sorted, rhs_sorted)



def integrate_over_domain_and_sort(sympy_equation: sympy.Eq, trial_function_u: sympy.Symbol):
    utils.print_space("Integrate over Domain and Sort")
    omega = sympy.Symbol('omega')

    lhs = integrate_summands(sympy_equation.lhs, omega)
    rhs = integrate_summands(sympy_equation.rhs, omega)

    sorted_equation = sort_equation(lhs, rhs, trial_function_u)
    print("Sorted and Integrated over domain:")
    sympy.pprint(sorted_equation)
    return sorted_equation, omega
