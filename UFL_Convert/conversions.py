from dolfinx import *
from ufl import *
from mpi4py import MPI
import sympy
import Util.utils as utils


def replace_trial_function(sympy_function, trial_function: sympy.Symbol, trial_function_name: str):
    unknown_symbolic_trial_function = sympy.Symbol(trial_function_name)
    substituted_term = sympy_function.subs(trial_function, unknown_symbolic_trial_function)
    return substituted_term

def replace_test_function(sympy_function, test_function: sympy.Symbol, test_function_name: str):
    unknown_symbolic_test_function = sympy.Symbol(test_function_name)
    substituted_term = sympy_function.subs(test_function, unknown_symbolic_test_function)
    return substituted_term

def convert_integral_to_ufl_string(integral: sympy.Integral, trial_function: sympy.Symbol, test_function: sympy.Symbol, trial_function_name: str, test_function_name: str):
    print("Converting integral ")
    integral_args = integral.args[0]
    integral_domain = integral.args[1][0]
    replaced_trial = replace_trial_function(integral_args, trial_function, trial_function_name)
    replaced_trial_and_test = replace_test_function(replaced_trial, test_function, test_function_name)
    print("ok")
    replaced_as_string = str(replaced_trial_and_test)
    print(type(integral_domain))
    print(integral_args)
    if integral_domain == sympy.Symbol("omega"):
        print("Omega found")
        string_as_integral_string = "(" + replaced_as_string + ") * dx"
    elif integral_domain == sympy.Symbol("surface"):
        print("surface found")
        string_as_integral_string = "(" + replaced_as_string + ") * ds"

    return string_as_integral_string

def convert(sympy_equation: sympy.Eq, symbolic_trial_function: sympy.Symbol, symbolic_test_function: sympy.Symbol, trial_function_name: str, test_function_name: str):
    # replace symbolic test-and trial-functions u(V), v(V) with unknown symbolic functions u, v
    utils.print_space("Convert to UFL")
    sympy.pprint(sympy_equation)
    lhs = sympy_equation.lhs
    rhs = sympy_equation.rhs
    rhs_ufl_string = ""
    # all atoms of the equation should be integrals
    rhs_integrals = rhs.atoms(sympy.Integral)

    for index, integral in enumerate(rhs_integrals):
        # replace the sympy integrals with 
        rhs_ufl_string = rhs_ufl_string + ("" if index == 0 else " + ") + convert_integral_to_ufl_string(integral, symbolic_trial_function, symbolic_test_function, trial_function_name, test_function_name)

    lhs_ufl_string = ""
    # all atoms of the equation should be integrals
    lhs_integrals = lhs.atoms(sympy.Integral)
    for index, integral in enumerate(lhs_integrals):
        lhs_ufl_string = lhs_ufl_string + ("" if index == 0 else " + ") + convert_integral_to_ufl_string(integral, symbolic_trial_function, symbolic_test_function, trial_function_name, test_function_name)
    
    #TODO sort equation here to lhs -> u rhs -> no u (bilinear form)
    utils.print_space("LHS")
    print(lhs_ufl_string)
    utils.print_space("RHS")
    print(rhs_ufl_string)
    return lhs_ufl_string, rhs_ufl_string
    