import Preprocessing.preprocessing as pre
import Main.transformations as transform
import UFL_Convert.conversions as cv
import sympy
from sympy.vector import CoordSys3D, Laplacian
import Util.utils as utils
from ufl import *
from typing import Optional
from Util.Boundaries import Boundaries


def get_weak_form(sympy_equation: sympy.Eq, trial_function_name: str, test_function_name: str, string_equation: Optional[str] = None, boundary_condition: Optional[Boundaries] = None, boundary_function: Optional[str] = None):
    surface = None
    if string_equation != None:
        sympy_equation = pre.parse_string_equation(string_equation)
    if(boundary_function != None):
        boundary_function = sympy.Symbol(boundary_function)
        surface = sympy.Symbol("surface")
    test_function_v = sympy.Symbol(test_function_name)
    trial_function_u = sympy.Symbol(trial_function_name)
    with_test_function = pre.multiply_with_test_function(sympy_equation, test_function_v)
    integrals, omega = pre.integrate_over_domain_and_sort(with_test_function, trial_function_u)
    partially_integrated = transform.integrate_by_parts(integrals, test_function_v, omega, boundary_condition, boundary_function, surface)
    partially_integrated_sorted = pre.sort_equation(partially_integrated.lhs, partially_integrated.rhs, trial_function_u)

    utils.print_space("Result")
    sympy.pprint(partially_integrated_sorted)
    a, L = cv.convert(partially_integrated_sorted, trial_function_u, test_function_v, trial_function_name, test_function_name)
    return a, L, partially_integrated_sorted