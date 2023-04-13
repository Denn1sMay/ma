import Preprocessing.preprocessing as pre
import Main.transformations as transform
import UFL_Convert.conversions as cv
import sympy
from sympy.vector import CoordSys3D, Laplacian
import Util.utils as utils
from ufl import *


def get_weak_form(sympy_equation: sympy.Eq, trial_function_name: str, test_function_name: str, string_equation=None):
    if string_equation != None:
        sympy_equation = pre.parse_string_equation(string_equation)
    test_function_v = sympy.Function(test_function_name)(sympy.Symbol("V"))
    trial_function_u = sympy.Function(trial_function_name)(sympy.Symbol("V"))
    with_test_function = pre.multiply_with_test_function(sympy_equation, test_function_v)
    integrals, omega = pre.integrate_over_domain(with_test_function)
    partially_integrated = transform.integrate_by_parts(integrals, test_function_v, omega)


    utils.print_space("Result")
    sympy.pprint(partially_integrated)
    a, L = cv.convert(partially_integrated, trial_function_u, test_function_v, trial_function_name, test_function_name)
    return a, L, partially_integrated