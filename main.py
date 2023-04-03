import Preprocessing.preprocessing as pre
import Main.transformations as tr
import sympy
from sympy.vector import CoordSys3D, Laplacian
import Util.utils as utils


def solve(sympy_equation: sympy.Eq, test_function_v: sympy.Function):
    with_test_function = pre.multiply_with_test_function(sympy_equation, test_function_v)
    integrals, omega = pre.integrate_over_domain(with_test_function)
    partially_integrated = tr.integrate_by_parts(integrals, test_function_v, omega)
    utils.print_space("Result")
    sympy.pprint(partially_integrated)
    return partially_integrated