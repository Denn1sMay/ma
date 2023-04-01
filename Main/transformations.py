import sympy
from sympy.vector import CoordSys3D, Laplacian, Del
import Util.utils as utils


def perform_integration_by_parts(integral: sympy.Integral, test_function_v: sympy.Function, omega: sympy.Symbol):
    utils.print_space("Performing integration by parts on current integral")
    sympy.pprint(integral)
    # Get the integrated function
    integral_args = integral.args[0]
    domain = integral.args[1]

    # Get the laplacian part of it
    a = integral_args.atoms(Laplacian)
    #TODO equation can contain multiple laplace operators
    laplacian_function = min(a)

    # Get the trial function (arguments/ parameters of the laplacian function)
    laplacian_args = laplacian_function.args
    trial_function_u = laplacian_args[0]
    print("First argument of Laplace Operator: ")
    print(trial_function_u)

    # Define the gradient of the trial function (argument of the laplace function)
    nabla_grad = sympy.Function("nabla_grad")

    nabla_trial_function_u = nabla_grad(trial_function_u)
    # Define the nabla operator on the test function
    nabla_test_function_v = nabla_grad(test_function_v)

    # Perform the integration by parts 
    #  -> replace the multiplication of laplace(u) * v
    trial_test_mult = laplacian_function * test_function_v
    inner_func = sympy.Function("inner_func")(nabla_trial_function_u, nabla_test_function_v)
    
    res_of_integration_by_parts = integral_args.subs(trial_test_mult, inner_func)
    integrated_parts = sympy.Integral(res_of_integration_by_parts, omega)
    print("")
    print("Transformed Integral: ")
    sympy.pprint(integrated_parts)
    return integrated_parts


def integrate_by_parts(sympy_equation: sympy.Eq, test_function_v: sympy.Function, omega: sympy.Symbol):
    all_integrals = sympy_equation.atoms(sympy.Integral)


    lhs_complete = sympy_equation.lhs
    lhs_integrals = lhs_complete.atoms(sympy.Integral)
    
    rhs_complete = sympy_equation.rhs
    rhs_integrals = rhs_complete.atoms(sympy.Integral)
    
    sympy.pprint(all_integrals)

    lhs = 0
    for i in lhs_integrals:
        if i.has(Laplacian):
            print("Laplacian found in current i")
            new_integral = perform_integration_by_parts(i, test_function_v, omega)
            lhs = lhs + new_integral
        else:
            lhs = lhs + i


    rhs = 0
    for i in rhs_integrals:
        if i.has(Laplacian):
            print("Laplacian found in current i")
            new_integral = perform_integration_by_parts(i, test_function_v, omega)
            rhs = rhs + new_integral
        else:
            rhs = rhs + i

    new_eq = sympy.Eq(lhs, rhs)

    return new_eq
