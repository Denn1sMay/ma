import sympy
from sympy.vector import CoordSys3D, Laplacian, Del
import Util.utils as utils
from enum import Enum
from Util.Boundaries import Boundaries
from typing import Optional

def perform_integration_by_parts_on_integral(integral: sympy.Integral, test_function_v: sympy.Symbol, omega: sympy.Symbol, boundary_condition: Optional[Boundaries] = None, boundary_function: Optional[sympy.Symbol] = None, surface: Optional[sympy.Symbol] = None):
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
    print("Argument of Laplace Operator: ")
    #TODO argument of laplace operator can be a complex expression
    print(trial_function_u)

    # Define the gradient of the trial function (argument of the laplace function)
    nabla_grad = sympy.Function("grad")

    nabla_trial_function_u = nabla_grad(trial_function_u)
    # Define the nabla operator on the test function
    nabla_test_function_v = nabla_grad(test_function_v)

    # Perform the integration by parts 
    #  -> replace the multiplication of laplace(u) * v
    laplacian_test_mult = laplacian_function * test_function_v
    inner_func = sympy.Function("inner")(nabla_trial_function_u, nabla_test_function_v)
    print("test")
    print("------------")
    print(integral_args)
    integral_over_domain = integral_args.subs(laplacian_test_mult, inner_func)
    print(integral_args)
    integrated_parts = sympy.Integral(-integral_over_domain, omega)
    if boundary_condition != None:
        if boundary_condition == Boundaries.neumann and boundary_function == None:
            raise Exception("Need to provide a boundary function symbol to use neumann boundaries")
        if boundary_condition == Boundaries.neumann:
            print("")
            integral_over_boundary = integral_args.subs(laplacian_test_mult, (boundary_function * test_function_v))
            integrated_parts = integrated_parts + sympy.Integral(integral_over_boundary, surface)


    print("")
    print("Transformed Integral: ")
    sympy.pprint(integrated_parts)
    return integrated_parts


def check_linearity(term):
    utils.print_space("Checking Linearity")
    print(type(term))
    exponential = term.atoms(sympy.Pow)
    if len(exponential) > 0:    
        for exponential_arg in exponential:
            laplace_in_exponential = exponential_arg.atoms(Laplacian)
            if len(laplace_in_exponential) > 0:
                print("Nonlinear PDE - path not programmed yet")
                raise Exception("Nonlinear PDE")

def integrate_by_parts(sympy_equation: sympy.Eq, test_function_v: sympy.Symbol, omega: sympy.Symbol, boundary_condition: Optional[Boundaries] = None, boundary_function: Optional[sympy.Symbol] = None, surface: Optional[sympy.Symbol] = None):
    all_integrals = sympy_equation.atoms(sympy.Integral)


    lhs_complete = sympy_equation.lhs
    lhs_integrals = lhs_complete.atoms(sympy.Integral)
    lhs_validation = sympy.Eq(sympy.Add(*lhs_integrals), lhs_complete)
    print("Validation")
    print(lhs_validation)
    if lhs_validation == False:
        print("current LHS:")
        sympy.pprint(lhs_complete)
        print("Sum of LHS Integrals:")
        sympy.pprint(sympy.Add(*lhs_integrals))
        raise Exception("Calculated Integrals do not match input equation")

    rhs_complete = sympy_equation.rhs
    rhs_integrals = rhs_complete.atoms(sympy.Integral)

    rhs_validation = sympy.Eq(sympy.Add(*rhs_integrals), rhs_complete)
    if rhs_validation == False:
        print("current RHS:")
        sympy.pprint(rhs_complete)
        print("Sum of RHS Integrals:")
        sympy.pprint(sympy.Add(*rhs_integrals))
        raise Exception("Calculated Integrals do not match input equation")
    
    utils.print_space("All Integrals of Equation: ")
    sympy.pprint(all_integrals)
    #TODO Sort equation by trial function u to lhs (a), other functions to rhs (L)
    lhs = 0
    for i in lhs_integrals:
        if i.has(Laplacian):
            print("Laplacian found in current integral:")
            sympy.pprint(i)
            check_linearity(i)
            new_integral = perform_integration_by_parts_on_integral(i, test_function_v, omega, boundary_condition, boundary_function, surface)
            lhs = lhs + new_integral
        else:
            lhs = lhs + i


    rhs = 0
    for i in rhs_integrals:
        if i.has(Laplacian):
            print("Laplacian found in current integral: ")
            sympy.pprint(i)
            check_linearity(i)
            new_integral = perform_integration_by_parts_on_integral(i, test_function_v, omega)
            rhs = rhs + new_integral
        else:
            rhs = rhs + i
            
    new_eq = sympy.Eq(lhs, rhs)

    return new_eq
