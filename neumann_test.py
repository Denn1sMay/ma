import numpy as np
import pyvista
from dolfinx import mesh, io, plot

from dolfinx.fem import (Constant, Function, FunctionSpace, 
                         assemble_scalar, dirichletbc, form, locate_dofs_geometrical)
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from ufl import SpatialCoordinate, TestFunction, TrialFunction, dot, ds, dx, grad, inner

from Util.Boundaries import Boundaries
import Util.utils as utils
import main


mesh = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)
V = FunctionSpace(mesh, ("CG", 1))
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v)) * dx

def u_exact(x):
    return 1 + x[0]**2 + 2*x[1]**2

def boundary_D(x):
    return np.logical_or(np.isclose(x[0], 0), np.isclose(x[0],1))

dofs_D = locate_dofs_geometrical(V, boundary_D)
u_bc = Function(V)
u_bc.interpolate(u_exact)
bc = dirichletbc(u_bc, dofs_D)

x = SpatialCoordinate(mesh)
g = -4 * x[1]
f = Constant(mesh, ScalarType(-6))
L = f * v * dx - g * v * ds

a_gen, L_gen, s = main.get_weak_form(None, "u", "v", "-Laplacian(u) = f", Boundaries.neumann, "g")

utils.print_space("Original Form:")
print("lhs")
print(a)
print("rhs")
print(L)
utils.print_space("Transformed Form:")
print("lhs")
print(eval(a_gen))
print("rhs")
print(eval(L_gen))


problem = LinearProblem(eval(a_gen), eval(L_gen), bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

V2 = FunctionSpace(mesh, ("CG", 2))
uex = Function(V2)
uex.interpolate(u_exact)
error_L2 = assemble_scalar(form((uh - uex)**2 * dx))
error_L2 = np.sqrt(MPI.COMM_WORLD.allreduce(error_L2, op=MPI.SUM))

u_vertex_values = uh.x.array
uex_1 = Function(V)
uex_1.interpolate(uex)
u_ex_vertex_values = uex_1.x.array
error_max = np.max(np.abs(u_vertex_values - u_ex_vertex_values))
error_max = MPI.COMM_WORLD.allreduce(error_max, op=MPI.MAX)
print(f"Error_L2 : {error_L2:.2e}")
print(f"Error_max : {error_max:.2e}")



# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
with io.XDMFFile(mesh.comm, "out/poisson.xdmf", "w") as file:
    file.write_mesh(mesh)
    file.write_function(uh)

try:
    import pyvista
    cells, types, x = plot.create_vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(cells, types, x)
    grid.point_data["u"] = uh.x.array.real
    grid.set_active_scalars("u")
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    warped = grid.warp_by_scalar()
    plotter.add_mesh(warped)
    if pyvista.OFF_SCREEN:
        pyvista.start_xvfb(wait=0.1)
        plotter.screenshot("uh_poisson.png")
    else:
        plotter.show()

except ModuleNotFoundError:
    print("'pyvista' is required to visualise the solution")
    print("Install 'pyvista' with pip: 'python3 -m pip install pyvista'")

