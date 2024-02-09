md"""
# Heat Equation

Using `MethodOfLines.jl` (https://github.com/SciML/MethodOfLines.jl/) to sumbolically define the PDE system using the [finite difference method](https://en.wikipedia.org/wiki/Finite_difference_method) (FDM).

## 2D steady-state heat equation

From the [MethodOfLines tutorial](https://docs.sciml.ai/MethodOfLines/stable/tutorials/heatss/).

$$
\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0
$$
"""

using ModelingToolkit
using MethodOfLines
using DomainSets
using NonlinearSolve
using DifferentialEquations
using Plots
PNG(img) = display("image/png", img) ## Force PNG output

#---
@variables x y
@variables u(..)

Dxx = Differential(x)^2
Dyy = Differential(y)^2

# PDE equation
eq = Dxx(u(x, y)) + Dyy(u(x, y)) ~ 0

# Boundary conditions
bcs = [
    u(0, y) ~ x * y,
    u(1, y) ~ x * y,
    u(x, 0) ~ x * y,
    u(x, 1) ~ x * y
]

# Space and time domains
domains = [
    x ∈ Interval(0.0, 1.0),
    y ∈ Interval(0.0, 1.0)
]

# PDE system
@named pdesys = PDESystem([eq], bcs, domains, [x, y], [u(x, y)])

# Discretize the 2D sapce
# Note that we pass in `nothing` for the time variable here,
# since we are creating a stationary problem without a dependence on time, only on space.
N = 8
discretization = MOLFiniteDifference([x=>N, y=>N], nothing, approx_order=2, grid_align=edge_align, should_transform=false)

# It's a `NonlinearProblem` for solving a steady-state
prob = discretize(pdesys, discretization)

# Solve the PDE
sol = NonlinearSolve.solve(prob, NewtonRaphson())

# Visualize the heat equation solution
fig = heatmap(
    sol[x], sol[y], sol[u(x, y)],
    xlabel="x values", ylabel="y values",
    title="Steady State Heat Equation",
    aspect_ratio=:equal,
    xlims=(0.0, 1.0), ylims=(0.0, 1.0), clims=(0.0, 1.0)
);

fig |> PNG
