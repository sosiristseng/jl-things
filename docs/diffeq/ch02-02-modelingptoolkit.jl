md"""
# ModelingToolkit.jl

- [Simulating Big Models in Julia with ModelingToolkit @ JuliaCon 2021 Workshop](https://youtu.be/HEVOgSLBzWA)
- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl): Symbolic representations for modeling numerical systems.

## First example
"""

using ModelingToolkit
using DifferentialEquations
using Plots
using LinearAlgebra
PNG(img) = display("image/png", img) ## Force PNG output

# Independent (time) and dependent (state) variables (x and RHS)
@variables t x(t) RHS(t)

# Setting parameters in the modeling
@parameters τ

# Differential operator w.r.t. time
D = Differential(t)

# Equations in MTK use the tilde character (`~`) as equality.
# Every MTK system requires a name. The `@named` macro simply ensures that the symbolic name matches the name in the REPL.
@named fol_separate = ODESystem([
    RHS  ~ (1 - x)/τ,
    D(x) ~ RHS
])

# `structural_simplify()` transforms simple DAEs with dependent terms to ODEs and reduces the number of state variables.

model = structural_simplify(fol_separate)
u0 = [x => 0.0]
tspan = (0.0, 10.0)
param = [τ => 3.0]
prob = ODEProblem(model, u0, tspan, param)
sol = solve(prob)

# The eliminated term (RHS in this example) is still tracible
fig = plot(sol, idxs=[x, RHS], legend=:right);
fig |> PNG

# ## Time-variant external force
# If the function is too complex and/or has discontinuity, one should apply `@register_symbolic` to the function to exclude it from symbolic transformations and use it as-is.
@variables t x(t) f(t)
@parameters τ
D = Differential(t)

value_vector = randn(10)

# Define a time-dependent random external force
f_fun(t) = t >= 10 ? value_vector[end] : value_vector[Int(floor(t))+1]

# "Register" arbitrary Julia functions to be excluded from symbolic transformations. Just use it as-is.
@register_symbolic f_fun(t)
@named fol_external_f = ODESystem([f ~ f_fun(t), D(x) ~ (f - x)/τ])

#---
prob = ODEProblem(structural_simplify(fol_external_f), [x => 0.0], (0.0,10.0), [τ => 0.75])
sol = solve(prob)
fig = plot(sol, idxs=[x, f]);
fig |> PNG

# ## Second order ODE system
# `ode_order_lowering(sys)` can automatically transform a second-order ODE into two first-order ODEs.

using Plots
using ModelingToolkit
using DifferentialEquations

@parameters σ ρ β
@variables t x(t) y(t) z(t)
D = Differential(t)

#---
eqs = [
    D(D(x)) ~ σ * (y-x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]

#---
@named sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

#---
u0 = [
    D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0
]

p = [
    σ => 28.0,
    ρ => 10.0,
    β => 8/3
]

tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, p, jac=true)
sol = solve(prob)
fig = plot(sol, idxs=(x, y, z), label="Trajectory", size=(600,600));
fig |> PNG

# ## Composing systems
# By defining connection equation(s) to couple ODE systems together, we can build component-based, hierarchical models.

using Plots
using ModelingToolkit
using DifferentialEquations

@parameters σ ρ β
@variables t x(t) y(t) z(t)

D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
       D(y) ~ x * (ρ - z) - y,
       D(z) ~ x * y - β * z]

@named lorenz1 = ODESystem(eqs)
@named lorenz2 = ODESystem(eqs)

# Define relations between the two systems
@variables a(t)
@parameters γ

connections = [0 ~ lorenz1.x + lorenz2.y + a * γ]

@named connLorenz = ODESystem(connections, t , [a], [γ], systems = [lorenz1, lorenz2])

# All state variables in the combined system
states(connLorenz)

#---
u0 = [
    lorenz1.x => 1.0, lorenz1.y => 0.0, lorenz1.z => 0.0,
    lorenz2.x => 0.0, lorenz2.y => 1.0, lorenz2.z => 0.0,
    a => 2.0
]

p = [
    lorenz1.σ => 10.0, lorenz1.ρ => 28.0, lorenz1.β => 8/3,
    lorenz2.σ => 10.0, lorenz2.ρ => 28.0, lorenz2.β => 8/3,
    γ => 2.0
]

tspan = (0.0, 100.0)
sol = solve(ODEProblem(connLorenz, u0, tspan, p, jac=true))
fig = plot(sol, idxs=(a, lorenz1.x, lorenz2.x), size=(600,600));
fig |> PNG

md"""
## Convert existing functions into MTK ones

`modelingtoolkitize(prob)`

And it can generate analytic Jacobin function for faster solving.

Example: **[DAE index reduction](https://mtk.sciml.ai/stable/mtkitize_tutorials/modelingtoolkitize_index_reduction/)** for the pendulum problem, which cannot be solved by regular ODE solvers.
"""

using Plots
using ModelingToolkit
using DifferentialEquations

function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx
    du[2] = T*x
    du[3] = dy
    du[4] = T*y - g
    ## Do not write your function like this after you've learned MTK
    du[5] = x^2 + y^2 - L^2
    return nothing
end

#---
pendulum_fun! = ODEFunction(pendulum!, mass_matrix = Diagonal([1, 1, 1, 1, 0]))
u0 = [1.0, 0.0, 0.0, 0.0, 0.0]
p = [9.8, 1.0]
tspan = (0.0, 10.0)
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)

# Convert the ODE problem to a MTK system.
tracedSys = modelingtoolkitize(pendulum_prob)

# `structural_simplify()` and `dae_index_lowering()` transform the index-3 DAE into an index-0 ODE.
pendulumSys = tracedSys |> dae_index_lowering |> structural_simplify
# The default `u0` is included in the system already so one can use an empty array `[]` as the initial conditions.
prob = ODAEProblem(pendulumSys, [], tspan)
sol = solve(prob, abstol=1e-8, reltol=1e-8)
fig = plot(sol, idxs=states(tracedSys));
fig |> PNG

fieldnames(typeof(sol))

# ## Solving nonlinear systems
# Use `NonlinearSolve.jl` and `NonlinearSystem()`

using ModelingToolkit
using NonlinearSolve
@variables x y z
@parameters σ ρ β

eqs = [
    0 ~ σ * (y-x),
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z
]

#---
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])

#---
guess = [x => 1.0, y => 0.0, z => 0.0]
ps = [σ => 10.0, ρ => 26.0, β => 8/3]
prob = NonlinearProblem(ns, guess, ps)
sol = solve(prob, NewtonRaphson()) ## Should be all zeroes

# Another nonlinear system example with `structural_simplify()`
@parameters t
@variables u1(t) u2(t) u3(t) u4(t) u5(t)

eqs = [
    0 ~ u1 - sin(u5)
    0 ~ u2 - cos(u1)
    0 ~ u3 - hypot(u1, u2)
    0 ~ u4 - hypot(u2, u3)
    0 ~ u5 - hypot(u4, u1)
]

@named sys = NonlinearSystem(eqs, [u1, u2, u3, u4, u5], [])

# You can simplify the problem using `structural_simplify()`
# There will be only one state variable left. The solve can solve the problem faster.
simple_sys = structural_simplify(sys)

prob = NonlinearProblem(simple_sys, [u5=>0.0])
sol = solve(prob, NewtonRaphson())

# The answer should be 1.6 and 1.0
@show sol[u5] sol[u1];  ## 1.6 and 1.0

# ## Stochastic Differential Equations (SDEs)
# Use `SDESystem(equations, noises, iv, dv, ps)`
using ModelingToolkit
using DifferentialEquations
using Plots

@parameters σ ρ β
@variables t x(t) y(t) z(t)

D = Differential(t)

eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]

# Add a diagonal noise with 10% of the magnitude.
noises = [0.1x, 0.1y, 0.1z]

@named sys = SDESystem(eqs, noises, t, [x, y, z], [σ,ρ,β])

#---
u0 = [x => 10.0, y => 10.0, z=> 10.0]
p = [σ => 10.0, ρ => 28.0, β => 8/3]
tspan = (0.0, 200.0)
prob = SDEProblem(sys, u0, tspan, p)
sol = solve(prob)
fig = plot(sol, idxs=(1, 2, 3));
fig |> PNG
