md"""
# Parameter estimation

`DiffEqParamEstim.jl` is not installed with `DifferentialEquations.jl`. You need to install it manually:

```julia
using Pkg
Pkg.add("DiffEqParamEstim")
using DiffEqParamEstim
```

## Estimate a single parameter from the data and the ODE model

Let's optimize the parameters of the Lotka-Volterra equation.
"""

using DifferentialEquations
using Plots
using DiffEqParamEstim
using Optim
PNG(img) = display("image/png", img) ## Force PNG output

#---
function lotka_volterra!(du, u, p, t)
    du[1] = dx = p[1]*u[1] - u[1]*u[2]
    du[2] = dy = -3*u[2] + u[1]*u[2]
end

u0 = [1.0;1.0]
tspan = (0.0, 10.0)
p = [1.5] ## The true value to be guessed
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Create a sample dataset with some noise.
ts = range(tspan[begin], tspan[end], 200)
data = [sol.(ts, idxs=1) sol.(ts, idxs=2)] .* (1 .+ 0.03 .* randn(length(ts), 2))

# Plotting the sample dataset and the true solution.
fig = plot(sol);
scatter!(fig, ts, data, label=["u1 data" "u2 data"]);
fig |> PNG

md"""
`DiffEqParamEstim.build_loss_objective()` builds a loss function for the ODE problem against the data.

We will minimize the mean squared error using `L2Loss()`.

Note that the data should be transposed.
"""

alg = Tsit5()
cost_function = build_loss_objective(prob, alg, L2Loss(collect(ts), transpose(data)), maxiters=10000, verbose=false)

fig = plot(
    cost_function, 0.0, 10.0,
    linewidth = 3, label=false, yscale=:log10,
    xaxis = "Parameter", yaxis = "Cost", title = "1-Parameter Cost Function"
);
fig |> PNG

# There is a dip (minimum) in the cost function at the true parameter value (1.5). We can use an optimizer, e.g., `Optim.jl`, to find the parameter value that minimizes the cost. (1.5 in this case)

result = Optim.optimize(cost_function, 0.0, 10.0)

# Note that the result is a vector because we used a different optimization algorithm

result.minimizer

# ## Estimate multiple parameters
# Let's use the Lotka-Volterra (Fox-rabbit) equations with all 4 parameters free.
function f2(du, u, p, t)
    du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
    du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end

u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
p = [1.5, 1.0, 3.0, 1.0]
prob = ODEProblem(f2, u0, tspan, p)
sol = solve(prob, Tsit5())

#---
ts = range(tspan[begin], tspan[end], 200)
data = [sol.(ts, idxs=1) sol.(ts, idxs=2)] .* (1 .+ 0.01 .* randn(length(ts), 2))

# Then we can find multiple parameters at once using the same steps.
cost_function = build_loss_objective(prob, Tsit5(), L2Loss(collect(ts), transpose(data)), maxiters=10000, verbose=false)
result_bfgs = Optim.optimize(cost_function, [1.3, 0.8, 2.8, 1.2], LBFGS())

# True parameters are `[1.5, 1.0, 3.0, 1.0]`
result_bfgs.minimizer
