md"""
# Solving Ordinary Differential Equations (ODEs)

- [DifferentialEquations.jl docs](https://diffeq.sciml.ai/dev/index.html)

## General steps to solve differential equations

- Define a model function representing the right-hand-side (RHS) of the sysstem.
  - Out-of-place form: `f(u, p, t)` where `u` is the state variable(s), `p` is the parameter(s), and `t` is the independent variable (usually time). The output is the right hand side (RHS) of the differential equation system.
  - In-place form: `f!(du, u, p, t)`, where the output is saved to `du`. The rest is the same as the out of place form. The in-place form has potential performance benefits since it allocates less arrays than the out-of-place form.
- Initial conditions (`u0`) for the state variable(s).
- (Optional) parameter(s) `p`.
- Define a problem (e.g. `ODEProblem`) using the modeling function (`f`), initial conditions (`u0`), simulation time span (`tspan == (tstart, tend)`), and parameter(s) `p`.
- Solve the problem by calling `solve(prob)`.
"""

md"""
## Radioactive decay

As a simple example, the decaying rate of a nuclear isotope is proprotional to its concentration:

$$
\frac{d}{dt}C(t) = - \lambda C(t)
$$

**State variable(s)**

- $C(t)$: The concentration of a decaying nuclear isotope.

**Parameter(s)**

- $\lambda$: The rate constant of decay. The half-life $t_{\frac{1}{2}} = \frac{ln2}{\lambda}$
"""

using DifferentialEquations
using Plots
PNG(img) = display("image/png", img) ## Force PNG output

# The exponential decay ODE model, out-of-place (3-parameter) form
expdecay(u, p, t) = p * u

#---
p = -1.0 ## Parameter
u0 = 1.0 ## Initial condition
tspan = (0.0, 2.0) ## Simulation start and end time points
prob = ODEProblem(expdecay, u0, tspan, p) ## Define the problem
sol = solve(prob) ## Solve the problem

# ### Visualize
# Using Plots.jl package.
fig = plot(sol);
fig |> PNG

# Solution handling options, see: https://diffeq.sciml.ai/stable/basics/solution/
# `sol(t)`: the solution at time `t` with interpolation.
sol(1.0)

# `sol.t`: time points of the solution. Notice *t=1* may not in one of the time points.
sol.t

# `sol.u`: the state variables at `sol.t`
sol.u

md"""
## The SIR model

A more complicated example is the [SIR model](https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model) describing infectious disease spreading. There are more state variables and parameters.

$$
\begin{align}
\frac{d}{dt}S(t) &= - \beta S(t)I(t)  \\
\frac{d}{dt}I(t) &= \beta S(t)I(t)  - \gamma I(t)  \\
\frac{d}{dt}R(t) &= \gamma I(t)
\end{align}
$$

**State variable(s)**

- $S(t)$ : the fraction of susceptible people
- $I(t)$ : the fraction of infectious people
- $R(t)$ : the fraction of recovered (or removed) people

**Parameter(s)**

- $\beta$ : the rate of infection when susceptible and infectious people meet
- $\gamma$ : the rate of recovery of infectious people
"""

using DifferentialEquations
using Plots

# SIR model, using the in-place form: f!(du, u, p ,t)
function sir!(du, u, p, t)
	s, i, r = u
	β, γ = p
	v1 = β * s * i
	v2 = γ * i
    du[1] = -v1
    du[2] = v1 - v2
    du[3] = v2
	return nothing
end

# Parameters, initial conditions, and time span of the SIR model
p = (β = 1.0, γ = 0.3)
u0 = [0.99, 0.01, 0.00]
tspan = (0.0, 20.0)

# Define the problem
prob = ODEProblem(sir!, u0, tspan, p)

# Solve the problem
sol = solve(prob)

# Visualize the solution
fig = plot(sol, labels=["S" "I" "R"], legend=:right);
fig |> PNG

# `sol[i]`: all components at timestep `i`
sol[2]

# `sol[i, j]`: `i`th component at timestep `j`
sol[1, 2]

# `sol[i, :]`: the timeseries for the `i`th component.
sol[1, :]

# `sol(t,idxs=1)`: the 1st element in time point(s) `t` with interpolation. `t` can be a scalar (single point) or an vector-like sequence. (multiple points)
sol(10, idxs=2)

#---
sol(0.0:0.1:20.0, idxs=2)

md"""
## Lorenz system

The Lorenz system is a system of ordinary differential equations having chaotic solutions for certain parameter values and initial conditions. ([Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system))

$$
\begin{align}
  \frac{dx}{dt} &=& \sigma(y-x) \\
  \frac{dy}{dt} &=& x(\rho - z) -y \\
  \frac{dz}{dt} &=& xy - \beta z
\end{align}
$$

In this example, we will use [LabelledArrays.jl](https://github.com/SciML/LabelledArrays.jl) to get DSL-like syntax.
"""
using LabelledArrays
using DifferentialEquations
using Plots

#---
function lorenz!(du,u,p,t)
    du.x = p.σ*(u.y-u.x)
    du.y = u.x*(p.ρ-u.z) - u.y
    du.z = u.x*u.y - p.β*u.z
end

#---
u0 = LVector(x=1.0, y=0.0, z=0.0)
p = (σ=10.0, ρ=28.0, β=8/3)
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob)

# `idxs=(1, 2, 3)` makes a phase plot with 1st, 2nd, and the 3rd state variable. With `LabelledArrays`, you can also use symbols instead of index numbers.

## Same as plot(sol, idxs=(1, 2, 3))
fig = plot(sol, idxs=(:x, :y, :z));
fig |> PNG

# The plot recipe is using the interpolation to smooth the plot. You can turn off `denseplot` to see the difference.
fig = plot(sol, idxs=(1, 2, 3), denseplot=false);
fig |> PNG

# The zeroth variable in `idxs` is the independent variable (usually time). The below command plots the time series of the second state variable (`y`).
fig = plot(sol, idxs=(0, 2));
fig |> PNG

md"""
## Non-autonomous ODEs

In non-autonomous ODEs, Some term(s) in the right-hadn-side (RHS) are time-dependent. For example, in this pendulum model with external, time-dependent, force.

$$
\begin{aligned}
\dot{\theta} &= \omega(t) \\
\dot{\omega} &= -1.5\frac{g}{l}sin(\theta(t)) + \frac{3}{ml^2}M(t)
\end{aligned}
$$

- $\theta$: pendulum angle
- $\omega$: angular rate
- M: time-dependent external torgue
- $l$: pendulum length
- $g$: gravitional acceleration
"""

using DifferentialEquations
using Plots

l = 1.0                             ## length [m]
m = 1.0                             ## mass [kg]
g = 9.81                            ## gravitational acceleration [m/s²]

#---
function pendulum!(du, u, p, t)
    du[1] = u[2]                    ## θ'(t) = ω(t)
    du[2] = -3g/(2l)*sin(u[1]) + 3/(m*l^2)*p(t) # ω'(t) = -3g/(2l) sin θ(t) + 3/(ml^2)M(t)
end

θ₀ = 0.01                           ## initial angular deflection [rad]
ω₀ = 0.0                            ## initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       ## initial state vector
tspan = (0.0,10.0)                  ## time interval

M = t->0.1sin(t)                    ## external torque [Nm] as the parameter for the pendulum model

prob = ODEProblem(pendulum!, u₀, tspan, M)
sol = solve(prob)

fig = plot(sol, linewidth=2, xaxis="t", label=["θ [rad]" "ω [rad/s]"], layout=(2,1));
fig |> PNG

md"""
## Linear ODE system

In this example, the ODE system is described by a matrix differential operator.

$\dot{u} = Au$
"""

using DifferentialEquations
using Plots

A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]

u0 = rand(4, 2)

tspan = (0.0, 1.0)
f = (u, p, t) -> A*u
prob = ODEProblem(f, u0, tspan)
sol = solve(prob)
fig = plot(sol);
fig |> PNG
