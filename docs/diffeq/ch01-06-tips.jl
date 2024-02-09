#===
# Tips and tricks

The following section contains tips and tricks in modeling and simulation.
===#
using Plots
PNG(img) = display("image/png", img) ## Force PNG output

#===
## Smooth and differentiable switch functions

Using a similar smooth function for discontinuity is more friendly to ODE solvers.

### Heaviside step function

A [Heaviside step function](https://en.wikipedia.org/wiki/Heaviside_step_function) (0 when x < a, 1 when x > a) could be approximated with a steep [logistic function](https://en.wikipedia.org/wiki/Logistic_function).
===#

logistic(x, k=100) = inv(1 + exp(-k * x))

# The function output switches from zero to one around `x=0`
plot(logistic, -1, 1, label="sigmoid approx.") |> PNG

# ### Single pulse
# A single pulse could be approximated with a product of two logistic functions
singlepulse(x, t0=0, t1=0.1, k=100) = logistic(x - t0, k) * logistic(t1 - x, k)

plot(singlepulse, -1, 1, label="sigmoid approx.") |> PNG

# ### Smooth minimal function
# From this discourse post: https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5
function smoothmin(x, k=100)
    ex = exp(-k)
    ey = exp(-k*(1+x))
    return (ex + (1+x)*ey)/(ex+ey)
end

plot(smoothmin, -1, 1, label="Exponential approx.") |> PNG

# ### Periodic pulses
# From: https://www.noamross.net/2015/11/12/a-smooth-differentiable-pulse-function/

function smoothpulses(t, tstart, tend, period=1, amplitude=period / (tend - tstart), steepness=1000)
    @assert tstart < tend < period
    xi = 3 / 4 - (tend - tstart) / (2 * period)
    p = inv(1 + exp(steepness * (sinpi(2 * ((t - tstart) / period + xi)) - sinpi(2 * xi))))
    return amplitude * p
end

plot(t->smoothpulses(t, 0.2, 0.3, 0.5), 0.0, 2.0, lab="pulses") |> PNG

#===
## Avoiding DomainErrors

Some functions like `sqrt(x)`, `log(x)`, and `pow(x)` will throw `DomainError` exceptions with negative `x`, interrupting differential equation solvers. One can use the respective functions in https://github.com/JuliaMath/NaNMath.jl, returning `NaN` instead of throwing a `DomainError`. Then, the differential equation solvers will reject the solution and retry with a smaller time step.
===#

import NaNMath as nm
nm.sqrt(-1.0) ## returns NaN
