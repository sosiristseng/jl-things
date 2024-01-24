md"""
# Optimization Problems

From: https://mtk.sciml.ai/dev/tutorials/optimization/

Find $(x, y)$ that minimizes the loss function $(a - x)^2 + b(y - x^2)^2$
"""

using ModelingToolkit
using Optimization
using OptimizationOptimJL
PNG(img) = display("image/png", img) ## Force PNG output

#---
@variables x y
@parameters a b

# Target (loss) function
loss = (a - x)^2 + b * (y - x^2)^2

# The model
@named sys = OptimizationSystem(loss, [x, y], [a, b])

#---
u0 = [
    x=>1.0
    y=>2.0
]
p = [
    a => 6.0
    b => 7.0
]

# Generate Gradient and Hessian for the optimizer
prob = OptimizationProblem(sys, u0, p, grad=true, hess=true)

## The true solution is (6.0, 36.0)
sol = solve(prob, Newton())
