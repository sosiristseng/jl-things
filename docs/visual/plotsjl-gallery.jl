md"""
# Plots.jl gallery

Sources:

- [PlotsGallery.jl](https://github.com/goropikari/PlotsGallery.jl) by goropikari
- [Plots.jl docs](https://docs.juliaplots.org/)
"""

using Plots
using Random
Random.seed!(2021)
PNG(img) = display("image/png", img) ## Force PNG output

# ## Attributes
# https://docs.juliaplots.org/stable/attributes/
# ### Ticks size and properties
fig = plot(
    sin, 0, 2π;
    xticks=0:0.5:2π,
    xrotation=60,
    xtickfontsize=25,
    bottom_margin=15Plots.mm
);
fig |> PNG

#---
fig = plot(
    sin, 0, 2π;
    xtick=(0:0.5:2π, ["$i a" for i in 0:0.5:2π]),
    ytick=-1:0.2:1,
    xrotation=60,
    yrotation=90,
);
fig |> PNG

# ### No axis
# `axis=false`
fig = plot(sin, 0, 2π, axis=false);
fig |> PNG

# ### Log scale for axes
# `xscale=:log10`, `yscale=:log10`
fig = plot(exp, -5, 5, yscale=:log10, title="semilogy", legend=nothing);
fig |> PNG

#---
fig = plot(log, 1e-5, 10, xscale=:log10, title="semilogx", legend=nothing);
fig |> PNG

#---
fig = plot(x->x^1.7, 1e-3, 3, scale=:log10, title="log-log", legend=nothing);
fig |> PNG

# ### Axis range
# `xlims` and `ylims`

fig = plot(sin, 0, 2π, xlims=(-10, 10), ylims=(-2,2));
fig |> PNG

#---
fig = plot(exp, 0, 10, yformatter=:scientific);
fig |> PNG

# ### Flip Axis
# `xflip=true` and/or `yflip=true`

fig = plot(identity, 0:0.01:2π, proj=:polar, xflip=true, yflip=true);
fig |> PNG

# ### Aspect ratio
# `aspect_ratio=:equal` or `aspect_ratio=<number>`
fig = heatmap(bitrand(10, 10), aspect_ratio=:equal, c=:blues, colorbar=false);
fig |> PNG

# ### Fonts
# LaTeX fonts are supported by the `LaTeXStrings.jl` package.

using Plots
using LaTeXStrings

fig = plot(sin, 0, 2π,
    title=L"y = \sin(x)",
    titlefont=font(40), ## title

    xlabel=L"x",
    ylabel="y",
    xguidefontsize=30, ## x-guide
    yguidefontsize=20, ## y-guide
    ## guidefontsize=20, ## both x,y-guide

    xtick=(0:0.5:2π, ["\$ $(i) \$" for i in 0:0.5:2π]),
    ytick=-1:0.5:1,
    xtickfontsize=15,
    ytickfontsize=20,
    ## tickfontsize=10, ## for both x and y

    label="Sin function",
    legendfontsize=12,

    xlims=(0,2π),
    ylims=(-1,1),
    bottom_margin=5Plots.mm,
    left_margin=10Plots.mm,
    top_margin=15Plots.mm
);
fig |> PNG

fib(x) = (((1+sqrt(5))/2)^x - ((1-sqrt(5))/2)^x)/sqrt(5)

ann = L"F_n = \frac{1}{\sqrt{5}} \left[\left( \frac{1+\sqrt{5}}{2} \right)^n - \left( \frac{1-\sqrt{5}}{2} \right)^n \right]"

fig = plot(fib, 1:12, marker=:circle, xlabel=L"n", ylabel=L"F_n", annotation=(5, 100, ann));
fig |> PNG

# ## Bar plots

using Plots
using StatsPlots
using StatsBase

# Prepare data
measles = [38556, 24472, 14556, 18060, 19549, 8122, 28541, 7880, 3283, 4135, 7953, 1884]
mumps = [20178, 23536, 34561, 37395, 36072, 32237, 18597, 9408, 6005, 6268, 8963, 13882]
chickenPox = [37140, 32169, 37533, 39103, 33244, 23269, 16737, 5411, 3435, 6052, 12825, 23332]
ticklabel = string.(collect('A':'L'))

md"""
### Grouped vertical bar plots

Requires the `StatsPlots.jl` package.
`groupedbar(data, bar_position = :dodge)`
"""

fig = groupedbar([measles mumps chickenPox], bar_position = :dodge, bar_width=0.7, xticks=(1:12, ticklabel), label=["measles" "mumps" "chickenPox"]);
fig |> PNG

md"""
### Stacked vertical bar plots

Requires `StatsPlots` package.
`groupedbar(data, bar_position = :stack)`
"""

fig = groupedbar([measles mumps chickenPox],
        bar_position = :stack,
        bar_width=0.7,
        xticks=(1:12, ticklabel),
        label=["measles" "mumps" "chickenPox"]
);
fig |> PNG

# ### Horizontal Bar Plot
# `bar(data, orientation=:h)`

fig = bar(1:12, orientation=:h, yticks=(1:12, ticklabel), yflip=true);
fig |> PNG

# ### Categorical Histogram Plot
status = ["Poor", "Fair", "Good", "Excellent"]
data = sample(status, Weights([1,1,2,2]), 100)
datamap = countmap(data)

fig = bar((x -> datamap[x]).(status), xticks=(1:4, status), legend=nothing);
fig |> PNG

# ## Histogram
# `histogram(data, bins=N)`

x = randn(1000)
y = randn(1000)
z = randn(1000)

fig = histogram(x, bins=20, alpha=0.4, label="A");
histogram!(fig, y, bins=20, alpha=0.4, label="B");
histogram!(fig, z, bins=20, alpha=0.4, label="C");
fig |> PNG

# ## Box plots
using Plots
using StatsPlots
using Statistics

n = 30
science = rand(1:10, n)
fig = boxplot(science, label="science");
fig |> PNG

english = rand(1:10, n)
fig = boxplot([science english], label=["science" "english"]);
fig |> PNG

md"""
## Contour Plots
### Over a function
`contour(xs, ys, f)` where `z = f(x, y)`
"""

using Plots

xs = range(0, stop=2, length=50)
ys = range(0, stop=2, length=50)
f = (x , y) -> x^2 + y^2

fig = contour(xs, ys, f);
fig |> PNG

md"""
### Nullclines
[Nullclines](https://en.wikipedia.org/wiki/Nullcline) (zero-growth isoclines) are curves where the derivative of one variable is zero. Nullclines are used to analyze evolution and stability of ODE systems.
"""

using Plots

dx = (x, y) -> 3x - 1.4x*y
dy = (x, y) -> -y + 0.8x*y

myrange = -1:0.01:4

fig = contour(myrange, myrange, dx, levels=[0], color=:red, legend=false);
contour!(fig, myrange, myrange, dy, levels=[0], color=:blue, legend=false);
fig |> PNG

# ### Contour plot over an array
# `contour(x1d, y1d, xy2d)`
## Notice that xs is transposed, making zz a 2D matrix
zz = f.(xs', ys)
fig = contour(xs, ys, zz);
fig |> PNG

md"""
### Filled Contour Plots
+ `contour(xs, ys, f, fill=true)`
+ `contourf(xs, ys, f)`
"""

fig = contour(0:0.01:5, 0:0.01:5, (x, y) -> sin(3x) * cos(x+y), xlabel="x", ylabel="y", fill=true);
fig |> PNG

md"""
## Datetime plot
- Use `Dates` package and `Data` data type
- Customize ticks
"""

using Plots
using Dates

days = 31
position = cumsum(randn(days))
x = Date(2018,1,1):Day(1):Date(2018,1,31)
ticks = [x[i] for i in 1:5:length(x)]

fig = plot(x, position,
     xlabel="Date",
     ylabel="Position",
     title="Track of random walker",
	 legend=nothing,
     xticks=ticks,
     xrotation=45,
     bottom_margin=10Plots.mm,
     left_margin=5Plots.mm
);
fig |> PNG

# ## Error bar
# `plots(..., xerr=xerr, yerr=yerr)`

using Plots

f = x -> 2 * x + 1
x = 0:0.1:2
n = length(x)
y = f.(x) + randn(n)

fig = plot(x, y,
    xerr=0.1 * rand(n),
    yerr=rand(n),
	legend=nothing
);
fig |> PNG

# ## Heatmap
# `heatmap(data)`
using Plots

a = rand(5,5)
xlabel = string.(collect('A':'E'))
ylabel = string.(collect('a':'e'))
fig = heatmap(a, xticks=(1:5, xlabel),
           yticks=(1:5, ylabel),
           aspect_ratio=:equal
);
fig |> PNG

# Add number annotations to plots
fontsize = 15
nrow, ncol = size(a)
ann = [(i,j, text(round(a[i,j], digits=2), fontsize, :white, :center)) for i in 1:nrow for j in 1:ncol]

annotate!(fig, ann, linecolor=:white);
fig |> PNG
