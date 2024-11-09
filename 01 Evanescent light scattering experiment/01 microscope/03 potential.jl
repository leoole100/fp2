using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit
using LinearAlgebra: norm, normalize
using CairoMakie
import CSV, Images
using StatsBase: fit, Histogram
using KernelDensity: kde
using Glob: glob
using MeanSquaredDisplacement: imsd
using Interpolations: interpolate, linear_interpolation
using Format: format

include("functions.jl")
include("../functions.jl")

# %% load the trajectories
cd(@__DIR__)
scale(t, s=0.13319672) = t .* s # px to μm
times(t) = 0:length(t)-1 ./ 10
df = load_trajectories()
c = mean(df.t[end], dims=1)
df

# %%
# calculate the probability distribution of the displacements
# center(t) = mean(t, dims=1)
radius(t, c=c) = [norm(t[j,:]' .- c) for j in 1:size(t,1)]
σ(t) = std(radius(t))

# a + b * (x - c)^2
model(x, p) = p[1] .+ p[2] .* (x .- p[3]).^2
model_string(p) = format(p[2], precision=2)*" ⋅ r²"


# %% plot the potential as a radial distribution
f = Figure()
a = Axis(f[1, 1], ylabel="pdf(x)")
ylims!(a, 0, nothing)
b = Axis(f[2, 1], xlabel="r in μm", ylabel= "V(x) in kT")
linkxaxes!(a, b)
hidexdecorations!(a, grid=false)
for i in eachrow(filter(x-> x.ot>0, df))
	r = scale(radius(i.t))
	k = dist(r, cutoff=0.1)
	x = k.x
	lines!(a, x, k.y, color=i.ot, colorrange=extrema(df.ot), label=format(i.ot, precision=2))
	lines!(b, x, potential(k.y), color=i.ot, colorrange=extrema(df.ot))

	m = curve_fit(model, x, potential(k.y), [0., 1., 0.]).param
	lines!(b, x, model(x, m), color=i.ot, colorrange=extrema(df.ot), linestyle=:dash, label=model_string(m))
end
axislegend(a, "Trap", position=:lt)
axislegend(b, "Fit", position=:lt)
save("../figures/01_03_1_potential.pdf", f)
f

# %% plot the measured spring constants
f = Figure()
a = Axis(f[1, 1], xlabel="ot", ylabel="k")

scatter!(df.ot, 1 ./σ.(df.t), label="σ")

k_fit = zeros(size(df, 1))
k_error = zeros(size(df, 1))
for i in 1:size(df, 1)
	k = dist(radius(scale(df.t[i])), cutoff=0.01)
	x = k.x .- mean(k.x)
	m = curve_fit(model, x, potential(k.y), [0., 1., 0.])
	k_fit[i] = m.param[2]
	k_error[i] = sqrt(m.jacobian[2,2])
end
# scatter!(df.ot, k_fit, label="fit")
errorbars!(df.ot, k_fit, k_error, label="fit")

axislegend(position=:lt)
f