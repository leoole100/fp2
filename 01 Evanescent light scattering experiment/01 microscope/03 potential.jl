using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit, stderror
import LsqFit
using LinearAlgebra: norm, normalize
using CairoMakie
import CSV, Images
using StatsBase: fit, Histogram
using KernelDensity: kde
using Glob: glob
using MeanSquaredDisplacement: imsd
using Interpolations: interpolate, linear_interpolation
using Format: format
using Measurements: measurement, value, uncertainty
using JLD2

include("functions.jl")
include("../functions.jl")

# load the trajectories
cd(@__DIR__)
scale(t, s=0.13319672) = t .* s # px to μm
times(t) = 0:length(t)-1 ./ 10
# df = load_trajectories()
df = load("../data/TLM/02.jld2")["df"]

c = mean(df.t[end], dims=1)

# %%
# calculate the probability distribution of the displacements
# center(t) = mean(t, dims=1)
radius(t) = [norm(t[j,:]') for j in 1:size(t,1)]
σ(t) = std(radius(t))

# a + b * (x - c)^2
model(x, p) = p[1] .+ p[2] .* (x .- p[3]).^2
model_string(p) = format(p[2], precision=2)*" ⋅ r²"


# %% bivariate
f = Figure(size=fullsize)
axs = [Axis(f[j, i]) for i in 1:2, j in 1:2]
for i in 1:size(df, 1)
	a = hcat(axs...)[i]
	a.title = format(df.ot[i], precision=2)
	# linkaxes!(a, axs[1,1])
	a.aspect = DataAspect()
	t = scale(df.t[i].-c)
	# t = scale(df.t[i].-mean(df.t[i], dims=1))
	k = kde(t)
	h = heatmap!(a, k.x, k.y, k.density, colormap=:binary)
	# h = heatmap!(a, k.x, k.y, potential(k.density))
	# draw a circle
	arc!(a, Point2f(0), 2.5, -π, π, color=:red, alpha=.5, linestyle=:dash)
end
linkaxes!(axs[2,2], axs[1,2])
save("../figures/01_03_1_bivariate.pdf", f)
f

# %% radial distribution
f = Figure()
aV = Axis(f[1, 1], xlabel="r in μm", ylabel="V(x) in kT")
aP = Axis(f[2, 1], xlabel="r in μm", ylabel="pdf(x)")
linkxaxes!(aV, aP)
hidexdecorations!(aV, grid=false)
for i in eachrow(filter(d -> d.ot>0, sort(df, :ot)))
	# r = scale(radius(i.t .-mean(i.t, dims=1)))
	r = scale(radius(i.t .-c))
	k = dist(r, cutoff=0.01, boundary=(0, maximum(r)))
	lines!(aP, k.x, k.y, color=i.ot, colorrange=extrema(df.ot), label=format(i.ot, precision=2))
	l = lines!(aV, k.x, potential(k.y), color=i.ot, colorrange=extrema(df.ot))
	# stephist!(aP, r, normalization=:pdf)
end
ylims!(aV, nothing, 5)
ylims!(aP, 0, nothing)
axislegend(position=:rt)
save("../figures/01_03_2_radial.pdf", f)
f

# %% group by the coordinates

f = Figure(size=halfsize)
a = Axis(f[1, 1], xlabel="x,y in μm", ylabel="V(x,y) in kT")
df.Vkx = fill(measurement(0., 0.), size(df, 1))
df.Vky = fill(measurement(0., 0.), size(df, 1))
i = df[2,:]
t = scale(i.t .-mean(i.t, dims=1))
for j in 1:2
	k = dist(t[:,j], cutoff=0.05)
	mdl = curve_fit(model, k.x, potential(k.y), [0., 1., 0.])
	offset = mdl.param[3]
	k = dist(t[:,j], cutoff=0)
	lines!(a, k.x .-offset, potential(k.y), color=i.ot, colorrange=extrema(df.ot), label=format(i.ot, precision=2)) 	
end
for m in 2:size(df, 1)
	i = df[m,:]
	t = scale(i.t .-mean(i.t, dims=1))
	fp = zeros(2)
	for j in 1:2
		k = dist(t[:,j], cutoff=0.05)
		mdl = curve_fit(model, k.x, potential(k.y), [0., 1., 0.])
		offset = mdl.param[3]
		lines!(a, k.x .-offset, potential(k.y), color=i.ot, colorrange=extrema(df.ot), label=format(i.ot, precision=2))
		lines!(a, 
			k.x .-offset, model(k.x, mdl.param), 
			color=i.ot, colorrange=extrema(df.ot), 
			linestyle=:dash,
			# label=model_string(mdl.param)
		)

		k = measurement(mdl.param[2], stderror(mdl)[2])
		println(k)
		if j == 1
			df[m ,:Vkx] = k
		else
			df[m ,:Vky] = k
		end
	end
end
i = df[1,:]
t = scale(i.t .-mean(i.t, dims=1))
for j in 1:2
	k = dist(t[:,j], cutoff=0.01)
	lines!(a, k.x, potential(k.y), color=i.ot, colorrange=extrema(df.ot), label=format(i.ot, precision=2))
end
xlims!(a, -5,5)
ylims!(a, 0, 5)
Colorbar(f[1, 2], limits=extrema(df.ot), label="Trap Strength")
colgap!(f.layout, 5)
# axislegend("Trap Stiffness", position=:lb, unique=true)
# Legend(f[1, 2], a, framevisible=false, unique=true)
save("../figures/01_03_3_axis.pdf", f)
f


# %% plot the measured spring constants
f = Figure(size=halfsize)
a = Axis(f[1, 1], xlabel="Trap Stiffness", ylabel="k in N/m")

function er(a, x, y, label, color=Nothing)
	errorbars!(a, x, value.(y), uncertainty.(y))
	scatter!(a, x, value.(y), label=label, markersize=7)
end

T = measurement(22.6, 0.1) + 273.15
kB = 1.38064852e-23
kT = kB * T

# plot fit results
p = er(a, df.ot[2:end], kT.*df.Vkx[2:end]*1e12, "V(x)")
er(a, df.ot[2:end], kT.*df.Vky[2:end]*1e12, "V(y)")

# plot the msd_inf
er(a, df.ot[3:end], kT./df.msd_inf[3:end] * 1e12, "MSD(∞)")

axislegend(position=:lt, framevisible=false)
# Legend(f[1, 2], a, framevisible = false)

save("../figures/01_03_4_spring_constants.pdf", f)
f