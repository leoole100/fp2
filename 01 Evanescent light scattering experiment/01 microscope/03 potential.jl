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

# %% load the trajectories
path = "../data/TLM/"
cd(@__DIR__)

function load_trajectory!(df, path)
	path = join(split(path, ".")[1:end-1], ".")
	f = DataFrame(CSV.File(path * ".csv"))
	for i in 1:ncol(f)÷2
		push!(df[!, :t], hcat(f[:,2*i-1], f[:,2*i]))
		n = split(path, "/")[end]
		n = join(split(n, ".")[1:end-1], ".")
		n = split(n, " ")[3]
		ot = parse(Float64, n[3:end])
		push!(df[!, :n], n)
		push!(df[!, :ot], ot)
	end
end

df = DataFrame(t=[], n=[], ot=[])
paths = glob(path * "*.trajectory.csv")
for p in paths
	load_trajectory!(df, p)
end

scale(t, s=0.13319672) = t .* s # px to μm
times(t) = 0:length(t)-1 ./ 10

sort!(df, :ot)
df

# %%
# calculate the probability distribution of the displacements
center(t) = mean(t, dims=1)
radius(t, c=center(t)) = [norm(t[j,:]' .- c) for j in 1:size(t,1)]
σ(t) = std(radius(t))
σ(df[1,:t]) # = dV/dx^2 / kT

potential(p) = -log.(p) .+ log(maximum(p))
function histogram(data)
	h = normalize(fit(Histogram, data), mode=:pdf)
	x = h.edges[1][1:end-1] .+ diff(h.edges[1]) ./ 2
	return (x=x[1:end-1], density=h.weights[1:end-1], bins=h.edges[1][1:end-1])
end

f = Figure()
a = Axis(f[1, 1], ylabel="pdf(x)")
ylims!(a, 0, nothing)
b = Axis(f[2, 1], xlabel="distance from center x", ylabel= "V(x) in kT")
linkxaxes!(a, b)
hidexdecorations!(a, grid=false)
for i in eachrow(df)
	r = radius(i.t)
	h = histogram(r)
	k = kde(r, boundary=(minimum(r), maximum(r)))

	l = stairs!(a, h.x, h.density, step=:center)
	scatter!(b, h.x, potential(h.density), color=l.color)

	lines!(a, k.x, k.density, color=l.color)
	lines!(b, k.x, potential(k.density), color=l.color, label=i.n)
end
axislegend(position=:lt)
save("../figures/01_potential.png", f)
f