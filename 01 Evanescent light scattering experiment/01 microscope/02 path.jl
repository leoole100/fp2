# find particles in images and track them

# using DataFrames, Statistics, LsqFit, LinearAlgebra
using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit
using LinearAlgebra: norm, normalize
using WGLMakie
import CSV, Images
using StatsBase: fit, Histogram
using KernelDensity: kde

# %%

path = "../data/test/particles"

# open the trajectory
cd(@__DIR__)
img = Images.load(path * ".mean.png");
f = DataFrame(CSV.File(path * ".trajectory.csv"))

df = DataFrame(t=[], n=[])

for i in 1:ncol(f)÷2
	push!(df[!, :t], hcat(f[:,2*i-1], f[:,2*i]))
	push!(df[!, :n], split(path, "/")[end] * " $i")
end
df

# plot the trajectories on the image
f = Figure()
image(f[1,1], img', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)
for p in eachrow(df)
	scatter!(p.t, label=p.n)
	# points = [Point2f(p.trajectory[j,:]...) for j in p.frames]
	# datashader!(points)
end
axislegend()
f

# %%
# calculate the mean squared displacement

reference(t) = t[1,:]
msd(t, r=reference(t)) = sum((t .- r').^2, dims=2)[:,1]
msd(t) = msd(t, reference(t))
msd(df[1,:t])

linear_model(x, p) = clamp.(p[1] .* x, 0, p[2]);
exp_model(x, p) = p[1] .- p[1] .* exp.(-p[2] .* x ./ p[1]);
fit_linear_model(t, times=1:size(t,1), p0=[1.0, 1.0]) = curve_fit(linear_model, times, msd(t), p0)
fit_exp_model(t, times=1:size(t,1), p0=[1000., 100.]) = curve_fit(exp_model, times, msd(t), p0)

f = Figure();
a = Axis(f[1, 1], xlabel="time", ylabel="mean squared displacement");
for i in eachrow(df)
	scatter!(msd(i.t), label=i.n)
	for (m, l) in zip([fit_linear_model, fit_exp_model], ["linear", "exp"])
		p = m(i.t)
		lines!(0:100, x->linear_model(x, p.param))
	end
end
# try models
lines!(0:100, x->linear_model(x, [0, 100., 1500.]), label="linear");
lines!(0:100, x->exp_model(x, [1000., 100.]), label="exp");
axislegend(a);
f

# %%
# calculate the probability distribution of the displacements
center(t) = mean(t, dims=1)
radius(t, c=center(t)) = [norm(t[j,:]' .- c) for j in 1:size(t,1)]
σ(t) = std(radius(t))
σ(df[1,:t]) # = dV/dx^2 / kT

potential(p) = -log.(p)
histogram(data) = normalize(fit(Histogram, data), mode=:probability)

f = Figure()
a = Axis(f[1, 1], ylabel="probability")
b = Axis(f[2, 1], xlabel="distance from center", ylabel="Potential in kT")
linkxaxes!(a, b)
hidexdecorations!(a, ticks = false)
for i in eachrow(df)
	r = radius(i.t)
	h = histogram(r)
	k = kde(r, boundary=(minimum(r), maximum(r)))

	l = stairs!(a, h.edges[1][1:end-1], h.weights, linestyle=:dash)
	scatter!(b, h.edges[1][1:end-1], potential(h.weights), color=l.color)

	lines!(a, k.x, k.density, color=l.color)
	lines!(b, k.x, potential(k.density), color=l.color, label=i.n)
end
axislegend()
save("../figures/01_potential.png", f)
f