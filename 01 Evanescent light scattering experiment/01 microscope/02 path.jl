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
using Glob: glob
using MeanSquaredDisplacement: imsd
# %%

path = "../data/TLM/"

# open the trajectory
cd(@__DIR__)
# img = Images.load(path * ".png");

function load_trajectory!(df, path)
	path = join(split(path, ".")[1:end-1], ".")
	f = DataFrame(CSV.File(path * ".csv"))
	for i in 1:ncol(f)÷2
		push!(df[!, :t], hcat(f[:,2*i-1], f[:,2*i]))
		n = split(path, "/")[end]
		n = join(split(n, ".")[1:end-1], ".")
		n = split(n, " ")[3]
		ot = parse(Float16, n[3:end])
		push!(df[!, :n], n)
		push!(df[!, :ot], ot)
	end
end

df = DataFrame(t=[], n=[], ot=[])
paths = glob(path * "*.trajectory.csv")
for p in paths
	load_trajectory!(df, p)
end
sort(df, :ot)
df

# plot the trajectories on the image
f = Figure()
image(f[1,1], img', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)
for p in reverse(eachrow(df))
	lines!(p.t, label=p.n, alpha=0.5)
	# points = [Point2f(p.trajectory[j,:]...) for j in p.frames]
	# datashader!(points)
end
axislegend()
f

# %%
# calculate the mean squared displacement
scale(t, s=0.13319672e-06) = t .* s
times(t) = 0:length(t)-1 ./ 10
function remove_drift(t)
	t = copy(t)
	t = t .- t[1,:]'
	s = t[end,:] ./ (size(t,1)-1)
	for i in 0:(size(t,1)-1)
		t[i+1,:] = t[i+1,:] .- i .* s
	end
	return t
end

# with p[1] Diffunsionskoeffizient, in dimension 2
model(t, p) = p[4] ./(p[2].*t.^p[1].+ p[3])
fit_model(t, m, p0=[-1, 2e13, 1.1e12,1]) = curve_fit(model, t, m, p0)

# %%
f = Figure();
a = Axis(f[1, 1], 
	xlabel="time in s", ylabel="mean squared displacement in m^2",
	xscale=log10, 
	yscale=log10,
)
for i in eachrow(filter(d -> d.ot<2, df))
	t=scale(i.t)
	# t = remove_drift(t)
	m = imsd(t)
	m = mean(m, dims=2)[:,1] # mean over x and y
	m = abs.(m)
	s = lines!(times(m), m, label=string(i.ot), color=i.ot, colorrange=(0, 1.1))
	# mdl = fit_model(times(m), m)
	# lines!(times(m), t->model(t, mdl.param), linestyle=:dash)

end
# lines!(1:1000 ./10, t->model(t, [-1, 2e13, 1.1e12,1]), color=:black, label="model")
axislegend(a, position=:lt);
ylims!(a, 5e-14, 1e-10)
xlims!(a, 1, nothing)
save("../figures/01_msd.pdf", f)
f

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