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

#%% plot the trajectories
f = Figure()
a = Axis(f[1,1], xlabel="x in μm", ylabel="y in μm")
for p in eachrow(df)
	t = scale(p.t)
	t = t .- mean(t, dims=1)
	lines!(t, label=p.n, alpha=0.5)
	# points = [Point2f(p.trajectory[j,:]...) for j in p.frames]
	# datashader!(points)
end
axislegend(position=:lt)
save("../figures/01_trajectories.pdf", f)
f

# %% calculate the mean squared displacement
function remove_drift(t)
	t = copy(t)
	t = t .- t[1,:]'
	s = t[end,:] ./ (size(t,1)-1)
	for i in 0:(size(t,1)-1)
		t[i+1,:] = t[i+1,:] .- i .* s
	end
	return t
end

function msd(t)
	t = scale(t)
	t = remove_drift(t)
	m = imsd(t)
	m = mean(m, dims=2)[:,1] # mean over x and y	
	m = abs.(m)
	return m
end

interp_msd = linear_interpolation(
	df.ot,
	msd.(df.t)
)

#%%

f = Figure()
a = Axis(f[1, 1], 
	xlabel="time in s", 
	ylabel=rich("msd in μm", superscript("2")),
	xscale=log10, 
	yscale=log10,
)

# interpolated values
# for ot in range(0.7, 1.0, 10)
# 	m = interp_msd(ot)
# 	s = lines!(times(m), m, color=ot, colorrange=extrema(df.ot))
# end

# measurements
for i in eachrow(filter(d -> d.ot<2, df))
	m = msd(i.t)
	s = scatter!(times(m), m, label=format(i.ot, precision=2), color=i.ot, colorrange=extrema(df.ot))
end

# Colorbar(f[1, 2], limits=extrema(df.ot), label="Optical trap strength")
axislegend(a, "Trap Strength", position=:lt)
ylims!(a, 1e-2, 1e2)
save("../figures/01_msd.pdf", f)
f
