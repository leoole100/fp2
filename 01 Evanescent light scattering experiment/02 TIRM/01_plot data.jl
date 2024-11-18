using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit
using LinearAlgebra: norm, normalize
using CairoMakie
using StatsBase: fit, Histogram
using KernelDensity: kde
import CSV, Images
using Glob: glob
using Format: format
using JLD2 
include("../functions.jl")

cd(@__DIR__)

# read all CSV files in the directory
path = "../data/TIRM/"
paths = glob(path * "*/*.dat")
paths = filter(p -> occursin("failed", lowercase(p)) == false, paths)
Ia = []
pa = []
for p in paths
	I = DataFrame(CSV.File(p, header=["t", "I"])).I
	push!(Ia, I)
  push!(pa, split(p, "/")[end-1])
end
df = DataFrame(I=Ia, p=pa)
df.l = map(p -> parse(Float64, split(p, " ")[2][1:end-2]), df[:,:p])
df.ot = map(p -> parse(Float64, split(p, " ")[3][3:end]), df[:,:p])
sort!(df, :ot)

save("../data/TIRM/data.jld2", "df", df)

df
#%%

# plot the I density
f = Figure()
a = Axis(f[1, 1], xlabel="intensity")
for i in eachrow(df)
	# I = remove_drift(i.I)[:]
	I = i.I
	k = kde(I, boundary=(minimum(I), maximum(I)))
	lines!(a, k.x, k.density, label=i.p, 
		color=i.ot, colorrange=extrema(df.ot),
	)
end
ylims!(0,16)
axislegend()
save("../figures/02_010_Intensity_Distribution.pdf", f)
f

#%%

# plot the ΔI density
f = Figure()
a = Axis(f[1, 1], xlabel="ΔI")
for i in eachrow(df)
	I = i.I
	I = diff(I)
	k = kde(I, boundary=(minimum(I), maximum(I)))
	lines!(a, k.x, k.density, label=i.p, 
		color=i.ot, colorrange=extrema(df.ot),
	)
end
xlims!(-.1, .1)

x = -.1:.001:.1
y = 1 .* exp.(-10 .* x)
y = 160 .* exp.(- x.^2 ./ .018^2)
lines!(x,y)

# save("../figures/02_020_IntensityDiff_Distribution.pdf", f)
f


# %%
background(I) = -.1
estimate_I0(I) = 2*maximum(I)-background(I)
calc_z(I, I0=estimate_I0(I), beta=1, b=background(I)) = log.(I0 ./ (I.-b)) ./ beta
potential(p) = -log.(p) .+ log(maximum(p))
function dist(data; cutoff=0)
	k = kde(data, boundary=(minimum(data), maximum(data)))
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

model(x, p) = p[1] .* exp.(-x .* p[2]) .+ p[3] .* x .- p[4]
fit_model(x, y) = curve_fit(model, x, y, [10, 1.0, 1.0, 1.0]).param

# model_ot(x, p) = p[1] .* exp.(-x .* p[2]) .+ p[3] .* x .- p[4] .+ p[5] .* (x.-p[6]).^2
# fit_model_ot(x, y) = curve_fit(model, x, y, [10, 1.0, 1.0, 1.0, .1, 1.5]).param

#%%

f = Figure()
# aV = Axis(f[1, 1], ylabel="V in kT")
# az = Axis(f[2, 1], ylabel="pdf", xlabel="z")
aI = Axis(f[2, 2], xlabel="intensity")
# linkxaxes!(az, aV)
# linkyaxes!(aI, az)
# ylims!(az, 0, nothing)
# hidexdecorations!(aV, grid=false)
# hideydecorations!(aI, grid=false)
# for i in eachrow(filter(d -> d.ot==1.5, sort(df, :l)))
for i in eachrow(filter(d -> d.ot<1, sort(df, :l)))
	pI = dist(i.I; cutoff=0.1)
	lines!(aI, pI.x, pI.y, label=i.l, color=i.ot, colorrange=extrema)
	z = calc_z(i.I)
	pZ = dist(z; cutoff=0.1)
	# lines!(az, pZ.x, pZ.y, label=i.l)
	v = potential(pZ.y)
	# l = lines!(aV, pZ.x, v, label=i.l)
	# lines!(aV, pZ.x, model(pZ.x, fit_model(pZ.x, v)), color=l.color, linestyle=:dash)
end
# lines!(aV, 0:.1:8, x -> model(x, [10, 1.0, 1.0, 1.0]), color=:black)
# lines!(aV, .5:.1:3, x -> model_ot(x, [10, 1.0, 1.0, 1.0, .1, 1.5]), color=:black)
# axislegend(aV, position=:rb)
# save("../figures/02_potential.png", f)
f