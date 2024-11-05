using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit
using LinearAlgebra: norm, normalize
using WGLMakie
using StatsBase: fit, Histogram
using KernelDensity: kde
import CSV, Images

cd(@__DIR__)

# %%
# generate some random data
df = DataFrame(I=[
	exp.(randn(10000000)),
	exp.(-exp.(randn(10000000)))
], l = ["A", "B"])

# read all CSV files in the directory
# path = "../data/02/"
# files = filter(x -> occursin(".csv", x), readdir(path))
# df = DataFrame(I=[])
# for file in files
# 	I = DataFrame(CSV.File(path * file))[:, :I]
# 	push!(df[:, :I], I)
#   push!(df[:, :l], split(file, "/")[1])
# end

# plot the data
# f = Figure()
# a = Axis(f[1, 1], xlabel="time", ylabel="intensity")
# for i in eachrow(df)
# 	lines!(i.I, label=i.l)
# end
# axislegend()
# f

# %%
background(I) = 0
estimate_I0(I) = maximum(I)-background(I)
calc_z(I, I0=estimate_I0(I), beta=1, b=background(I)) = log.(I0 ./ (I.-b)) ./ beta
potential(p) = -log.(p) .+ log(maximum(p))
function dist(data; cutoff=0)
	k = kde(data, boundary=(minimum(data), maximum(data)))
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

model(x, p) = p[1] .* exp.(-x .* p[2]) .+ p[3] .* x .- p[4]
fit_model(x, y) = curve_fit(model, x, y, [10, 1.0, 1.0, 1.0])

f = Figure()
aV = Axis(f[1, 1], ylabel="V in kT")
# az = Axis(f[2, 1], ylabel="pdf", xlabel="z")
# aI = Axis(f[2, 2], xlabel="intensity")
# linkxaxes!(az, aV)
# linkyaxes!(aI, az)
# ylims!(az, 0, nothing)
# hidexdecorations!(aV, grid=false)
# hideydecorations!(aI, grid=false)
for i in eachrow(df)
	# pI = dist(i.I; cutoff=0.1)
	# lines!(aI, pI.x, pI.y, label=i.l)
	z = calc_z(i.I)
	pZ = dist(z; cutoff=0.01)
	# lines!(az, pZ.x, pZ.y, label=i.l)
	v = potential(pZ.y)
	l = lines!(aV, pZ.x, v, label=i.l)
	fit = fit_model(pZ.x, v)
	lines!(aV, pZ.x, model(pZ.x, fit.param), color=l.color, linestyle=:dash)
end
# lines!(aV, 0:.1:8, x -> model(x, [10, 1.0, 1.0, 1.0]), color=:black)
axislegend(aV)
f