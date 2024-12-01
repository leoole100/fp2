#=
Script to analyse the amplified Johnson noise of a Resistor to estimate T.
amplified with preamp G1 and amp G2, with a assumed constant noise on G1.
=#

using DataFrames
using CSV: CSV
using CairoMakie
using Measurements: measurement, value, uncertainty
using LsqFit: curve_fit, stderror
using Format: format
using StatsBase
using JLD2
include("functions.jl")

# %%
# Load the data
cd(@__DIR__)
df = DataFrame(CSV.File("../data/06 temperatures.csv"))

df.V = measurement.(df.Vsq, load("../data/gen/04 meter uncertainty.jld2")["std"]) 
df.T = measurement.(df.T, 3)
df.V = df.V .* 10 ./ (600 .* df.G2).^2 # scale measurements to volts²
df.Δf = 1e3*df.Δf # convert to Hz
df.S = df.V ./ df.Δf
df.S = df.S .- load("../data/gen/01 Johnson noise RT.jld2")["S0"]

sort!(df, [:T, :Δf, :R])
select(df, [:T, :Δf, :R, :S])

# %%
# fit lines
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [0., 1.0]
sort!(df, [:R, :Δf], rev=[true, false])
groups = groupby(df, [:R, :Δf], sort=false)
fits = DataFrame(
	Δf = getindex.(keys(groups), 2),
	R = getindex.(keys(groups), 1),
	p = [
		curve_fit(mdl, value.(g.T), value.(g.S), p0)
		for g in groups
	]
)
fits.p = [p.param for p in fits.p]

f = Figure(size=fullsize)
a = Axis(f[1, 1]; 
	ylabel="S in nV²/Hz",
	xlabel="T in K",
)
for (d, f) in zip(groups, eachrow(fits))
	s = 1e9
	scatter!(value.(d.T), value.(d.S).*s)
	x = range(-30, 300, 100)
	lines!(x, mdl(x, f.p).*s, label="$(f.R)Ω, $(f.Δf)Hz")
end
ylims!(low=0)
Legend(f[1,2], a, framevisible=false)
save("../figures/02 temperature.pdf", f)
f

# %%
# calculate where the lines intersect
# T = (S1 - S2) / (m2 - m1)

T = [
	(f1.p[1] - f2.p[1]) / (f2.p[2] - f1.p[2])
	for f1 in eachrow(fits)
	for f2 in eachrow(fits)
]
T = T[.!isnan.(T)]
T =  T[-100 .< T .< 100]	# remove outliers

f = Figure(size=halfsize)
a = Axis(f[1,1],
	xlabel="T in K"
)
# density!(T)
hist!(T)
vlines!([0], color=:black)
# hideydecorations!(a)
ylims!(low=0)
save("../figures/02 temperature distribution.pdf", f)
f

#%%
# calculate the mean and uncertainty
T_mean = measurement(mean(T), std(T))