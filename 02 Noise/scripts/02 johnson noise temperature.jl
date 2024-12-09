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
using Format
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
df.S = df.S .- load("../data/gen/05 Johnson noise RT.jld2")["S0"]

sort!(df, [:T, :Δf, :R])
select(df, [:T, :Δf, :R, :S])

# %%
# fit lines
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [0., 1.0]
sort!(df, [:R, :Δf], rev=[false, false])
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
a = Axis(f[1:2, 1]; 
	ylabel="S in (μV)²/Hz",
	xlabel="T in K",
)
linestyles = Dict(
	unique(df.R) .=> [:solid, :dash, :dot]
)
for (d, f) in zip(groups, eachrow(fits))
	s = 1e12
	scatter!(value.(d.T), value.(d.S).*s, alpha=0.7,  
	color=f.Δf, colorrange=extrema(df.Δf)
	)
	x = range(-30, 300, 100)
	lines!(x, mdl(x, f.p).*s, alpha=0.7, 
		label="$(format(f.R, autoscale=:metric))Ω, \t $(format(f.Δf/1000, precision=0)) kHz",
		color=f.Δf, colorrange=extrema(df.Δf),
		linestyle=linestyles[f.R]
	)
end
ylims!(low=0)
Legend(
	f[1,2],
	[
		PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)
	],
	format.(unique(df.Δf)./1000, precision=0),
	"Δf / kHz",
	framevisible=false
)
Legend(
	f[2,2],
	[
		LineElement(linestyle=linestyles[c]) for c in unique(df.R)
	],
	format.(unique(df.R), autoscale=:metric),
	"R / Ω",
	framevisible=false
)

colgap!(f.layout, 5)
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