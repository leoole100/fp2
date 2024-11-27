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

groups = groupby(df, [:Δf, :R])
fits = DataFrame(
	Δf = getindex.(keys(groups), 1),
	R = getindex.(keys(groups), 2),
	p = [
		curve_fit(mdl, value.(g.T), value.(g.S), p0)
		for g in groups
	]
)
fits.p = [p.param for p in fits.p]


# %%
f = Figure()
s = 1e15
a = Axis(f[1, 1]; 
	ylabel="S in 10^$(format(log10(s))) V²/Hz",
	xlabel="T in K",
)
for (d, f) in zip(groups, eachrow(fits))
	scatter!(value.(d.T), value.(d.S).*s)
	x = range(-30, 300, 100)
	lines!(x, mdl(x, f.p).*s, label="$(f.R)Ω, $(f.Δf)Hz")
	# errorbars!(value.(d.T), value.(d.S).*s, uncertainty.(d.S).*s, uncertainty.(d.T), color=d.Δf)
end
ylims!(low=0)
axislegend(position=:lt)
save("../figures/02 temperature.pdf", f)
f