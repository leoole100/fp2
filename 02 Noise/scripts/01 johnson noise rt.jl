#=
Script to analyse the amplified Johnson noise of a Resistor to estimate (kb T).
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
df = DataFrame(CSV.File("../data/01 johnson noise rt.csv"))

# df.V = measurement.(df.Vsq, df.VsqU) 
df.V = measurement.(df.Vsq, load("../data/gen/04 meter uncertainty.jld2")["std"]) 
df.V = df.V .* 10 ./ (df.G1 .* 100 .* df.G2).^2 # scale measurements to volts²
df.Δf = 1e3*df.Δf # convert to Hz
df.S = df.V ./ df.Δf
select!(df, Not([:VsqU, :f1, :f2, :G1, :G2, :Vsq])) # remove unused columns

sort(df, [:Δf, :R])


# %%
# fit a line for each group of Δf
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [1.0, 0.0]

# create new columns for the fit parameters
fits = DataFrame(
	p = [
		curve_fit(mdl, value.(g.R), value.(g.S), p0)
		for g in groupby(df, :Δf)
	],
	Δf = unique(df.Δf)
)
fits.p = [measurement.(p.param, stderror(p)) for p in fits.p]
fits.kT = [d.p[2] / 4 for d in eachrow(fits)]
T = measurement(22.0, 3) + 273.15
fits.k = [kT / T for kT in fits.kT]
fits.S0 = [p[1] for p in fits.p]
fits


# %%
f = Figure(size=fullsize)
s = 1e15
a = Axis(f[1, 1]; 
	ylabel="S in f V²/Hz",
	# yscale=log10,
	# xscale=log10,
)
for d in groupby(df, :Δf)
	scatter!(d.R, value.(d.S).*s, color=d.Δf, colorrange=extrema(df.Δf))
	errorbars!(d.R, value.(d.S).*s, uncertainty.(d.S).*s, color=d.Δf)
end
for f in eachrow(fits)
	x = range(1, 10e3, 100)
	lines!(x, mdl(x, value.(f.p)).*s, color=f.Δf, colorrange=extrema(df.Δf))
end
axislegend(a,
	[PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)],
	format.(round.(unique(df.Δf)/1e3)),
	"Δf / kHz",
	orientation=:horizontal,
	position=:lt
)

ar = Axis(f[2,1];
	ylabel="Residuals",
	xlabel="Resistance in Ω", 
	# xscale=log10,
)
for (d, fit) in zip(groupby(df, :Δf), eachrow(fits))
	y = s.*(value.(d.S) - mdl(value.(d.R), value.(fit.p)))
	errorbars!(d.R, y, s.*uncertainty.(d.S); color=:gray)
	scatterlines!(d.R, y, color=fit.Δf, colorrange=extrema(df.Δf))
end
linkxaxes!(a, ar)
hidexdecorations!(a; grid=false)
rowgap!(f.layout, 5)
rowsize!(f.layout, 1, Relative(3/4))
save("../figures/01 johnson noise.pdf", f)
f


# %%
save("../data/gen/01 Johnson noise RT.jld2", Dict(
	"S0"=>mean(fits.S0)
))

# %%
# collect results as a DataFrame
results = [fits.kT, fits.k, fits.S0]
DataFrame(
	:parameter => ["kT", "k", "S0"],
	:mean => [mean(value.(r)) for r in results],
	:accuracy => [std(values.(r)) for r in results],
	:precision => [uncertainty(mean(r)) for r in results],
	:measurement => [measurement(mean(value.(r)), std(value.(r))) for r in results],
)

# mean for k
k = measurement(
	mean(value.(fits.k)),
	sqrt(mean(uncertainty.(fits.k))^2 + std(value.(fits.k))^2)
)