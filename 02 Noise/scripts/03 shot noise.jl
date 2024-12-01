#=
Script to analyse the amplified Shot noise to estimate e.
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

function estimate_e(df)
	# df.V = measurement.(df.Vsq, df.VsqU) 
	df.V = measurement.(df.Vsq, load("../data/gen/04 meter uncertainty.jld2")["std"]) 
	df.V = df.V .* 10 ./ (100 .* df.G2 .*1000).^2 # scale measurements to volts²
	df.V *= 1/(10000)^2 	# convert to current, with a 10kΩ resistor
	df.Δf = 1e3*df.Δf # convert to Hz
	df.I *= 1e-6 # to A
	df.S = df.V ./ df.Δf
	sort!(df, [:I, :Δf])

	return vcat([diff(d.S) ./ diff(d.I) / 2 for d in groupby(df, :Δf)] ...)
end

f = Figure(size=halfsize)
s = 1e19
a = Axis(f[1,1],
	xlabel="measured e in 10^$(format(log10(1/s))) C"
)

dfs = [
	DataFrame(CSV.File("../data/04 photocurrent.csv")),
	DataFrame(CSV.File("../data/05 transimpedance amplifier.csv"))
]

for (df, l) in zip(dfs, ["I", "TIA"])
	e = estimate_e(df)

	# remove outliers
	e = e[.1e-19 .<e.<10e-19]

	hist!(value.(e).*s,
		label=l,
	)
	
	println(
		"$(l):\t",
		measurement(
			mean(value.(e)),
			std(value.(e))
		)
	)
end

vlines!([1.602e-19].*s, color=:black)
ylims!(low=0)
# hideydecorations!(a)
# axislegend(a, position=:rt)
save("../figures/03 shot noise.pdf", f)
f