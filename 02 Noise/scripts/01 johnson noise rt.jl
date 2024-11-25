#=
Script to analyse the amplified Johnson noise of a Resistor.
amplified with preamp G1 and amp G2, with a assumed constant noise on G1.
=#

using DataFrames: DataFrame, groupby, eachcol, select, combine
using CSV: CSV
using CairoMakie
using Measurements: measurement, value, uncertainty
using LsqFit: curve_fit, stderror

# %%
# Load the data
cd(@__DIR__)
df = DataFrame(CSV.File("../data/01 johnson noise rt.csv"))

# create the measurement (the scaled V^2)
df.V2 = measurement.(df.Vsq, df.VsqU)
# remove the gain
df.V2 = df.V2 .* 10 # 10 is a additional gain
df.V2 = df.V2 ./ (df.G1 .* df.G2).^2

# estimate the noise for each preamp configuration independently
function estimate_noise(df)
	return minimum(df.V2)	
end
# remove the noise from the Measurements
grouped = groupby(df, :G1)
for g in grouped
	noise = estimate_noise(g)
	g.noise .= noise
	g.VJ2 .= g.V2 .- noise
end

# fit a line to estimate k_b=1/(4TΔf) VJ²/R
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [1.0, 0.0]
T = measurement(25.0, 0.1) + 273.15
grouped = groupby(df, :Δf)
n = length(grouped)
fit_params = DataFrame(
	Δf = unique(df.Δf),
	kb = fill(measurement(0.0, 0.0), n),
	fit = fill(p0, n)
)
for (g, f) in zip(grouped, eachrow(fit_params))
	p = curve_fit(mdl, value.(g.R), value.(g.VJ2), p0)
	f.kb = 1 / (4T * f.Δf) * measurement(p.param[2], stderror(p)[2])
	f.fit = p.param
end

# plot the noise over the resistance
f = Figure()
a = Axis(f[1, 1];
	xlabel="Resistance [Ω]", ylabel="V² Junction [V²]",
	xscale=log10, yscale=log10
)
# marker = [MarkerElement(marker=m) for m in [:x, :o]]
marker = [:x, :o]
colors = [PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)]
for d in eachrow(df)
	scatter!(d.R, value(d.V2), 
		marker=marker[1],
		color=d.Δf, colorrange=extrema(df.Δf)
	)
	scatter!(d.R, value(d.VJ2), 
		marker=marker[2],
		color=d.Δf, colorrange=extrema(df.Δf)
	)
end
for f in eachrow(fit_params)
	x = logrange(1, 1e3, 100)
	lines!(x, mdl(x, f.fit), color=f.Δf, colorrange=extrema(df.Δf))
end

axislegend(a,
	[MarkerElement(marker=m) for m in marker],
	["V²", "VJ²"],
	position=:rb
)

axislegend(a,
	[PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)],
	string.(unique(df.Δf)),
	"Δf",
	position=:lt
)
f
