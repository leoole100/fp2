#=
Script to analyse the amplified Johnson noise of a Resistor to estimate (kb T).
amplified with preamp G1 and amp G2, with a assumed constant noise on G1.
=#

using DataFrames: DataFrame
using CSV: CSV
using CairoMakie
using Measurements: measurement, value, uncertainty
using LsqFit: curve_fit, stderror
include("functions.jl")

# %%
# Load the data
cd(@__DIR__)
df = DataFrame(CSV.File("../data/01 johnson noise rt.csv"))

scale_measurements!(df) # estimates V^2 from measured values
estimate_noise_VJ2!(df)	# estimates VJ^2 by removing the amplifier noise

df.Δf = 1e3*df.Δf # convert to Hz

# fit a line to estimate k_b=1/(4TΔf) VJ²/R
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [1.0, 0.0]
T = measurement(25.0, 0.1) + 273.15
fit_params = fit_groups(df, :Δf, :R, :VJ2, mdl, p0)
fit_params.Δf = unique(df.Δf)
fit_params.kb = 1 ./ (4T*measurement.(fit_params.Δf, .04*fit_params.Δf)) .* hcat(fit_params.p...)[2, :]

# plot the noise over the resistance
f = Figure()
a = Axis(f[1, 1];
	xlabel="Resistance [Ω]", ylabel="V² Junction [V²]",
	# xscale=log10, yscale=log10
)
# marker = [MarkerElement(marker=m) for m in [:x, :o]]
colors = [PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)]
for d in eachrow(df)
	scatter!(d.R, value(d.VJ2), 
		color=d.Δf, colorrange=extrema(df.Δf)
	)
end
for f in eachrow(fit_params)
	x = range(0, 10e3, 100)
	lines!(x, mdl(x, value.(f.p)), color=f.Δf, colorrange=extrema(df.Δf))
end

axislegend(a,
	[PolyElement(color=c, colorrange=extrema(df.Δf)) for c in unique(df.Δf)],
	string.(unique(df.Δf)),
	"Δf",
	position=:lt
)
xlims!(low=0)
print(fit_params)
f
