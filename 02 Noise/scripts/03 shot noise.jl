#=
Script to analyse the amplified Shot noise to estimate e.
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
# df = DataFrame(CSV.File("../data/02 johnson noise temperature.csv"))

scale_measurements!(df) # estimates V^2 from measured values
estimate_noise_VJ2!(df)	# estimates VJ^2 by removing the amplifier noise

# fit lines to estimate T0
groups = groupby(df, :Δf)
mdl(x, p) = p[1] .+ p[2] .* x
p0 = [1.0, 0.0]
fit_params = DataFrame(Dict(
	:p => fill(measurement.(p0), length(groups)),
	:Δf => [g.Δf[1] for g in groups],
))
for (g, f) in zip(groups, eachrow(fit_params))
	p = curve_fit(mdl, value.(g.T), value.(g.VJ2), p0)
	f.p = measurement.(p.param, stderror(p))
end
fit_params[:, :T0] = [-p[1]/p[2] for p in fit_params.p]


# plot the VJ2 over T fir different R and Δf
f = Figure()
a = Axis(f[1, 1];
	xlabel="Temperature in K", 
	ylabel="V² Junction",
)
for g in groupby(df, :Δf)
		scatter!(
			value.(g.T), value.(g.VJ2),
			color=g.Δf[1], colorrange=extrema(df.Δf),
			label="Δf=$(g.Δf[1]) Hz"
		)
end
for f in eachrow(fit_params)
	x = range(0, 300, length=100)
	lines!(x, mdl(x, value.(f.p)), 
		color=f.Δf, colorrange=extrema(df.Δf)
	)
end
Legend(f[1,2], a, framevisible=false)
xlims!(low=0)
ylims!(low=0)
f