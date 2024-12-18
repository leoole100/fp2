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
# fit a plane for V over R and Δf
plane(x, p) = p[1].*x[:, 2] .+ p[2] .* x[:, 1] .* x[:, 2]
p = curve_fit(
	plane, [df.R df.Δf], value.(df.V), [0.0, 1.0]
)
p = measurement.(p.param, stderror(p))



# %%
# Plot in 3d
f = Figure(size=fullsize)
a = Axis3(f[1, 1]; 
	xlabel="R in kΩ",
	ylabel="Δf in kHz",
	zlabel="(μV)²",
)

x = range(0, stop=maximum(df.R), length=20)
y = range(0, stop=maximum(df.Δf), length=20)
z = [p[1] .* y .+ p[2] .* x .* y for x in x, y in y]
contour3d!(x./1e3, y./1e3, value.(z).*1e12, levels=20, linewidth=2)
# wireframe!(x, y, value.(z))
surface!(x./1e3, y./1e3, value.(z).*1e12, alpha=.5)

# plot the data
residuals = value.(df.V) .- plane([df.R df.Δf], p)
s = stem!(
	df.R./1e3, df.Δf./1e3, value.(df.V).*1e12,
	strokecolor=:black, strokewidth=1.5,
	markersize=12,
	color=value.(residuals).*1e12, 
	colorrange=value.((-maximum(abs.(residuals)), maximum(abs.(residuals)))).*1e12,
	colormap=:coolwarm
	# color= value.(df.V), colorrange=colorrange
)
Colorbar(f[1, 2], s, label="Residuals")

xlims!(low=0); ylims!(low=0); zlims!(low=0)

save("../figures/05 johnson noise rt plane.pdf", f)
f

# %%
# calculate kb T from the slope
T = measurement(22.0, 3) + 273.15
kb = p[2] / 4 / T

save("../data/gen/05 Johnson noise RT.jld2", Dict(
	"S0"=>p[1],
))