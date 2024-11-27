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
df = DataFrame(CSV.File("../data/04 photocurrent.csv"))

# df.V = measurement.(df.Vsq, df.VsqU) 
df.V = measurement.(df.Vsq, load("../data/gen/04 meter uncertainty.jld2")["std"]) 
df.V = df.V .* 10 ./ (100 .* df.G2 .*1000).^2 # scale measurements to volts²
df.V *= 1/(10000)^2 	# convert to current, with a 10kΩ resistor
df.Δf = 1e3*df.Δf # convert to Hz
df.I *= 1e-6 # to A
df.S = df.V ./ df.Δf
sort!(df, [:I, :Δf])

select(df, [:I, :Δf, :S])


#%%

e = vcat([diff(d.S) ./ diff(d.I) / 2 for d in groupby(df, :Δf)] ...)

f = Figure()
a = Axis(f[1,1],
	xlabel="measured e in 10^-19 C"
)
density!(value.(e).*1e19)
vlines!([1.602e-19].*1e19, color=:black)
ylims!(low=0)
hideydecorations!(a)
f