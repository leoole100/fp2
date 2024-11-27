# evaluate the precision of the multimeter

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
df = DataFrame(CSV.File("../data/07 unceritenty multimeter.csv"))
df.V = df.Voltage
select!(df, :V)

# %%
f = Figure()
a = Axis(f[1, 1]; xlabel="Voltage in V", ylabel="Frequency in Hz")
density!(df.V)
ylims!(low=0)
f

# %%
meter_std = std(df.V)
save("../data/gen/04 meter uncertainty.jld2", Dict("std"=>meter_std))
meter_std