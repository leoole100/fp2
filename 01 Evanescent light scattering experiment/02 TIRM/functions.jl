include("../functions.jl")

using JLD2
using DataFrames

# Constants
# from previous analysis
# D0=0.122
# R=1.86

# from matlab script
R=2.15
D0=0.11195

Δt=1e-3

# Define functions
# remove_offset(I) = I .- minimum(vcat(df.I...))
remove_offset(I) = abs.(I)
# calculate z
z_estimate(I; β=1, I0=1) = -log.(remove_offset(I)./I0) .* β
Dz_theoretical(z; D0=D0, R=R) = D0 ./ (R./z + 0.2 .* log.(R./z) + 0.9712)

# Load data
cd(@__DIR__)
df = load("../data/TIRM/data.jld2")["df"]
df = filter(i -> i.ot<1.1 || i.ot==1.5, df)
df