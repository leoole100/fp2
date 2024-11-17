# find the distribution of z from the intensity data
# using following assumptions
# I(z) = I0 * exp(- z / β)
# Δz = N(0, σ^2)

using Turing
using CairoMakie
using Distributions
using KernelDensity
using StatsBase

# %% generate test data
I0 = 1
β = 100
sigma = .1
z0 = 0.1


dz = rand(Normal(z0, sigma), 1000)
z = cumsum(dz)
I = I0 .* exp.(-z ./ β)
round_to_number(x, n) = round.(Int, x ./ n) .* n
# I = round_to_number(I, .5)

f = Figure()
a = Axis(f[1, 1])
density!(I)
ylims!(low=0)
f

# %% define a model
@model function model(I_dist)
	σ ~ truncated(Normal(0, 1), 0, Inf)
	β ~ truncated(Normal(1, 1), 0, Inf)
	I0 ~ truncated(Normal(10, 1), 0, Inf)
	dZ ~ Normal(0, σ)

	Z ~ truncated(Normal(1, 1),0, Inf)
	I_mdl = I0 .* exp.(-Z ./ β)
	I_mdl ~ I_dist
end

I_dist = fit(Normal, I)
k = kde(I)
I_dist = pdf(k)
kde_model

chain = sample(model(I_dist), NUTS(), 1000)


# %% plot the sampled density

