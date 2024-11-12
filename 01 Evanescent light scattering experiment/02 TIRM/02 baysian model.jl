# find the distribution of z from the intensity data
# using following assumptions
# I(z) = I0 * exp(- z / β)
# Δz = N(0, σ^2)

using Turing

# %% define a model
@model function model(I)
	# Priors for the parameters
	I0 ~ Truncated(Normal(0, 10), 0, Inf)
	beta ~ LogNormal(0, 1)
	σdz ~ truncated(Normal(0, 1), 0, Inf)

	dZ = Normal(0, σdz)

	

	for i in eachindex(I)
		I[i] 
	end

	# apply a random walk model
	# z = Vector{Real}(undef, length(I))
	# z[1] = 0
	# for i in 2:length(I)
	# 	z[i] ~ z[i-1] + Normal(0, σdz)
	# end

	# # apply the model
	# for i in eachindex(I)
	# 	I[i] ~ Normal(I0 * exp(-z[i] / beta), .1)
	# end
end


# test the model
dz = rand(Normal(5, .5), 100)
z = cumsum(dz)
I = 10 .* exp.(-z ./ 1.5)

chain = sample(model(I), NUTS(), 100)

# %%
chain