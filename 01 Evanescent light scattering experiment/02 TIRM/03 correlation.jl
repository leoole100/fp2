using CairoMakie
using StatsBase: moment, std, autocor
using KernelDensity: kde
using Optim: optimize, GradientDescent
using DSP: conv

# load and filter data
include("functions.jl")

D_estimate(std_z; Δt=Δt) = std_z.^2 ./ (2*Δt)  # from script

#%% calculate using derived correlation theorem
weights(z, z0; bandwidth=.1) = exp.(-(z .- z0).^2 ./ 2*bandwidth^2) ./ sqrt(2*π*bandwidth^2)
function Dz_measure(z, z0; bandwidth=.1)
	w = weights(z, z0, bandwidth=bandwidth)
	c = conv(w.*z, reverse(z))
	return sum(c)./ length(c)./Δt
end

Dz_list(z, z0; bandwidth=.1) = [Dz_measure(z, i, bandwidth=bandwidth) for i in z0]

i = df[1,:]
z = z_estimate(i.I, I0=2)
z0 = 0:.5:3
f = Figure()
a = Axis(f[1,1], xlabel="z", ylabel="D(z)")
plot!(z0, Dz_list(z, z0, bandwidth=.5))
# plot!(Dz_estimate(i.I)..., color=:red)
f

# %% 
# look at the auto correlation of z
f = Figure()
a = Axis(f[1,1], xlabel="τ in s", ylabel="C(τ)")
for i in eachrow(df)
	z = z_estimate(i.I)
	cor = autocor(z)
	lines!(
		range(1, length(cor)) * Δt,
		cor,
		color=i.ot, colorrange=extrema(df.ot)
	)
end
f