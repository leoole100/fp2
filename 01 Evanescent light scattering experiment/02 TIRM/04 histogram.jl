using CairoMakie
using StatsBase: moment, std, mean
using KernelDensity: kde

# load and filter data
include("functions.jl")

# %% look at same z
f = Figure()
a = Axis(f[1,1],
xlabel="dz", ylabel="pdf"
)
for i in eachrow(df)
	z = z_estimate(i.I)
	dz = diff(z)
	k = kde(dz, boundary=(-.5,.5),bandwidth=.025)
	lines!(k.x, k.density,
	color=i.ot, colorrange=extrema(df.ot)
	)
end
f

# %%
diff_z(z, d) = z .- circshift(z, d)

# the main assumption is that D(z0) = 1/2 dσ²(z,τ)/dτ
# ∫dτ da(τ)/dτ is computed as mean(diff(a))
function dσdt(z; z_bins=10, τ=1:10)
	z0 = range(minimum(z), maximum(z), length=z_bins+1)
	dzdt = zeros(length(z0)-1)
	dz = [diff_z(z, i) for i in τ]
	for i in 1:length(z0)-1
		z_mask = z0[i] .< z .< z0[i+1]
		if sum(z_mask) > 3
			std_dz = [std(dz_i[z_mask]) for dz_i in dz]
			dzdt[i] = mean(std_dz)
		else
			dzdt[i] = NaN
		end
	end
	return z0[1:end-1], dzdt
end
function D(dσdt)
	z, σ = dσdt
	D = σ.^2 ./ (2*Δt)
	return z, D
end

f = Figure()
a = Axis(f[1,1], 
	xlabel="z in μm", ylabel="D(z) in μm²/s",
	yscale=log10,
)
for i in eachrow(df)
	# z = z_estimate(i.I, I0=2.5, β=0.5)
	z = z_estimate(i.I, I0=2, β=0.5)
	z, dz = D(dσdt(z))
	lines!(a, z, dz.^2, 
		color=i.ot, colorrange=extrema(df.ot),
		label=i.p
	)
end
lines!(
	.01:.01:1, 
	z -> Dz_theoretical(z),
	color=:black,
	label="Theoretical"
)
# Legend(f[1,2], a, label="Trap Strength", framevisible=false)
# axislegend(position=:lt)
colgap!(f.layout, 1)
f