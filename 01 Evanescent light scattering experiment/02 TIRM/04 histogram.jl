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
# similar to the provided matlab code
function dσdt(z; z_bins=100, τ=1:10, cutoff=1/100)
	z0 = range(minimum(z), maximum(z), length=z_bins+1)
	dzdt = zeros(length(z0)-1)
	count = zeros(length(z0)-1)
	dz = [diff_z(z, i) for i in τ]
	for i in 1:length(z0)-1
		z_mask = z0[i] .< z .< z0[i+1]
		count[i] = sum(z_mask)
		if sum(z_mask) > length(z).*cutoff
			std_dz = [std(dz_i[z_mask]) for dz_i in dz]
			dzdt[i] = mean(std_dz)
		else
			dzdt[i] = NaN
		end
	end
	return z0[1:end-1], dzdt, count
end
function D(dσdt)
	z, σ, c = dσdt
	D = σ.^2 ./ (2*Δt)
	return z, D, c
end
D(dσdt(z_estimate(df[1,:].I)))

#%%
p = (I0=.85, β=.016)

f = Figure()
a = Axis(f[1,1], 
	xlabel="z in μm", ylabel="D(z) in μm²/s",
	# yscale=log10,
	# xscale=log10,
	title="I₀="*string(p.I0)*", β="*string(p.β)*" μm"
)
lines!(
	0:.01:.1, 
	z -> Dz_theoretical(z),
	color=:black,
	label="Theoretical",
)
# for i in eachrow(df[end:end, :])
for i in eachrow(df)
	# z = z_estimate(i.I, I0=2.5, β=0.3) # expected
	z = z_estimate(i.I, I0=p.I0, β=p.β)
	z, d, c = D(dσdt(z, cutoff=1/40))
	scatter!(z, d, 
		color=i.ot, colorrange=extrema(df.ot),
		label=i.p,
		markersize=c./maximum(c).*10 .+ 5,
		alpha=.7
	)
	# fit a line
	j = argmax(c)
	m = diff(d)[j] / diff(z)[j]
	b = d[j] - m*z[j]
	lines!(0:.01:z[j+2], z->m.*z.+b, color=i.ot, colorrange=extrema(df.ot))
end
# Legend(f[1,2], a, label="Trap Strength", framevisible=false)
# axislegend(position=:lt)
colgap!(f.layout, 1)
f