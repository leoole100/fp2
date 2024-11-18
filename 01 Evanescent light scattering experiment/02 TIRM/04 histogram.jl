using CairoMakie
using StatsBase: moment, std, mean, var
using KernelDensity: kde

# load and filter data
include("functions.jl")

diff_z(z, d) = circshift(z, -d) .- z
mean_slope(y) = mean(diff(y))

# similar to the provided matlab code
function dvdt(z; z_bins=100, t_bins=5, cutoff=1/100)
	z0 = range(extrema(z)..., length=z_bins+1)
	z_step = diff(z0)[1]
	dzdt = zeros(length(z0))
	count = zeros(length(z0))
	dz = [diff_z(z, i) for i in 1:t_bins]
	for i in 1:length(z0)
		z_mask = z0[i]-z_step/2 .< z .< z0[i]+z_step/2
		count[i] = sum(z_mask)
		if sum(z_mask) > length(z).*cutoff
			std_dz = [var(dz_i[z_mask]) for dz_i in dz]
			dzdt[i] = mean_slope(std_dz)
		else
			dzdt[i] = NaN
		end
	end
	return z0, dzdt, count
end
function D(dσdt)
	z, v, c = dσdt
	D = v ./ (2*Δt)
	mask = z .< .4
	return z[mask], D[mask], c[mask]
end

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

p = (I0=1.2, β=.156)
# p = (I0=.85, β=.016) #fits but only on linear part
# p = (I0=2.5, β=0.3) # expected

f = Figure()
a = Axis(f[1,1], 
	xlabel="z in μm", ylabel="D(z) in μm²/s",
	# yscale=log10,
	# xscale=log10,
	title="I₀="*string(p.I0)*", β="*string(p.β)*" μm"
)
lines!(
	0:.01:.7, Dz_theoretical,
	color=:black,
	label="Theoretical",
)
for i in eachrow(df[4:4, :])
# for i in eachrow(df)
	z = z_estimate(i.I, I0=p.I0, β=p.β)
	z, d, c = D(dvdt(z, cutoff=1/50))
	scatter!(z, d, 
		color=i.ot, colorrange=extrema(df.ot),
		label=i.p,
		markersize=10 .*c./maximum(c).+2,
		alpha=.7
	)
end
# Legend(f[1,2], a, label="Trap Strength", framevisible=false)
axislegend(position=:lt)
colgap!(f.layout, 1)
f