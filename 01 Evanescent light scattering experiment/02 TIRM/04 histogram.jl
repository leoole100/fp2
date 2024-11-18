using CairoMakie
using StatsBase: moment, std, mean, var
import StatsBase
using KernelDensity: kde
using Optim
using Format: format

# load and filter data
include("functions.jl")

diff_z(z, d) = circshift(z, -d) .- z

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

function D(dσdt; mask=true)
	z, v, c = dσdt
	D = v ./ (2*Δt)
	if mask
		mask = z .< 0.4 .&& isnan.(D) .== false
		return z[mask], D[mask], c[mask]
	end
	return z, D, c
end


# %%
df[:, :fit] .= fill([0.0, 0.0], size(df, 1))

# add beta from df.l column
l_beta = Dict(5.0=>.454, 10.0=>.322, 15.0=>.265, 20.0=>.230)
df[:, :β] = zeros(size(df, 1))
for i in 1:size(df, 1)
	df[i, :β] = l_beta[df[i, :l]+5.0]
	df[i, :fit] = [1.0, df[i, :β]]
end
df[4, :fit] = [1.3, .2]
df[5, :fit] = [.9, .27]


df[:, :Brenner] = fill(false, size(df, 1))
df[4, :Brenner] = true
df[5, :Brenner] = true

f = Figure(size=fullsize)
a = Axis(f[1:2,1], 
	xlabel="z in μm", ylabel="D(z) in μm²/s",
	# yscale=log10,
	# xscale=log10,
)
lines!(
	0:.01:.5, Dz_theoretical,
	color=:black,
	label="Model",
)
for (i,r) in enumerate(eachrow(df[[4,5, 1,2,3, 6],:]))
# for i in eachrow(df)
	# z = z_estimate(i.I, I0=p[1], β=p[2])
	z = z_estimate(r.I, I0=r.fit[1], β=r.fit[2])
	z, d, c = D(dvdt(z, cutoff=1/200))
	if r.fit[2] == r.β
		β_string = "β = "*string(r.β)
	else
		β_string = "β = "*string(r.β)*" ("*format(r.fit[2])*")"
	end
	scatterlines!(z, d,	
		color=r.ot, colorrange=extrema(df.ot),
		label="I₀ = "*string(r.fit[1])*", "*β_string,
		# markersize=7 .*c./maximum(c) .+ 3,
		# alpha=c./maximum(c),
		marker=[:circle, :rect, :diamond][mod(i, 3)+1],
	)
end
Legend(f[1,2], a, framevisible=false)
Colorbar(f[2,2], limits=extrema(df.ot), label="Trap Strength", vertical=false, flipaxis=false)
resize_to_layout!(f)
save("../figures/02_04_01_diffusion.pdf", f)
f


# %% look at dz distribution
f = Figure(size=fullsize)
aI = Axis(f[1,1], xlabel="I", ylabel="pdf")
adI = Axis(f[1,2], xlabel="dI/dt")
az = Axis(f[2,1], xlabel="z", ylabel="pdf")
adz = Axis(f[2,2], xlabel="dz/dt")

k =nothing
for (i, r) in enumerate(eachrow(df[[5,6], :]))
	l ="β = "*format(r.β, precision=3)*" μm"
	l *= "\nI₀ = "*format(r.fit[1], precision=1)
	if r.Brenner 
		l *= " (fit)"
	end

	k = dist(r.I, cutoff=0.05)
	lines!(aI, k.x, k.y, label=l)
	k=kde(diff(r.I), boundary=(-.05, .05))
	lines!(adI, k.x, k.density)
	
	z = z_estimate(r.I, I0=r.fit[1], β=r.fit[2])
	k = dist(z, cutoff=0.05)
	lines!(az, k.x, k.y)

	k=kde(diff(z), boundary=(-.03, .03), bandwidth=.001)
	lines!(adz, k.x, k.density)
end
for a in [aI, adI, az, adz]
	ylims!(a, low=0)
end
Legend(f[:,3], aI, framevisible=false)
resize_to_layout!(f)
save("../figures/02_04_02_hist.pdf", f)
f