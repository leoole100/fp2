using CairoMakie
using StatsBase: moment, std, mean, var
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

# %% look at dz distribution
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

# %% implement optimization procedure
function error(p; I)
	z = z_estimate(I, I0=p[1], β=p[2])
	z, d, c = D(dvdt(z, cutoff=1/100), mask=false)
	theo = Dz_theoretical.(clamp.(z, 0, Inf))
	d = d .* c
	d = filter(!isnan, d)
	if length(d) == 0
		return Inf
	end
	return sum(d.^2)
end

function fit(I)
	result = optimize(p->error(p, I=I), [1.2, .2])
	return result
end

result = fit(df.I[4])
result.minimizer


# %%
df[:, :fit] .= fill([1, .2], size(df, 1))
l_beta = Dict(5.0=>.454, 10.0=>.322, 15.0=>.265, 20.0=>.230)

# add beta from df.l column
df[:, :β] = zeros(size(df, 1))
for i in 1:size(df, 1)
	df[i, :β] = l_beta[df[i, :l]+5.0]
	df[i, :fit] = [1.2, df[i, :β]]
end

df[4, :fit] = [1.3, .2]
df[5, :fit] = [.9, .27]

f = Figure(size=fullsize)
a = Axis(f[1,1], 
	xlabel="z in μm", ylabel="D(z) in μm²/s",
	# yscale=log10,
	# xscale=log10,
)
lines!(
	0:.01:.5, Dz_theoretical,
	color=:black,
	label="Theoretical",
)
for (i,r) in enumerate(eachrow(df))
# for i in eachrow(df)
	# z = z_estimate(i.I, I0=p[1], β=p[2])
	z = z_estimate(r.I, I0=r.fit[1], β=r.fit[2])
	z, d, c = D(dvdt(z, cutoff=0))
	scatter!(z, d,
		color=r.ot, colorrange=extrema(df.ot),
		label="Trap Strength: "*format(r.ot, precision=2)*" β_set="*format(r.β)*"\n"*"I₀="*string(r.fit[1])*", β_exp="*string(r.fit[2])*" μm",
		markersize=7 .*c./maximum(c) .+ 3,
		alpha=.7,
		marker=[:circle, :rect, :diamond][mod(i, 3)+1]
	)
end
Legend(f[1,2], a, framevisible=false)
colgap!(f.layout, 1)
resize_to_layout!(f)
save("../figures/02_04_01_diffusion.pdf", f)
f