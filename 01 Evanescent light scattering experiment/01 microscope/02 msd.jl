using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit, stderror, LsqFitResult
using LinearAlgebra: norm, normalize
using CairoMakie
import CSV, Images
using KernelDensity: kde
using Glob: glob
using MeanSquaredDisplacement: imsd
using Interpolations: interpolate, linear_interpolation
using Format: format
using Measurements: measurement

include("functions.jl")
include("../functions.jl")

# %% load the trajectories
cd(@__DIR__)
scale(t, s=0.13319672) = t .* s # px to μm
times(t) = (0:size(t,1)-1) ./ 10
df = load_trajectories()
c = mean(df.t[end], dims=1)
df

#%% plot the trajectories in different plots
f = Figure(size=fullsize)
axs = [Axis(f[i,j]) for i in [1,2] for j in [1,2]]
l = Nothing
for (a, p) in zip(axs, eachrow(df))
	t = p.t .- c
	t = scale(t)
	l = lines!(a, t, color=times(t)./60, alpha=.5)
	a.title = format(p.ot, precision=2)
end
Colorbar(f[1:2,3], l, label="Time in min")
save("../figures/01_02_11_trajectories_time.pdf", f)
f
#%% plot the trajectorie
f = Figure(size=halfsize)
a = Axis(f[1,1], xlabel="x in μm", ylabel="y in μm", aspect = DataAspect())
for p in eachrow(df)
	t = p.t .- c
	t = scale(t)
	# t = t .- mean(t, dims=1)
	# t = t .- t[1,:]'
	lines!(t, label=format(p.ot, precision=2), alpha=.75, color=p.ot, colorrange=extrema(df.ot))
end
# axislegend(position=:lt, "Trap Strength")
# Colorbar(f[1,2], limits=extrema(df.ot), label="Trap Strength")
Legend(f[1,2], a, label="Trap Strength", framevisible=false)
colgap!(f.layout, 1)
resize_to_layout!(f)
save("../figures/01_02_1_trajectories.pdf", f)
f

# %% plot the time series
f = Figure()
a = Axis(f[1,1], 
	# yscale=log10, 
	ylabel="Center distance μm", 
	xlabel="time s"
)
for p in eachrow(reverse(df))
	t = p.t .- c
	t = scale(t)
	lines!([norm(t[i,:]) for i in 1:size(t,1)], color=p.ot, colorrange=extrema(df.ot), label=format(p.ot, precision=2))
end
axislegend()
ylims!(.1, nothing)
save("../figures/01_02_2_center_distances.pdf", f)
f


# %% calculate the mean squared displacement
function msd(t)
	t = scale(t)
	t = remove_drift(t)
	m = imsd(t)
	m = sqrt.(m[:,1] .^2 .+ m[:,2] .^2)
	return m[1:100]
end

# model
# D₀, clip
diffusion_model(x, p) = 1 ./(1 ./p[1] .* (x ).^-1 .+ 1 ./p[2])
# diffusion_model(x, p) = 1 ./((p[1].*x).^(-p[3]) .+ p[2].^(-p[3])).^(1/p[3])
# diffusion_model(x, p) = p[1] .* x .+ p[2]

# plot the mean squared displacement
f = Figure(size=halfsize)
a = Axis(f[1, 1], 
	xlabel="τ in s", 
	ylabel=rich("MSD(τ) in μm", superscript("2")),
	xscale=log10, 
	yscale=log10,
)

# measurements
msd_fit = []
for j in 1:size(df, 1)
	i = df[j, :]
	m = msd(i.t)
	s = scatter!(times(m), m, label=format(i.ot, precision=2), color=i.ot, colorrange=extrema(df.ot), markersize=5, alpha=.5)
	d0 = mean(diff(m[1:5]) ./ diff(times(m)[1:5]))
	mf = curve_fit(diffusion_model, times(m), m, [d0, 1., 1], lower=[0., 0., 0.])
	push!(msd_fit, mf)
	ml = lines!(
		a, times(m), clamp.(diffusion_model(times(m), mf.param), 0, Inf), 
		color=s.color, colorrange=extrema(df.ot),
	)
	println(2*kT./last(m) .* 1e12 .* 1e9)
end
Legend(f[1,2], a, framevisible=false)
colgap!(f.layout, 0)
xlims!(a, 1e-1, 1e1)
ylims!(a, .05, 1e1)
save("../figures/01_02_2_msd.pdf", f)

function fit_param(m, i)
	try
		return measurement(m.param[i], stderror(m)[i])
	catch
		return measurement(m.param[i], NaN)
	end
end

df.msd_fit = msd_fit
df.msd_inf = [fit_param(m, 2) for m in msd_fit]
df.msd_D0 = [fit_param(m, 1) for m in msd_fit]

T = measurement(22.6, 0.1) + 273.15
kB = 1.38064852e-23
kT = kB * T
df[:, :k_msd] = 2*kT./df.msd_inf .* 1e12 .* 1e9

f
#%%
using JLD2 
save("../data/TLM/02.jld2", "df", df)