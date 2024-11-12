using DataFrames: DataFrame, nrow, eachrow, ncol
using Statistics: mean, std
using LsqFit: curve_fit
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
f = Figure()
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
f = Figure()
a = Axis(f[1,1], xlabel="x in μm", ylabel="y in μm", aspect = DataAspect())
for p in eachrow(df)
	t = p.t .- c
	t = scale(t)
	# t = t .- mean(t, dims=1)
	# t = t .- t[1,:]'
	lines!(t, label=format(p.ot, precision=2), alpha=.75, color=p.ot, colorrange=extrema(df.ot))
end
axislegend(position=:lt, "Trap Strength")
resize_to_layout!(f)
save("../figures/01_02_1_trajectories.pdf", f)
f

# %% plot the time series
f = Figure()
a = Axis(f[1,1], 
	yscale=log10, 
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
f

# %% define msd function
i = df[2, :]
norm(i.t .- c)

# %% calculate the mean squared displacement
function msd(t)
	t = scale(t)
	t = remove_drift(t)
	m = imsd(t)
	m = mean(m, dims=2)[:,1] # mean over x and y	
	m = abs.(m)
	return m[1:1000]
end

# model
# D₀, clip
diffusion_model(x, p) = 1 ./(1 ./p[1] .* x.^-1 .+ 1 ./p[2])
function diffusion_label(p)
	if p[2] > 100
		return format(p[1]*1e3, precision=0)*" nm²/s"
	end
	return format(p[1]*1e3, precision=0)*" nm²/s, "*format(p[2], precision=2)*" nm²"
end


#%% plot the mean squared displacement
f = Figure()
a = Axis(f[1, 1], 
	xlabel="τ in s", 
	ylabel=rich("msd in μm", superscript("2")),
	xscale=log10, 
	yscale=log10,
)

# measurements
plots_s = []
plots_f = []
df.msd_inf = fill(measurement(0., 0.), size(df, 1))
df.msd_D0 = fill(measurement(0., 0.), size(df, 1))
for j in 1:size(df, 1)
	i = df[j, :]
	m = msd(i.t)
	s = scatter!(times(m), m, label=format(i.ot, precision=2), color=i.ot, colorrange=extrema(df.ot))
	mf = curve_fit(diffusion_model, times(m), m, [1e-2, 1])
	ml = lines!(
		a, times(m), diffusion_model(times(m), mf.param), 
		color=s.color, linestyle=:dash, colorrange=extrema(df.ot),
		label=diffusion_label(mf.param)
	)
	df.msd_inf[j] = measurement(mf.param[2], sqrt(mf.jacobian[2,2]))
	df.msd_D0[j] = measurement(mf.param[1], sqrt(mf.jacobian[1,1])) * 1e3
	push!(plots_s, s)
	push!(plots_f, ml)
end
axislegend(a, plots_s, [p.label for p in plots_s], "Trap Strength", position=:lt)
axislegend(a, plots_f, [p.label for p in plots_f], "Fits", position=:rb)
ylims!(a, 1e-2, 1e2)
xlims!(1, 1e2)
save("../figures/01_02_2_msd.pdf", f)
f

#%% save the data
using JLD2 
save("../data/TLM/02.jld2", "df", df)