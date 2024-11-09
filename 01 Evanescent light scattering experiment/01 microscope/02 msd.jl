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

# %% load the trajectories
path = "../data/TLM/"
cd(@__DIR__)

function load_trajectory!(df, path)
	path = join(split(path, ".")[1:end-1], ".")
	f = DataFrame(CSV.File(path * ".csv"))
	for i in 1:ncol(f)÷2
		push!(df[!, :t], hcat(f[:,2*i-1], f[:,2*i]))
		n = split(path, "/")[end]
		n = join(split(n, ".")[1:end-1], ".")
		n = split(n, " ")[3]
		ot = parse(Float64, n[3:end])
		push!(df[!, :n], n)
		push!(df[!, :ot], ot)
	end
end

df = DataFrame(t=[], n=[], ot=[])
paths = glob(path * "*.trajectory.csv")
for p in paths
	load_trajectory!(df, p)
end

scale(t, s=0.13319672) = t .* s # px to μm
times(t) = 0:length(t)-1 ./ 10

sort!(df, :ot)
df

#%% plot the trajectories
f = Figure()
a = Axis(f[1,1], xlabel="x in μm", ylabel="y in μm")
for p in eachrow(df)
	t = scale(p.t)
	t = t .- mean(t, dims=1)
	lines!(t, label=p.n, alpha=0.5)
	# points = [Point2f(p.trajectory[j,:]...) for j in p.frames]
	# datashader!(points)
end
axislegend(position=:lt)
save("../figures/01_trajectories.pdf", f)
f

# %% calculate the mean squared displacement
function remove_drift(t)
	t = copy(t)
	t = t .- t[1,:]'
	s = t[end,:] ./ (size(t,1)-1)
	for i in 0:(size(t,1)-1)
		t[i+1,:] = t[i+1,:] .- i .* s
	end
	return t
end

function msd(t)
	t = scale(t)
	t = remove_drift(t)
	m = imsd(t)
	m = mean(m, dims=2)[:,1] # mean over x and y	
	m = abs.(m)
	return m[1:1000]
end

interp_msd = linear_interpolation(
	df.ot,
	msd.(df.t)
)

# model
# nΓ, clip
diffusion_model(x, p) = clamp.(p[1] * x, 0, p[2])
diffusion_label(p)= "min("*format(mf.param[1]*1e3, precision=0)*" nm²/s ⋅ τ, "*format(mf.param[2], precision=2)*")"

#%% plot the mean squared displacement
f = Figure()
a = Axis(f[1, 1], 
	xlabel="τ in s", 
	ylabel=rich("msd in μm", superscript("2")),
	xscale=log10, 
	yscale=log10,
)

# interpolated values
# for ot in range(0.7, 1.0, 10)
# 	m = interp_msd(ot)
# 	s = lines!(times(m), m, color=ot, colorrange=extrema(df.ot))
# end

# measurements
plots_s = []
plots_f = []
for i in eachrow(filter(d -> d.ot<2, df))
	m = msd(i.t)
	s = scatter!(times(m), m, label=format(i.ot, precision=2), color=i.ot, colorrange=extrema(df.ot))
	mf = curve_fit(diffusion_model, times(m), m, [1e-2, 1])
	ml = lines!(
		a, times(m), diffusion_model(times(m), mf.param), 
		color=s.color, linestyle=:dash, colorrange=extrema(df.ot),
		label=diffusion_label(mf.param)
	)
	push!(plots_s, s)
	push!(plots_f, ml)
end
axislegend(a, plots_s, [p.label for p in plots_s], "Trap Strength", position=:lt)
axislegend(a, plots_f, [p.label for p in plots_f], "Fits", position=:rb)

# # free diffusion model
# m = msd(df.t[1])
# t = times(m)
# mf = curve_fit(diffusion_model, t, m, [1e-2, 1, 100])
# ml = lines!(a, t, diffusion_model(t, mf.param), 
# 	color=:black, linestyle=:dash,
# )
# axislegend(a, [ml], 
# 	[diffusion_label_free(mf.param)],
# 	position=:rb,
# )

# Colorbar(f[1, 2], limits=extrema(df.ot), label="Optical trap strength")
ylims!(a, 1e-2, 1e2)
save("../figures/01_msd.pdf", f)
f
