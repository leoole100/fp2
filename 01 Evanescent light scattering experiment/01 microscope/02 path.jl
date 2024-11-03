# find particles in images and track them

using DataFrames, Statistics, LsqFit, LinearAlgebra
using WGLMakie
import CSV, Images

path = "../data/test/particles"

# %%
# open the trajectory
cd(@__DIR__)
img = Images.load(path * ".mean.png");
f = DataFrame(CSV.File(path * ".trajectory.csv"))

df = DataFrame(trajectory=[])

trajectory = []
for i in 1:ncol(f)÷2
	push!(df[!, :trajectory], hcat(f[:,2*i-1], f[:,2*i]))
end
df[!, :frames] = [1:size(p.trajectory)[1] for p in eachrow(df)]

df

# %%
# plot the trajectories on the image
f = Figure()
image(f[1,1], img', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)
for p in eachrow(df)
	scatter!( p.trajectory)
	# points = [Point2f(p.trajectory[j,:]...) for j in p.frames]
	# datashader!(points)
end
f

# %%
# calculate the mean squared displacement

# df[!, :reference] = [mean(p.trajectory, dims=1) for p in eachrow(df)]
df[!, :reference] = [p.trajectory[1, :] for p in eachrow(df)];
df[!, :msd] = [[sum((p.trajectory[j,:] .- p.reference).^2) for j in p.frames] for p in eachrow(df)];

linear_model(x, p) = clamp.(p[2] .* x, 0, p[3]);
df[!, :linear_model] = [curve_fit(linear_model, p.frames, p.msd, [0.0, 0.0, 1.0]) for p in eachrow(df)];

exp_model(x, p) = p[1] .- p[1] .* exp.(-p[2] .* x ./ p[1]);
df[!, :exp_model] = [curve_fit(exp_model, p.frames, p.msd, [50.,1.]) for p in eachrow(df)];


f = Figure()
a = Axis(f[1, 1], xlabel="time", ylabel="mean squared displacement")
for i in eachrow(df)
	plot!(a, 1:size(i.trajectory,1), i.msd)
	if i.linear_model.converged
		# lines!(a, 1:size(i.trajectory,1), linear_model(1:size(i.trajectory,1), i.linear_model.param))
	end
	if i.exp_model.converged
		lines!(a, 1:size(i.trajectory,1), exp_model(1:size(i.trajectory,1), i.exp_model.param))
	end
end
# try models
lines!(0:100, x->linear_model(x, [0, 100., 1000.]), label="linear")
lines!(0:100, x->exp_model(x, [1000., 100.]), label="exp")
axislegend(a)
f

# %%
# calculate the probability distribution of the displacements
df[!, :center] = [mean(p.trajectory, dims=1) for p in eachrow(df)]
df[!, :distance] = [[norm(p.trajectory[j,:]' .- p.center) for j in p.frames] for p in eachrow(df)]
df[!, :σ] = [std(p.distance) for p in eachrow(df)]

f = Figure()
a = Axis(f[1, 1], xlabel="distance from center", ylabel="probability")
for i in 1:nrow(df)
	p = df[i, :]
	offset = i
	hist!(p.distance, normalization=:probability, offset=offset, scale_to=0.6, )
	lines!(a, 0:0.1:maximum(p.distance), x->exp(-x^2/(2*p.σ^2))*.6 +offset,)
end
f

# %%
# calculate the probability distribution for a given potential
