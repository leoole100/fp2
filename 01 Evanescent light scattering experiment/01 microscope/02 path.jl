# find particles in images and track them

using DataFrames, Statistics, LsqFit
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

df


f = Figure()
image(f[1,1], img', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)
for p in eachrow(df)
	scatter!(f[1,1], p.trajectory)
end
f

# %%
# calculate the mean squared displacement
model(x, p) = clamp.(p[2] .* x, 0, p[3])

# df[!, :reference] = [mean(p.trajectory, dims=1) for p in eachrow(df)]
df[!, :reference] = [p.trajectory[1, :] for p in eachrow(df)]
df[!, :msd] = [[sum((p.trajectory[j,:] .- p.reference).^2) for j in 1:size(p.trajectory,1)] for p in eachrow(df)]
df[!, :fit] = [LsqFit.curve_fit(model, 1:size(p.trajectory,1), p.msd, [0.0, 0.0, 1.0]) for p in eachrow(df)]

f = Figure()
a =Axis(f[1, 1])
for i in eachrow(df)
	plot!(a, 1:size(i.trajectory,1), i.msd)
	if i.fit.converged
		lines!(a, 1:size(i.trajectory,1), model(1:size(i.trajectory,1), i.fit.param))
	end
end

f
