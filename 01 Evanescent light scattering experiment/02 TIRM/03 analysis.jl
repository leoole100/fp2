using CairoMakie
using JLD2
using StatsBase: moment, std
using KernelDensity: kde
using Optim: optimize, GradientDescent
using DataFrames

# load data
cd(@__DIR__)
df = load("../data/TIRM/data.jld2")["df"]
df

D0=0.122
R=1.86
β=0.454

# %% Define functions
remove_offset(I) = I .- minimum(vcat(df.I...))
# calculate z
z_estimate(I; β=1, I0=maximum(I)) = log.(I0 ./ remove_offset(I)) ./ β
dz_estimate(I; β=1, I0=maximum(I)) = diff(z_estimate(I, β, I0))
D_estimate(std_z; Δt=1e-3) = std_z.^2 ./ (2*Δt)  # from script
function Dz_estimate(I; I0=maximum(I), β=1)
	z = z_estimate(I, I0=I0, β=β)
	dz = diff(z)
	z = z[2:end]
	
	zdz = DataFrame(z=z, dz=dz)
	groups =  groupby(zdz, :z)
	groups = filter(x -> nrow(x) > 100, groups)
	std_dz = combine(groups, :dz => std => :std_dz)
	filter!(i->i.z<3, std_dz)
	return  std_dz.z, D_estimate(std_dz.std_dz)
end

Dz_theoretical(z; D0=D0, R=R) = D0 ./ (R./z + 0.2 .* log.(R./z) + 0.9712)

# %%
f = Figure()
a = Axis(f[1,1], 
	xlabel="z/β", ylabel="D(z) in β^2/s",
	# xscale=log10,
	# yscale=log10
)
Dz = Nothing
for i in eachrow(df)
	Dz = Dz_estimate(i.I, I0=1)
	scatter!(reverse(Dz)...,
		color=i.ot, colorrange=extrema(df.ot),
		label=i.p,
		alpha=.5
	)
end
lines!(
	0:.1:10,
	z -> Dz_theoretical(z)*20,
	label="model",
	color=:black
)

Legend(f[1,2], a, framevisible=false)
colgap!(f.layout, 0)
f

