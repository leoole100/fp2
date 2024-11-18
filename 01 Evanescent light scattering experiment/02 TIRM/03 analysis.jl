using CairoMakie
using JLD2
using StatsBase: moment, std, autocor
using KernelDensity: kde
using Optim: optimize, GradientDescent
using DataFrames
using DSP: conv

# load data
cd(@__DIR__)
df = load("../data/TIRM/data.jld2")["df"]
df = filter(i -> i.ot<1.1 || i.ot==1.5, df)

D0=0.122
# R=1.86
R=2.15
β=0.454
Δt=1e-3

# Define functions
# remove_offset(I) = I .- minimum(vcat(df.I...))
remove_offset(I) = abs.(I)
# calculate z
z_estimate(I; β=β, I0=1) = log.(I0 ./ remove_offset(I)) .* β
D_estimate(std_z; Δt=Δt) = std_z.^2 ./ (2*Δt)  # from script
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

Dz_theoretical(z; D0=D0/β^2, R=R/β) = D0 ./ (R./z + 0.2 .* log.(R./z) + 0.9712)


#%%
weights(z, z0; bandwidth=.1) = exp.(-(z .- z0).^2 ./ 2*bandwidth^2) ./ sqrt(2*π*bandwidth^2)
function Dz_measure(z, z0; bandwidth=.1)
	w = weights(z, z0, bandwidth=bandwidth)
	c = conv(w.*z, reverse(z))
	return sum(c)./ length(c)./Δt
end

Dz_list(z, z0; bandwidth=.1) = [Dz_measure(z, i, bandwidth=bandwidth) for i in z0]

i = df[1,:]
z = z_estimate(i.I, I0=2)
z0 = 0:.5:3
f = Figure()Dz_list(z, z0, bandwidth=.5)
a = Axis(f[1,1], xlabel="z", ylabel="D(z)")
plot!(z0, Dz_list(z, z0, bandwidth=.5))
# plot!(Dz_estimate(i.I)..., color=:red)
f

#%%
I = df.I[end]
z = z_estimate(I, I0=1, β=1)
dz = diff(z)
z = z[2:end]

zdz = DataFrame(z=z, dz=dz, t=1:length(dz))
# groups =  groupby(zdz, :z)
# groups = filter(x -> nrow(x) > 100, groups)
# std_dz = combine(groups, :dz => std => :std_dz)

# %% 
# look at the auto correlation of z
f = Figure()
a = Axis(f[1,1], xlabel="τ in s", ylabel="C(τ)")
for i in eachrow(df)
	z = z_estimate(i.I)
	cor = autocor(z)
	lines!(
		range(1, length(cor)) * Δt,
		cor,
		color=i.ot, colorrange=extrema(df.ot)
	)
end
f
# %%
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
f = Figure()
a = Axis(f[1,1], 
	xlabel="z/β", ylabel="D(z) in β^2/s",
	# xscale=log10,
	# yscale=log10
)
Dz = Nothing
for i in eachrow(df)
	Dz = Dz_estimate(i.I, I0=2.4)
	scatter!(Dz...,
		color=i.ot, colorrange=extrema(df.ot),
		label=i.p,
		alpha=.5
	)
end
lines!(
	0:.1:3,
	z -> Dz_theoretical(z),
	label="model",
	color=:black
)

Legend(f[1,2], a, framevisible=false)
colgap!(f.layout, 0)
f

