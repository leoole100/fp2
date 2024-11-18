using CairoMakie
using StatsBase: moment, std, autocor
using KernelDensity: kde
using Optim: optimize, GradientDescent
using DSP: conv

# load and filter data
include("functions.jl")

D_estimate(std_z; Δt=Δt) = std_z.^2 ./ (2*Δt)  # from script

# %% look at same z
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

