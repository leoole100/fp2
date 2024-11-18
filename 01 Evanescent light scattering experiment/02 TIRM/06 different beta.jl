using CairoMakie
using StatsBase: moment, std, mean, var
import StatsBase
using KernelDensity: kde
using Optim
using Format: format
using LsqFit: curve_fit

# load and filter data
include("functions.jl")

filter!(r->r.ot==1.5, df)

potential(p) = -log.(p) .+ log(maximum(p))

function dist(data; cutoff=0)
	k = kde(data, boundary=(minimum(data), maximum(data)))
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

model(x, p) = p[1] .* exp.(-x .* p[2]) .+ p[3] .* x .+ p[4]
fit_model(x, y) = curve_fit(model, x, y, [20.0, 20.0, 1.0, 0]).param
model_minimum(p) = 1/p[2] * log(p[2]p[1]/p[3])

background(I) = -.1
estimate_I0(I) = 2*maximum(I)-background(I)
calc_z(I, I0=estimate_I0(I), beta=1, b=background(I)) = log.(I0 ./ (I.-b)) ./ beta

# add beta from df.l column
l_beta = Dict(5.0=>.454, 10.0=>.322, 15.0=>.265, 20.0=>.230)
df[:, :β] = zeros(size(df, 1))
for i in 1:size(df, 1)
	df[i, :β] = l_beta[df[i, :l]+5.0]
end

# %%
f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="z in μm", ylabel="V(z) in kT",
	# yscale=log10,
	# xscale=log10,
)
df[:, :fit] .= fill(zeros(4), size(df, 1))
for (i,r) in enumerate(eachrow(df))
	I = r.I
	# I = z_estimate(I, I0=2*maximum(I), β=r.β)
	I = calc_z(I)
	k = dist(I, cutoff=0.05)
	v = potential(k.y)
	z = k.x
	m = fit_model(z, v)
	df[i, :fit] = m
	mdl = model(z,m)

	z = z./model_minimum(m) .-1

	l = lines!(z, v, label=format(r.β, precision=3))	
	lines!(z, mdl, linestyle=:dash, color=l.color)
end
Legend(f[1,2], a, "β in μm",framevisible=false)
colgap!(f.layout, 0)
save("../figures/02_06_01_different_beta.pdf", f)
f

