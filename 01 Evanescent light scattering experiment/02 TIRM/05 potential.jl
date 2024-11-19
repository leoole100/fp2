using CairoMakie
using StatsBase: moment, std, mean, var
import StatsBase
using KernelDensity: kde
using Optim
using Format: format
using LsqFit: curve_fit, stderror
using Measurements: measurement

# load and filter data
include("functions.jl")

filter!(r->r.ot<1.1, df)

potential(p) = -log.(p) .+ log(maximum(p))

function dist(data; cutoff=0)
	k = kde(data, boundary=(minimum(data), maximum(data)))
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

model(x, p) = p[1] .* exp.(-x .* p[2]) .+ p[3] .* x .+ p[4]
fit_model(x, y) = curve_fit(model, x, y, [20.0, 20.0, 1.0, 0])
model_minimum(p) = 1/p[2] * log(p[2]p[1]/p[3])

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
df[:, :fit] .= fill(fill(measurement(0,0), 4), size(df, 1))
for (i,r) in enumerate(eachrow(df))
	I = r.I
	I = z_estimate(I, I0=1, β=r.β)
	k = dist(I, cutoff=0.05)
	v = potential(k.y)
	z = k.x
	m = fit_model(z, v)
	df[i, :fit] = measurement.(m.param, stderror(m))
	m = m.param
	mdl = model(z,m)

	z = z.-model_minimum(m)

	l = lines!(z, v, 
		label=format(r.ot, precision=2), 
		color=r.ot, colorrange=extrema(df.ot),
	)	
	lines!(z, mdl,
		color=r.ot, colorrange=extrema(df.ot), 
		linestyle=:dash
	)
end
Legend(f[1,2], a, "Trap Stiffness", framevisible=false)
colgap!(f.layout, 0)
save("../figures/02_05_01_potential.pdf", f)
f

# %%
s = [f[3]*kT*1e6*1e15 for f in df.fit] # fN
s = s .* measurement(1, .1)

k = [1/f[2]*1e3* measurement(1, .1) for f in df.fit] # fN

f = Figure(size=halfsize)
a = Axis(f[1,1], xlabel="Trap Stiffness", ylabel=" fN")
plot!(df.ot, value.(s))
errorbars!(df.ot, value.(s), uncertainty.(s))
line(x,p) = p[1].*x .+ p[2]
p = curve_fit(line, df.ot, value.(s), [1., 0.]).param
lines!(.618:.001:1, x->line(x, p), color=:black)
print(line(.618, p))
xlims!(a, low=.618)
save("../figures/02_05_02_gravity.pdf", f)
f
