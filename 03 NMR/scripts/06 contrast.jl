include("00 functions.jl")
using Measurements, LsqFit
cd(@__DIR__)

#%%
d = DataFrame(
	c = [0, 0.25, 0.5, 1, 2],
	T1 = [1800±100, 1200±020, 860±30, 590±40, 306±8],
	T2 = [2400±60, 1320±40, 950±50, 600±40, 330±40],
)
d.T1 /= 1000
d.T2 /= 1000

d = d[2:end,:]



f = Figure(size=halfsize)
a = Axis(f[1,1];
 xlabel="Concentration mMol", ylabel="Time Constant s",
 xscale=log10, yscale=log10
)

for (t,l) in zip([d.T1, d.T2], ["T1", "T2"])
	scatter!(d.c, value.(t), label=l)
	errorbars!(d.c, value.(t), uncertainty.(t))
end
axislegend()
save("../figures/06 contrast.pdf", f)
f