using CairoMakie
using Measurements: value, uncertainty, ±
using DataFrames
cd(@__DIR__)

#%%

d = DataFrame(
	c = [0, 0.25, 0.5, 1, 2],
	T2 = [2400±60, 1320±40, 950±50, 600±40, 330±40],
	T1 = [1800±100, 1200±020, 860±30, 590±40, 306±8],
)

f = Figure(size=(400,300))
a = Axis(f[1,1];
 xlabel="Concentration mMol", ylabel="Time Constant ms",
#  xscale=log10, yscale=log10
)

for (t,l) in zip([d.T1, d.T2], ["T1", "T2"])
	scatter!(d.c, value.(t)./1e3, label=l)
	errorbars!(d.c, value.(t)./1e3, uncertainty.(t)./1e3)
end
axislegend()
save("../figures/01 concentration.pdf")
f