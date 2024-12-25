include("00 functions.jl")
using CSV, LsqFit, Measurements

# %%

at = DataFrame(CSV.File("../data/06 T1Bp Water/Delay ms, Attenuation with Ebars.txt"), [:d, :A, :Amin, :Amax])
at.d = at.d/1000

# %%
mdl(x, p) = 1 .- exp.(-1 .* x ./ p[1])
fitp = curve_fit(mdl, at.d, at.A, [2.])
T1 = measurement(fitp.param[1], stderror(fitp)[1])
print(T1)

f = Figure(size=halfsize)
a = Axis(f[1,1],
	ylabel = "Attenuation", xlabel="τₚ in S"
)

scatter!(at.d, at.A)

x = 0:.1:4
lines!(x, mdl(x, fitp.param))

xlims!(low=0)
ylims!(0, 1)

save("../figures/03 T1.pdf", f)
f

