using CSV, DataFrames
include("00 functions.jl")
cd(@__DIR__)

#%%
v = DataFrame(CSV.File("../data/01 Monitor_Noise/T,A in uV.txt"), ["T", "V"])
s = spectrum(v)

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="f in Hz", ylabel="Amplitude in uV/Hz"
)
vlines!(1977.6, color=:black)
lines!(s.f, s.A)
xlims!(1500, 2500)
# ylims!(low=0)
save("../figures/01 noise.pdf", f)
f