using CSV, DataFrames, RollingFunctions
include("00 functions.jl")
cd(@__DIR__)

#%%
fid = DataFrame(CSV.File("../data/04 FID/t,S uV.txt"), ["t", "V"])
s_fid = spectrum(fid.t, fid.V)
c = DataFrame(CSV.File("../data/05 C/t,S uV.txt"), ["t", "V"])
s_c = spectrum(c.t, c.V)

n = DataFrame(CSV.File("../data/01 Monitor_Noise/T,A in uV.txt"), ["t", "V"])
s_n = spectrum(n.t, n.V)

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="f in Hz", ylabel="FID in uV/Hz"
)


lines!(s_fid.f, s_fid.A, label="FID")
lines!(s_c.f, s_c.A, label="with decay")
# lines!(s_n.f, s_n.A, label="Noise")

axislegend(a, framevisible=false, position=:lt, padding=0)

c = s_fid.f[argmax(s_fid.A)]
w = 200
xlims!(c-w, c+w)

save("../figures/09 resonance.pdf", f)
f