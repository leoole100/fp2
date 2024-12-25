using CSV, DataFrames, RollingFunctions
include("00 functions.jl")
cd(@__DIR__)

#%%
fid = DataFrame(CSV.File("../data/04 FID/t,S uV.txt"), ["t", "V"])
s = spectrum(fid.t, fid.V)

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="f in Hz", ylabel="FID in uV/Hz"
)

c = s.f[argmax(s.A)]
vlines!([c], color=:black)

lines!(s.f, s.A, label="Abs")
lines!(s.f, -s.Re, label="Re")
lines!(s.f, s.Im, label="Im")

# Legend(f[1,2], a, framevisible=false)
# colgap!(f.layout, 4)

axislegend(a, framevisible=false, position=:lt, padding=0)

w = 20
xlims!(c-w, c+w)

save("../figures/02 fid.pdf", f)
f