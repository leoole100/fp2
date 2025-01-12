using CSV, DataFrames, RollingFunctions
include("00 functions.jl")
cd(@__DIR__)

#%%
c = DataFrame(CSV.File("../data/02 Analyse_Coil/C,F.txt"), ["C", "f"])
v = DataFrame(CSV.File("../data/01 Monitor_Noise/T,A in uV.txt"), ["T", "V"])
println("RMS: ", sqrt(mean(v.V .^2)))

s = spectrum(v.T, v.V)
fid = DataFrame(CSV.File("../data/04 FID/f,A.txt"), ["f", "A"])

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="f in Hz", ylabel="Noise in uV/Hz"
)
hideydecorations!(a; label=false, ticklabels=false, ticks=false)
ylims!(0, .5)

# lines!(rollmean(fid.f,20), rollmean(fid.A, 20), color=:gray)
lines!(s.f, s.A)


# plot capacitance
color = Makie.wong_colors()[2]
aC = Axis(f[1,1],
	yaxisposition=:right, yticklabelcolor = color,
	ylabel="C in nF"
)
hidexdecorations!(aC)
hideydecorations!(aC; label=false, ticklabels=false, ticks=false)
linkxaxes!(a, aC)

lines!(c.f, c.C, color=color)
scatter!([1982], [11.2], color=color)

xlims!(1000, 3000)
save("../figures/01 noise.pdf", f)
f