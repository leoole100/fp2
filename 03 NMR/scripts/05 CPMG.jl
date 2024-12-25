include("00 functions.jl")
using CSV, Glob, Format, LsqFit, Measurements

# %%

paths = sort(glob("../data/*CPMG*/"))[5:end]
filter!(p->!occursin("contrast", p), paths)

traces = [DataFrame(CSV.File(p*"t ms, S uV.txt"), ["t", "S"]) for p in paths]
attenuation = [DataFrame(CSV.File(p*"t ms, Attenuation.txt"), ["t", "A"]) for p in paths]

data = collect(zip(paths, traces, attenuation))

# %%
mdl(x, p) = exp.(-1 .* x ./ p[1])

fits = [
	curve_fit(mdl, a.t, a.A, [2000.])
	for (p, t, a) in data
]

T2 = [measurement(f.param[1], stderror(f)[1]) for f in fits]


f = Figure(size=(7inch, fullsize[2]))

axs = [Axis(f[i,1]) for i in 1:length(data)]

for ((p, t, a), ax, c) in zip(data, axs, Makie.wong_colors())
	lines!(ax,
		t.t, t.S,
		color=c
	)
	linkxaxes!(ax, axs[1])
	linkyaxes!(ax, axs[1])
end

for a in axs[1:end-1]
	hidexdecorations!(a, grid=false)
end
rowgap!(f.layout, 0)
axs[end].xlabel = "Time in ms"
Label(f[:, 0], "Signal in uV", rotation = pi/2)


a = Axis(f[1:length(data),2],
		xlabel="Time in ms", ylabel="Attenuation"
)

for ((p, t, a), f, t2) in zip(data, fits, T2)
	local s = scatter!(a.t, a.A, label=split(p, "/")[end-1]*"\nT2=$(t2) ms")
	local x = 1:7000
	lines!(x, mdl(x, f.param), color=s.color)
end

xlims!(low=0)
ylims!(0, 1)
# axislegend()
Legend(f[:,3], a, framevisible=false)

save("../figures/05 CPMG.pdf", f)
f