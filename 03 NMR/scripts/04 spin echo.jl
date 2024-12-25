include("00 functions.jl")
using CSV, Glob, Format

# %%

paths = sort(glob("../data/*pin echo*/t*.txt"))
shims = [0, 0.59, 0.295]

traces = [DataFrame(CSV.File(p), ["t", "S"]) for p in paths]

data = collect(zip(paths, shims, traces))

sort!(data; by= i->i[2])

# %%
f = Figure(size=(halfsize[1],2.2inch))
axs = [Axis(f[i,1]) for i in 1:3]

for a in axs[1:end-1]
	hidexdecorations!(a, grid=false)
end

for ((p, s, t), a) in zip(data, axs)
	lines!(a,
		t.t, t.S,
		color=s, colorrange=extrema(shims),
	)
	linkyaxes!(a, axs[1])
	# a.ylabel = "shim:\n$(format(s))"
end

rowgap!(f.layout, 0)
axs[end].xlabel = "Time in S"
Label(f[:, 0], "Signal in uV", rotation = pi/2)

save("../figures/04 spin echo shims.pdf", f)

f