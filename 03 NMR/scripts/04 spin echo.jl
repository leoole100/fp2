include("00 functions.jl")
using CSV, Glob, Format
cd(@__DIR__)

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

# %%
f = Figure(size=(halfsize[1],2.2inch))
axs = [Axis(f[i,1]) for i in 1:3]
axs_f = [Axis(f[i,2]) for i in 1:3]

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

save("../figures/04 spin echo shims combined.pdf", f)

f

# %%
f = Figure(size=halfsize)
a = Axis(f[1,1], ylabel="A in uV/Hz", xlabel="frequency in Hz")

for (p, s, t) in data
	sp = spectrum(t.t, t.S)
	c = sp.f[argmax(sp.A)]
	lines!(
		sp.f.-c, sp.A,
		color=s, colorrange=extrema(shims),
	)
	# a.ylabel = "shim:\n$(format(s))"
end

c = 1981.6
c=0
w = 4
xlims!(c-w, c+w)

save("../figures/04 spin echo shims spectrum.pdf", f)

f