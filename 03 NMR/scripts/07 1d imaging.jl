include("00 functions.jl")
using Glob, CSV
cd(@__DIR__)

# %%
data = DataFrame(
	paths = glob("../data/*1D*/*Image.txt"),
	direction = ["Z", "X", "Y"]
)

data[:, :images] =  [DataFrame(CSV.File(p), ["f", "S"]) for p in data.paths]
sort!(data, :direction)

data

# %%
f = Figure(size=halfsize)
axs = [Axis(f[i,1]) for i in 1:3]

for (d, a) in zip(eachrow(data), axs)
	lines!(a,
		d.images.f, d.images.S
	)
	a.yticks = [0]
	a.ylabel = d.direction
	ylims!(a, low=0)
end

for a in axs[1:end-1]
	linkaxes!(a, axs[end])
	hidexdecorations!(a, grid=false)
end

axs[end].xlabel = "Frequency shift in Hz"

rowgap!(f.layout, 0)
save("../figures/07 1d imaging.pdf", f)
f