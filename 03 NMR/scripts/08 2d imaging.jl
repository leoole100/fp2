include("00 functions.jl")
using Glob, CSV
cd(@__DIR__)

# %%
data = DataFrame(
	path = glob("../data/*2D*/*.txt")[2:end],
	shape=[ # in cm
		[12, 12],
		[10, 10],
		[10, 10],
		[10, 10]
	],
	duration = [4, .64, 4, 4], 	# Polarizing duration
	echo = [.2, .2, .5, .3]		# echo time
)
data = data[[2,1,3], :]


data[:, :image] =  [CSV.read(p, CSV.Tables.matrix; header=false, delim=" ") for p in data.path]./1000
data.image = transpose.(data.image)
ext = extrema(vcat(hcat(data.image ...)...))

data



# %%
f = Figure(size=(6inch, fullsize[2]))
axs = [Axis(f[1,i]) for i in 1:size(data,1)]

for (i, (d, a)) in enumerate(zip(eachrow(data), axs))
	local x = (-d.shape[1]/2, d.shape[1]/2)
	local y = (-d.shape[2]/2, d.shape[2]/2)
	
	im = image!(a,
		x, y, d.image,
		# colorrange=ext,
		colormap=[:transparent, :black]
	)

	contour!(a, x, y, d.image)

	for p in [[1.7, -1.9], [.9, 2]]
		arc!(a, p, 1, 0, 2π, color=:red, linestyle=:dash)
	end

	# Colorbar(f[2,i], im, vertical=false, flipaxis=false)
	linkaxes!(a, axs[1])
	a.aspect = DataAspect()
	a.title = "Polarizing $(d.duration) s, Echo $(d.echo) s"

	xlims!(-2, 4)
	ylims!(-4, 4)
end

axs[1].ylabel = "Z in cm"
Label(f[2,:], "y in cm")
rowgap!(f.layout, 5)


for a in axs[2:end]
	hideydecorations!(a, grid=false)
end

resize_to_layout!(f)
save("../figures/08 2d imaging.pdf", f)
f