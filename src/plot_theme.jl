using CairoMakie, Colors

# setup how figures are saved
figsize(width) = (width, width/1.618).*72
size_small = figsize(3)
size_large = figsize(5)
savefig(path::String, fig::Figure; kwargs...) = save(path, fig, pt_per_unit=1; kwargs...)
savelarge(path::String, fig::Figure; kwargs...) = save(path, fig, size=size_large, pt_per_unit=1; kwargs...)
savesmall(path::String, fig::Figure; kwargs...) = save(path, fig, size=size_small, pt_per_unit=1; kwargs...)

seeblau = colorant"#00A9E0"
seegrau = colorant"#9AA0A7"
seeblau_list = (colorant"#CCEEF9", colorant"#A6E1F4", colorant"#59C7EB", colorant"#00A9E0", colorant"#008ECE")
colors = [(0.0, 0.6627450980392157, 0.8784313725490196),
(0.6039215686274509, 0.6274509803921569, 0.6549019607843137),
(0.24313725490196078, 0.32941176470588235, 0.5882352941176471),
(0.5568627450980392, 0.12549019607843137, 0.2627450980392157),
(0.8784313725490196, 0.3764705882352941, 0.49411764705882355),
(0.996078431372549, 0.6274509803921569, 0.5647058823529412),
(0.0392156862745098, 0.5647058823529412, 0.5254901960784314),
(0.027450980392156862, 0.44313725490196076, 0.5294117647058824)]
colors = [RGB(t...) for t in colors]

theme = Theme(
	fontsize=10,
	size=size_large,
	fonts = (;
		regular="Roboto"
	),
	palette=(color=colors,),
	Axis = (
		leftspinecolor = seegrau,
		rightspinecolor = seegrau,
		bottomspinecolor = seegrau,
		topspinecolor = seegrau,
		xtickcolor=seegrau,
		ytickcolor=seegrau
	),
	Legend = (
		framecolor=seegrau,
	)
)
# set_theme!(theme)