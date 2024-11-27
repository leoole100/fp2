using DataFrames: DataFrame, groupby
using Measurements: measurement

function estimate_noise(df)
	# return minimum(df.V2)	
	return 0
end

function estimate_noise_VJ2!(df)
	grouped = groupby(df, :G1)
	for g in grouped
		noise = estimate_noise(g)
		g.noise .= noise
		g.VJ2 .= g.V2 .- noise
	end
end

function fit_groups(
	df::DataFrame,
	g::Symbol,
	x::Symbol,
	y::Symbol,
	mdl::Function,
	p0::Vector,
) 
	grouped = groupby(df, g)
	n = length(grouped)

	fit_params = DataFrame(Dict(
		g => unique(grouped.groups),
		:p => fill(measurement.(p0), n),
	))

	for (g, f) in zip(grouped, eachrow(fit_params))
		p = curve_fit(mdl, value.(g[!, x]), value.(g[!, y]), p0)
		f.p = measurement.(p.param, stderror(p))
	end

	return fit_params
end