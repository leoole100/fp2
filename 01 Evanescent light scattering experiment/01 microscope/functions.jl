include("../functions.jl")

function load_trajectory!(df, path)
	path = join(split(path, ".")[1:end-1], ".")
	f = DataFrame(CSV.File(path * ".csv"))
	for i in 1:ncol(f)÷2
		push!(df[!, :t], hcat(f[:,2*i-1], f[:,2*i]))
		n = split(path, "/")[end]
		n = join(split(n, ".")[1:end-1], ".")
		n = split(n, " ")[3]
		ot = parse(Float64, n[3:end])
		push!(df[!, :n], n)
		push!(df[!, :ot], ot)
	end
end

function load_trajectories(path = "../data/TLM/")
	df = DataFrame(t=[], n=[], ot=[])
	paths = glob(path * "*.trajectory.csv")
	for p in paths
		load_trajectory!(df, p)
	end
	sort!(df, :ot)
	return df
end