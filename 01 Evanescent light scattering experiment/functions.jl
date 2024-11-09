function dist(data; cutoff=0)
	k = kde(data, boundary=(minimum(data), maximum(data)))
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

potential(p) = -log.(p)