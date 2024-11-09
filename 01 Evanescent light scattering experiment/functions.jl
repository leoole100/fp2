function dist(data; cutoff=0, boundary=(minimum(data), maximum(data)))
	k = kde(data, boundary=boundary)
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

potential(p) = -log.(p)

# make plots consistent
# will be moved to a top level file
inch = 96
pt = 4/3
cm = inch/2.54
update_theme!(
	# fontsize=10,
	# fonts = (;
	# 	regular="Roboto"
	# ),
)