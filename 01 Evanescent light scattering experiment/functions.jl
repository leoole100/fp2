function dist(data; cutoff=0, boundary=(minimum(data), maximum(data)))
	k = kde(data, boundary=boundary)
	x, y = k.x, k.density
	x, y = x[y .> cutoff], y[y .> cutoff]
	return (x=x, y=y)
end

potential(p) = -log.(p)

function remove_drift(t)
	t = copy(t)
	t = t .- t[1,:]'
	s = t[end,:] ./ (size(t,1)-1)
	for i in 0:(size(t,1)-1)
		t[i+1,:] = t[i+1,:] .- i .* s
	end
	return t
end

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