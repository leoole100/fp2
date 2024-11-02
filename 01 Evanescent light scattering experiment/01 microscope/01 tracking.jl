# find particles in images and track them

using Images, ImageFeatures, DataFrames, Statistics
using WGLMakie
import CSV, VideoIO

path = "../data/test/particles.avi"

# %%
# open video file
cd(@__DIR__)
video = VideoIO.openvideo(path)
img = read(video)
# image(img',axis=(aspect=DataAspect(),yreversed = true,))

#%%
# make normalization images
white = fill(maximum(Gray{Float64}.(img)), size(img))
black = fill(minimum(Gray{Float64}.(img)), size(img))

normalize(img) = (clamp01.(Gray.(img)) .- black) ./ (white .- black)
# normalize(img)

# %%
# transform the images
trans(x) = x^4
transform(img) = Gray.(trans.(1 .- normalize(img)))
# transform(img)

# f = Figure()
# image(f[1:2,1], transform(img)', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false);
# lines(f[1,2], 0:0.05:1, trans);
# hist(f[2,2], reduce(vcat, Float64.(transform(img))), normalization=:probability, axis=(yscale=log10,));
# f

# %%
# find particles in the images
# using Laplacian of Gaussian blob detection
blobs = blob_LoG(transform(img),  range(2,10, length=10));
blobs = DataFrame(
	c = [(blob.location[2], blob.location[1]) for blob in blobs],
	σ = [blob.σ[1] for blob in blobs],
	a = [blob.amplitude for blob in blobs],
);
sort!(blobs, :a, rev=true)
# blobs = blobs[ blobs.a .> 0.03, :]
blobs = blobs[1:2, :]
blobs.r = blobs.σ .* 2;
# blobs

# image(img', axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)
# scatter!(blobs.c, markersize=blobs.r.*2, color=:red, marker=:x, markerspace=:data)
# current_figure()

# %%
# last positions
last_positions = hcat([Float64.([first(i), last(i)]) for i in blobs.c]...)'
radii = blobs.r

coords = [Float64.([x, y ]) for y in 1:size(img, 1), x in 1:size(img, 2)];
com(i, c) = sum(c .* Float64.(i))/sum(Float64.(i)) # center of mass
function crop(img, c, r)
	x, y = round.(Int, c)
	r = round(Int, r)
	x1, x2 = max(1, x-r), min(size(img, 2), x+r)
	y1, y2 = max(1, y-r), min(size(img, 1), y+r)
	return @view img[y1:y2, x1:x2]	
end
# crop(transform(img), last_positions[1, :], radii[1])
cropped_com(img, c, r) = com(crop(img, c, r), crop(coords, c, r))
# cropped_com(transform(img), last_positions[1, :], radii[1])

function update_positions!(last_positions, radii, img)
	for i in 1:length(radii)
		last_positions[i, :] = cropped_com(img, last_positions[i, :], radii[i] .+ 5)
	end
end

function update_positions(last_positions, radii, img)
	new_positions = copy(last_positions)
	update_positions!(new_positions, radii, img)
	return new_positions	
end

# @time update_positions(last_positions, radii, transform(img)) # <1 ms

# %%
# track the particles
last_positions = transpose(hcat([Float64.([first(i), last(i)]) for i in blobs.c]...))
last_positions = reinterpret(Float64, last_positions)
positions = Array{Float64, 2}[last_positions]

video = VideoIO.openvideo(path)
img = read(video)
transformed = transform(img)
meanImg = Gray{Float64}.(transformed)
while !eof(video)
	img = read(video)
	transformed = transform(img)
	update_positions!(last_positions, radii, transformed)
	push!(positions, last_positions)
	meanImg = meanImg .+ img
end

meanImg = meanImg ./ length(positions)
positions

# %%
df = DataFrame()

# add a column for each particle
for i in 1:2
	df[!, "x$i"] = [pos[i, 1] for pos in positions]
	df[!, "y$i"] = [pos[i, 2] for pos in positions]
end

df

basepath = join(split(path, '.')[1:end-1], ".")
CSV.write(basepath * ".trajectory.csv", df)
save(basepath * ".mean.png", clamp01.(Gray.(meanImg)))
for i in 1:length(radii)
	save(basepath * ".particle$i.png", crop(img, last_positions[i, :], radii[i] .+ 5))
end