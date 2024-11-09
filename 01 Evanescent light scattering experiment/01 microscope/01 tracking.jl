# find particles in images and track them

using DataFrames: DataFrame, sort!
using Statistics: mean
using ImageFeatures: blob_LoG
using Colors: Gray
using Images: clamp01
using WGLMakie
import CSV, VideoIO
using FileIO: load

path = "../data/TLM/006 pt02 ot0.avi"

#%%
# open video file
cd(@__DIR__)
video = VideoIO.openvideo(path)
img = read(video)
# image(img',axis=(aspect=DataAspect(),yreversed = true,))

# make normalization images
white = load("../data/TLM/white.bmp")
black = load("../data/TLM/black.bmp")
white, black = float.(Gray.(white)), float.(Gray.(black))

normalize(img) = (clamp01.(Gray.(float(img))) .- black) ./ (white .- black)
normalize(img)

# transform the images
trans(x) = x^4
transform(img) = Gray.(trans.(1 .- normalize(img)))
transform(img)

image(transform(img)'./maximum(transform(img)), axis=(aspect=DataAspect(), yreversed=true,), interpolate=false)


# %%

coords = [Float64.([x, y ]) for y in 1:size(img, 1), x in 1:size(img, 2)];
com(i, c) = sum(c .* Float64.(i))/sum(Float64.(i)) # center of mass

# find single particle
last_positions = hcat(com(transform(img), coords))'
radii = [40]

function crop(img, c, r)
	x, y = round.(Int, c)
	r = round(Int, r)
	x1, x2 = max(1, x-r), min(size(img, 2), x+r)
	y1, y2 = max(1, y-r), min(size(img, 1), y+r)
	return @view img[y1:y2, x1:x2]	
end
cropped_com(img, c, r) = com(crop(img, c, r), crop(coords, c, r))
# cropped_com(transform(img), last_positions[1, :], radii[1])

# @time com(transform(img), coords)

function update_positions!(last_positions, radii, img)
	for i in 1:length(radii)
		last_positions[i, :] = cropped_com(img, last_positions[i, :], radii[i] .+ 5)
	end
end

update_positions!(last_positions, radii, transform(img))
@time update_positions!(last_positions, radii, transform(img)) # <1 ms

crop(normalize(img), last_positions[1, :], radii[1]) # show the particle

# %%
# track the particles
# last_positions = transpose(hcat([Float64.([first(i), last(i)]) for i in blobs.c]...))
last_positions = hcat(com(transform(img), coords))'
last_positions = reinterpret(Float64, last_positions)
positions = Array{Float64, 2}[last_positions]

video = VideoIO.openvideo(path)
img = read(video)
transformed = transform(img)
meanImg = Gray{Float64}.(transformed)
i = 0
while !eof(video)
	img = read(video)
	transformed = transform(img)
	update_positions!(last_positions, radii, transformed)
	push!(positions, last_positions)
	# meanImg = normalize(meanImg) .+ img
	i += 1
	if i % 100 == 0
		println(i)
	end
end

# meanImg = meanImg ./ length(positions)
positions

# %%
df = DataFrame()

# add a column for each particle
df[!, "x"] = [p[1] for p in positions]
df[!, "y"] = [p[2] for p in positions]

df

basepath = join(split(path, '.')[1:end-1], ".")
CSV.write(basepath * ".trajectory.csv", df)
save(basepath * ".png", clamp01.(Gray.(img)))
save(basepath * ".particle.png", clamp01.(crop(normalize(img), last_positions[1, :], radii[1])))