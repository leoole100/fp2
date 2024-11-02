# find particles in images and track them

using Images, ImageFeatures, DataFrames, Statistics
using WGLMakie
import CSV, VideoIO

# %%
# open video file
cd(@__DIR__)
function load_image_gray(path)
	images = VideoIO.load(path)
	images = map(x -> Gray.(x), images)
	return images
end

images = load_image_gray("../data/test/particles.avi")
mosaicview(first(images), mean(images), last(images), nrow=1)

# %%
# normalize images

# get a reference white image
# white = mean(VideoIO.openvideo("../data/test/white.avi"))
# white = Gray.(white)
# black = mean(VideoIO.openvideo("../data/test/black.avi"))
# black = Gray.(black)

# generate one using the maximum intensity
white = fill(maximum(first(images)), size(first(images)))
black = fill(minimum(first(images)), size(first(images)))

# normalize images
normalize(img) = (Gray{Float64}.(img) .- black) ./ (white .- black)

custom_clamp(x) = x<.5 ? 0. : 1
preprocess(img) = Gray.(custom_clamp.((1 .- normalize(img)).^2))
preprocess(img)

images = map(preprocess, images)



# %%
# find the particles

img = first(images);
img_blurred = imfilter(img, Kernel.gaussian(5))
minima = findlocalmaxima(img_blurred)
particles = DataFrame(
	c = [(m[1], m[2]) for m in minima],
	a = Float64.(img_blurred[minima])
)
sort!(particles, :a, rev=true)
particles = particles[1:2, :]

particles.size = fill(10, nrow(particles))
particles

image(img; axis=(
	aspect = DataAspect(), yreversed = true
))
for p in eachrow(particles)
	scatter!(p.c..., color=:red)
	poly!(
		Point2f[
			(p.c[1] - p.size, p.c[2] - p.size),
			(p.c[1] + p.size, p.c[2] - p.size),
			(p.c[1] + p.size, p.c[2] + p.size),
			(p.c[1] - p.size, p.c[2] + p.size),
		],
		strokecolor = :red, strokewidth = 1, color = :transparent
	)
end
current_figure()

# %%
# track the particles
positions = fill((0., 0.), (length(images), nrow(particles)))

x_coords = [x for y in 1:size(img, 1), x in 1:size(img, 2)]
y_coords = [y for y in 1:size(img, 1), x in 1:size(img, 2)]

function crop_image(
	img::AbstractArray, 
	center::Tuple{Integer, Integer},
	radius::Integer, 
	padding::Integer
)
	y,x = center
	r = radius
	x1, x2 = max(1, x-r-padding), min(size(img, 2), x+r+padding)
	y1, y2 = max(1, y-r-padding), min(size(img, 1), y+r+padding)
	return @view img[y1:y2, x1:x2]
end

function center_of_mass(
	img_cropped::AbstractArray,
	c::Tuple{Int, Int},
	r::Int,
	padding::Int=0
)
	x_cropped = crop_image(x_coords, c, r, padding)
	y_cropped = crop_image(y_coords, c, r, padding)

	# find the center of mass
	tot = sum(img_cropped)
	x = sum(x_cropped .* img_cropped) / tot
	y = sum(y_cropped .* img_cropped) / tot

	x, y = Float64(x), Float64(y)

	return x, y
end

# allocate memory
processed = preprocess(first(images))


for i in 1:length(images)
	preprocessed = preprocess(images[i])
	for j in 1:nrow(particles)
		p = particles[j, :]
		img_cropped = crop_image(preprocessed, p.c, p.size, 0)
		c = center_of_mass(img_cropped, p.c, p.size, 0)
		positions[i, j] = c
		particles[j, :c] = round.(c)
	end

end


#%%
# plot the tracking results
image(
	mean(images); axis=(
	aspect = DataAspect(), yreversed = true
))
lines!(
	positions[:, 1],
)