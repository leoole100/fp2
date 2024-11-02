# generate test images

using ImageDraw, Colors, Images, DataFrames
import VideoIO, CSV

#%%
# functions for generating test images
function draw_particle!(img::AbstractArray, center::Point, radius::Int, color::Colorant = Gray(0.1))
	draw!(img, Ellipse(CirclePointRadius(center, radius)), color)
end

function draw_scratch!(img::AbstractArray, p1::Point, p2::Point)
	draw!(img, LineSegment(p1, p2), Gray(0.0))
end

function test_img(
	particles,
	isize::Tuple{Int, Int} = (200, 200)
)

	# create a single test image
	image = fill(Gray(0.7), isize)

	# draw the particles
	for (center, radius) in particles
		draw_particle!(image, center, radius)
	end

	# add readout noise
	image += 0.01*randn(size(image)...)
	# add poison noise
	image = image .+ .2 .* image .* randn(size(image)...)

	draw_scratch!(image, Point(20, 20), Point(170, 170))

	# blur the image
	image = imfilter(image, Kernel.gaussian(1))

	return image
end

test_img([
	(Point(100, 100), 10),
	(Point(50, 100), 12)
])

#%%
# create a image sequence with particles moving
# make two particles move linearly

radii = [10, 12]
frames = 100
positions = Vector{Vector{Point}}(undef, frames)

for i in 1:frames
	positions[i] = [
		Point(100 + round(50*i/frames), 100),
		Point(50, 100 + round(50*i/frames))
	]
end

particles = [zip(positions[i], radii) for i in 1:frames]


# create the image sequence
images = [test_img(p) for p in particles]

mosaicview(first(images), last(images); ncol=2)

# %%
# save the images in a avi file

# convert images to 8-bit
images_8bit = [Float64.(clamp01.(img)) for img in images]
images_8bit = [round.(UInt8, 255*img) for img in images_8bit]

cd(@__DIR__)
VideoIO.save("../data/test/particles.avi", images_8bit)

# %%
# make a dataframe with positions for each frame and particles
particles = DataFrame()
for i in 1:length(radii)
	particles[:, Symbol(radii[i])] = [p[i] for p in positions]
end

particles
CSV.write("../data/test/particles.csv", particles)