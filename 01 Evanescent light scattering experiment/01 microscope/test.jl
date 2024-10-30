using Images, ImageDraw, ImageFeatures, DataFrames
using Revise
includet("functions.jl")

#%%
# generate a test image
img = fill(Gray(0.7), 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, Point(150, 150), 40)
# draw_circle!(img, Point(100, 30), 50)
draw_line!(img, Point(20, 20), Point(170, 170))

img

#%%
# find circles in the image
centers, radii = find_circles(img, 40:60)

particles = DataFrame(c=Tuple{Int, Int}[], r=Int[])
for (center, radius) in zip(centers, radii)
	push!(particles, (center, radius))
end
particles

# %%
# crop the image for the particles
padding = 1

x_coords = [x for y in 1:size(img, 1), x in 1:size(img, 2)]
y_coords = [y for y in 1:size(img, 1), x in 1:size(img, 2)]

# add columns to the particles DataFrame
particles.c_opt = fill((0., 0.), nrow(particles))

for p in eachrow(particles)
	img_cropped = crop_image(img, p.c, p.r, padding)
	c = center_of_mass(1 .- img_cropped, x_coords, y_coords, p.c, p.r, padding)
	p.c_opt = c
end

particles