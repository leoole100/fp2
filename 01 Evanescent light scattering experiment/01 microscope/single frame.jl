using Images, ImageDraw, ImageFeatures, DataFrames, CairoMakie
using Revise
includet("functions.jl")

#%%
# generate a test image
img = fill(Gray(0.7), 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, Point(150, 150), 40)
draw_circle!(img, Point(100, 30), 10)
draw_line!(img, Point(20, 20), Point(170, 170))

img

#%%
# find circles in the image
centers, radii = find_circles(img, 7:60;
	min_dist=20,
	vote_threshold=25,
)

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
particles.img = Vector(undef, nrow(particles))

for p in eachrow(particles)
	img_cropped = crop_image(img, p.c, p.r, padding)
	c = center_of_mass((1 .- img_cropped), x_coords, y_coords, p.c, p.r, padding)
	p.img = img_cropped
	p.c_opt = c
end

particles

#%%
f = Figure()

image(f[1, 1], img,
	axis = (aspect = DataAspect(), yreversed = true,)
)

scatter!(f[1, 1], last.(particles.c), first.(particles.c), color = :red)
scatter!(f[1, 1], last.(particles.c_opt), first.(particles.c_opt), color = :yellow)

f