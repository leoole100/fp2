using Images, ImageDraw, ImageFeatures, DataFrames, CairoMakie
using Revise
includet("functions.jl")

#%%
# generate a test image
img = fill(Gray(0.7), 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, Point(150, 150), 30)
draw_circle!(img, Point(100, 30), 20)
draw_line!(img, Point(20, 20), Point(170, 170))
img .= clamp01.(img .+ 0.1*randn(200, 200))

img

#%%
# find circles in the image using the Hough transform
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
# optimize using the center of mass
padding = 1

x_coords = [x for y in 1:size(img, 1), x in 1:size(img, 2)]
y_coords = [y for y in 1:size(img, 1), x in 1:size(img, 2)]

# add columns to the particles DataFrame
particles.c_opt = fill((0., 0.), nrow(particles))
particles.img = Vector(undef, nrow(particles))

for p in eachrow(particles)
	img_cropped = crop_image(img, p.c, p.r, padding)
	c = center_of_mass((1 .- img_cropped).^2, x_coords, y_coords, p.c, p.r, padding)
	p.img = img_cropped
	p.c_opt = c
end

particles


# %%
# using local minimum search
# blur the image
img_blurred = imfilter(img, Kernel.gaussian(20))
1 .- img_blurred
minima = findlocalminima(img_blurred)
minima = DataFrame(c = [(c[2], c[1]) for c in minima])
minima.strength = [Float64(img_blurred[c[2], c[1]]) for c in minima.c]
sort!(minima, :strength, rev=false)
minima = minima[ minima.strength .< 0.5, :]

minima

# %%
# plot the results
f = Figure()

image(f[1, 1], img,
	axis = (aspect = DataAspect(), yreversed = true,)
)

scatter!(f[1, 1], last.(particles.c), first.(particles.c), color = :red, label = "Hough transform")
scatter!(f[1, 1], last.(particles.c_opt), first.(particles.c_opt), color = :yellow, label = "Center of mass")
scatter!(f[1, 1], last.(minima.c), first.(minima.c), color = :green, label = "Local minima")

axislegend()

f