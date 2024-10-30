using Images, ImageDraw, ImageFeatures, DataFrames
import CairoMakie
using Revise
includet("functions.jl")

#%%
# generate a test image
img = fill(Gray(0.7), 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, Point(150, 150), 30)
draw_circle!(img, Point(100, 30), 20)
draw_line!(img, Point(20, 20), Point(170, 170))
img = imfilter(img, Kernel.gaussian(1))
img .= clamp01.(img .+ 0.1*randn(200, 200))

img

#%%
# find circles in the image using the Hough transform
centers, radii = find_circles(img, 7:60;
	min_dist=20,
	vote_threshold=35,
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
# particles.img = Vector(undef, nrow(particles))

img_processed = copy(img)
img_processed = img_processed .- minimum(img_processed)
img_processed = img_processed ./ maximum(img_processed)
img_processed = imfilter(img_processed, Kernel.gaussian(3))
img_processed = (1 .- img_processed).^4
img_processed

for p in eachrow(particles)
	img_cropped = crop_image(img_processed, p.c, p.r, padding)
	c = center_of_mass(img_cropped, x_coords, y_coords, p.c, p.r, padding)
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
# using Laplacian of Gaussian blob detection
blobs = blob_LoG((1 .- img).^2, range(10,20, length=10))
blobs[1]
blobs = DataFrame(
	c = [(blob.location[2], blob.location[1]) for blob in blobs],
	σ = [blob.σ[1] for blob in blobs],
	a = [blob.amplitude for blob in blobs],
)
sort!(blobs, :a, rev=true)
blobs = blobs[ blobs.a .> 0.03, :]
blobs.r = blobs.σ .* sqrt(2)
blobs


# %%
# plot the results
f = CairoMakie.Figure()

CairoMakie.image(f[1, 1], img,
	axis = (aspect = CairoMakie.DataAspect(), yreversed = true,)
)

CairoMakie.scatter!(f[1, 1], last.(particles.c), first.(particles.c), color = :red, label = "Hough transform")
CairoMakie.scatter!(f[1, 1], last.(particles.c_opt), first.(particles.c_opt), color = :yellow, label = "Center of mass")
CairoMakie.scatter!(f[1, 1], last.(minima.c), first.(minima.c), color = :green, label = "Local minima")
CairoMakie.scatter!(f[1, 1], last.(blobs.c), first.(blobs.c), color = :blue, label = "LoG")

CairoMakie.axislegend()

cd(@__DIR__)
save("../figures/single frame.pdf", f)
f