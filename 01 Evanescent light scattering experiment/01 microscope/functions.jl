using Images, ImageDraw, ImageFeatures

function draw_circle!(img::AbstractArray, center::Point, radius::Int, color::Colorant = Gray(0.0))
	draw!(img, Ellipse(CirclePointRadius(center, radius)), color)
end

function draw_line!(img::AbstractArray, p1::Point, p2::Point)
	draw!(img, LineSegment(p1, p2), Gray(0.1))
end

function test_image(points::Vector{Point}, radii::Vector{Int64})
	# make the background
	img = fill(Gray(0.7), 200, 200)
	img .= clamp01.(img .+ 0.1*randn(200, 200))

	# draw the circles
	for (center, radius) in zip(points, radii)
		draw_circle!(img, center, radius)
	end

	# draw scratches
	draw_line!(img, Point(20, 20), Point(170, 170))

	# add noise
	img .= clamp01.(img .+ 0.1*randn(200, 200))
	img = imfilter(img, Kernel.gaussian(1))

	return img
end


"""
returns the centers and radii of circles found in the image, in integer pixel precision
"""
function find_circles(
	img::AbstractArray, 
	radii::UnitRange;
	min_dist::Number=minimum(radii), 
	vote_threshold::Number=minimum(radii),
	canny_threshold::Tuple{Number, Number}=(Percentile(99), Percentile(80))
)
	img = Gray.(img)
	img_edges = canny(img, canny_threshold)
	dx, dy = imgradients(img, KernelFactors.ando5)
	img_phase = phase(dx, dy)
	centers, radii = hough_circle_gradient(img_edges, img_phase, radii, min_dist = min_dist, vote_threshold=vote_threshold)
	centers = [(center[2], center[1]) for center in centers]
	return centers, radii
end

function crop_image(
	img::AbstractArray, 
	center::Tuple{Integer, Integer},
	radius::Integer, 
	padding::Integer
)
	x, y = center
	r = radius
	x1, x2 = max(1, x-r-padding), min(size(img, 2), x+r+padding)
	y1, y2 = max(1, y-r-padding), min(size(img, 1), y+r+padding)
	return @view img[y1:y2, x1:x2]
end

function center_of_mass(
	img_cropped::AbstractArray,
	x_coords::AbstractArray{Int, 2},
	y_coords::AbstractArray{Int, 2},
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