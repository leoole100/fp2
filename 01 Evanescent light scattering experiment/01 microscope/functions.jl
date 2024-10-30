using Images, ImageDraw, ImageFeatures

function draw_circle!(img::AbstractArray, center::Point, radius::Int, color::Colorant = Gray(0.0))
	draw!(img, Ellipse(CirclePointRadius(center, radius)), color)
end

function draw_line!(img::AbstractArray, p1::Point, p2::Point)
	draw!(img, LineSegment(p1, p2), Gray(0.1))
end

function find_circles(
	img::AbstractArray, 
	radii::UnitRange, 
	min_dist::Number=minimum(radii), 
	vote_threshold::Number=minimum(radii);
	canny_threshold::Tuple{Number, Number}=(Percentile(99), Percentile(80))
)
	img = Gray.(img)
	img_edges = canny(img, canny_threshold)
	dx, dy = imgradients(img, KernelFactors.ando5)
	img_phase = phase(dx, dy)
	centers, radii = hough_circle_gradient(img_edges, img_phase, radii, min_dist = min_dist, vote_threshold=vote_threshold)
	return centers, radii
end