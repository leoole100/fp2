using Images, ImageDraw

#%%

white = Gray(0.7)
black = Gray(0.1)

function draw_circle!(img, p, r, color::Colorant = Gray(0.0))
	draw!(img, Ellipse(CirclePointRadius(p..., r)), color)
end
function draw_line!(img, p1, p2)
	draw!(img, LineSegment(p1, p2), black)
end

#%%

img = fill(white, 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, (150, 150), 40)
# draw_circle!(img, (50, 50), 10)
draw_circle!(img, (100, 30), 50)
draw_line!(img, Point(20, 20), Point(170, 170))
img = imfilter(img, Kernel.gaussian(1))
save("data/test_images/image.png", img)
img