using Images, ImageDraw, ImageFeatures
using Revise
includet("functions.jl")

#%%
# generate a test image
img = fill(Gray(0.7), 200, 200)
img .= clamp01.(img .+ 0.1*randn(200, 200))
draw_circle!(img, Point(150, 150), 40)
draw_circle!(img, Point(100, 30), 50)
draw_line!(img, Point(20, 20), Point(170, 170))

img

#%%
# find circles in the image
centers, radii = find_circles(img, 40:60)
