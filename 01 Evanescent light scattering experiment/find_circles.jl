using ImageFeatures, Images, DataFrames, ImageDraw

#%%
img = load("data/test_images/image.png")
img = Gray.(img)

img_edges = canny(img, (Percentile(99), Percentile(80)));
# image(img_edges)
dx, dy = imgradients(img, KernelFactors.ando5);
img_phase = phase(dx, dy);

centers, radii = hough_circle_gradient(img_edges, img_phase, 40:60,
	min_dist = 50, vote_threshold=20
)
centers = Tuple.(centers)

#%%

# draw the detected circles
img_result = fill(RGB(0.0, 0.0, 0.0), size(img))
img_result += img
img_result

for (c, r) in zip(centers, radii)
	draw!(img_result, Ellipse(CirclePointRadius(c[2], c[1], 3)), RGB(1.0, 0.0, 0.0))
	draw!(img_result, Ellipse(CirclePointRadius(c[2], c[1], r; fill=false, thickness=20)), RGB(1.0, 0.0, 0.0))
end

save("data/test_images/image_result.png", img_result)
img_result