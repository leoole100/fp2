using AlgebraOfGraphics, CairoMakie, DataFrames
using Revise
using fp2

# %%
set_theme!(fp2.theme)
cd(@__DIR__)


fig = Figure()
ax = Axis(fig[1,1])
plot!(ax, rand(10), rand(10), label="1")
plot!(ax, rand(10), rand(10), label="2")
plot!(ax, rand(10), rand(10), label="3")
fig[1, 2] = Legend(fig, ax, framevisible = false)
fig

fp2.savesmall("fig/2024-10-18 points.pdf", fig)
fp2.savelarge("fig/2024-10-18 points large.pdf", fig)