using FFTW, DataFrames, StatsBase, Glob
function spectrum(t, V)
	fs = 1/mean(diff(t))
	A = fft(V)
	freq = fftfreq(size(V)[1], fs)
	A = A[1:length(A)÷2]
	freq = freq[1:length(A)]
	A = A .* 2/fs
	# return DataFrame([freq, abs.(A)], ["f", "A"])
	return DataFrame(Dict(
		:f=>freq,
		:A=>abs.(A),
		:Re=>real.(A),
		:Im=>imag.(A)
	))
end

function runall()
	cd(@__DIR__)
	for p in glob("*.jl")
		println("\n"*p)
		include(p)
	end
end

# make plots consistent
using CairoMakie
inch = 96
pt = 4/3
cm = inch/2.54
update_theme!(
	fontsize=10pt,
	fonts = (;
		regular="Roboto Regular",
		bold="Roboto Bold",
	),
)
halfsize = (3inch, (3/1.618)inch)
fullsize = (5inch, (5/1.618)inch)