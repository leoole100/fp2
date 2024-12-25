using FFTW, DataFrames, StatsBase
function spectrum(V::DataFrame)
	fs = 1/mean(diff(V.T))
	A = fft(V.V)
	freq = fftfreq(size(V, 1), fs)
	A = A[1:length(A)÷2]
	freq = freq[1:length(A)]
	A = A .* 2/fs
	return DataFrame([freq, abs.(A)], ["f", "A"])
end

# make plots consistent
using CairoMakie
inch = 96
pt = 4/3
cm = inch/2.54
update_theme!(
	fontsize=10,
	fonts = (;
		regular="Roboto Regular",
		bold="Roboto Bold",
	),
)
halfsize = (3inch, (3/1.618)inch)
fullsize = (5inch, (5/1.618)inch)