# %%

execfile("00 functions.py")

start = 0
step_size = 0.01
N_t = 10
N = 500000

def step(x=0, step_size=step_size):
	x = x + np.random.choice([-1.,1.])*step_size
	if x < 0:
		return 0.
	elif x > 1:
			return 1.
	else:
		return x
	
# %%
# run simulation
data = np.zeros((N_t, N))
for i in range(0,N_t):
	x = 0.5
	for j in range(0, N):
		data[i, j] = x
		x = step(x)

# %%
for d in data:
	plt.plot(d)

plt.show()

# %%
hist, bins = np.histogram(data.ravel(), bins=50, density=True)
plt.bar(bins[:-1], hist, align="edge", width=1/len(hist))
plt.plot([0,1], [1,1], color="black")
plt.xlabel("x")
plt.ylabel("Density")
plt.savefig("../figures/02 potential well.pdf")
plt.show()

# %%
