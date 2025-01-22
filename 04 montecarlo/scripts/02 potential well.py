# %%

execfile("00 functions.py")

start = 0
step_size = 0.01
N_t = 10
N = 500000

def step(x=0, step_size=step_size):
	x_new = x + np.random.choice([-1,1])*step_size
	if x_new < 0 or x_new > 1:
		return x
	else:
		return x_new
	
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
hist, bins = np.histogram(data, bins=100, density=True, range=(0,1))
plt.bar(bins[:-1], hist)
plt.xlabel("x")
plt.ylabel("Density")
plt.savefig("../figures/02 potential well.pdf")
plt.show()
