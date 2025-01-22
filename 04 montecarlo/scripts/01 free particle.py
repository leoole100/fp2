# %%

execfile("00 functions.py")

start = 0
step_size = 0.01
N_t = 50
N = 1000

def step(x=0, step_size=step_size):
	return x + np.random.choice([-1,1])*step_size


# %%
# run simulation
data = np.zeros((N_t, N))
for i in range(0,N_t):
	x = 0
	for j in range(0, N):
		data[i, j] = x
		x = step(x)

# %%
for d in data:
	plt.plot(d)

plt.show()

# %%
# msd
msd = np.cumsum(data**2, axis=1)
msd_mean = np.average(msd, axis=0)

plt.plot(msd_mean)
plt.xlabel("Step")
plt.ylabel("MSD")
plt.xscale("log")
plt.yscale("log")
plt.savefig("../figures/01 free particle.pdf")
plt.show()

# %%
