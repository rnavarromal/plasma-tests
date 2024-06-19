import numpy as np
import matplotlib.pyplot as plt

nsteps = 64000
dt = 0.01

x = np.empty(shape=(nsteps,3))
v = np.empty(shape=(nsteps,3))

dtype=np.float64

with open("output.dat", "rb") as f:
    for i in range(nsteps):
        x[i] = np.fromfile(f, dtype=dtype, count=3)
        v[i] = np.fromfile(f, dtype=dtype, count=3)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*x.T)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
# plt.gca().set_aspect("equal")
plt.show()

