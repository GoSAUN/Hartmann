import numpy as np
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate

x1 = np.loadtxt("x1.dat")
x2 = np.loadtxt("x2.dat")
x3 = np.loadtxt("x3.dat")
x4 = np.loadtxt("x4.dat")
x5 = np.loadtxt("x5.dat")

x1T = integrate.simps(x1)


fig, axs = plt.subplots(2, 2, figsize=(8, 8))
fig.suptitle("Finite Mode", fontsize=14)

axs[0, 0].plot(x1,x2,"b-")
axs[0, 0].set_xlabel("$x_{1}$")
axs[0, 0].set_ylabel("$x_{2}$")
axs[0, 0].grid(True)
axs[0, 0].set(aspect='equal')
axs[1, 0].plot(x1,x3,"b-")
axs[1, 0].set_xlabel("$x_{1}$")
axs[1, 0].set_ylabel("$x_{3}$")
axs[1, 0].grid(True)
axs[1, 0].set(aspect='equal')
axs[0, 1].plot(x1,x4,"b-")
axs[0, 1].set_xlabel("$x_{1}$")
axs[0, 1].set_ylabel("$x_{4}$")
axs[0, 1].grid(True)
axs[0, 1].set(aspect='equal')
axs[1, 1].plot(x1,x5,"b-")
axs[1, 1].set_xlabel("$x_{1}$")
axs[1, 1].set_ylabel("$x_{5}$")
axs[1, 1].grid(True)
axs[1, 1].set(aspect='equal')
plt.savefig("FiniteMode")
plt.show()


plt.plot(x1,x2,"b-")
plt.xlabel("$x_{1}$")
plt.ylabel("$x_{2}$")
plt.title("Modo finito")
plt.grid(True)
plt.axis("equal")
plt.savefig("FiniteMode1")
plt.show()

print(x1T)
plt.plot(x1T,"r--")
plt.show()

