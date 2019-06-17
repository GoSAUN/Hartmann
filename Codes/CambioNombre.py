import numpy as np 

k = input()

x = np.loadtxt("x.dat")
y = np.loadtxt("y.dat")
x = np.loadtxt("x.dat")

u = np.loadtxt("vx.dat")
v = np.loadtxt("vy.dat")
w = np.loadtxt("vz.dat")

Bx = np.loadtxt("Bx.dat")
By = np.loadtxt("By.dat")
Bz = np.loadtxt("Bz.dat")

np.savetxt("x"+str(k)+".dat", x)
np.savetxt("y"+str(k)+".dat", y)
np.savetxt("z"+str(k)+".dat", u)
np.savetxt("vx"+str(k)+".dat", v)
np.savetxt("vy"+str(k)+".dat", x)
np.savetxt("vz"+str(k)+".dat", y)
np.savetxt("Bx"+str(k)+".dat", u)
np.savetxt("By"+str(k)+".dat", v)
np.savetxt("Bz"+str(k)+".dat", v)
