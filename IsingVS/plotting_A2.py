import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 
import scipy.special as scp

def mean_abs_mag_analytical(beta):
    return (1-(np.sinh(2*beta))**(-4))**(1/8)


def mean_mag_analytical(beta,h):
    return (np.sinh(beta * h)/np.sqrt(np.sinh(beta * h)**2 + np.exp(-4*beta)))

def mean_energy_density(beta):
    chi = 2*(np.tanh(2*beta)/np.cosh(2*beta))
    return   -1/np.tanh(2*beta) * (1 + (2* np.tanh(2*beta)**2 - 1)*(2/np.pi)*scp.ellipk(chi))

betas_1 = np.linspace(0.440687, 0.5, 1000)
betas_2 = np.linspace(0.5, 1, 100)
betas_3 = np.linspace(0, 1, 1000)


 
fig, ax = plt.subplots()
#ax.plot(betas_1, mean_abs_mag_analytical(betas_1))
#ax.plot(betas_2, mean_abs_mag_analytical(betas_2))
#ax.plot(betas_3, mean_mag_analytical(betas_3))
ax.plot(betas_3, mean_energy_density(betas_3), color = "r", label = "analytical solution")

ax.set_ylabel(r'mean energy density $\langle e \rangle$')
ax.set_xlabel(r'temperature $\beta$')
plt.title("explicit mean energy density of small planar lattices")
plt.grid()
#ax.scatter(betas_meas, mea_mag_ana)

LS = [2,3,4]
colors = ["b", "g", "orange"]
for i,L in enumerate(LS):
    df = pd.read_table("output_files/explicit_vals_L=" + str(L) + ".txt")
    betas_meas = np.array(df.iloc[:,0])
    abs_mag_meas = np.array(df.iloc[:,3])
    e_meas = np.array(df.iloc[:,1])
    ax.scatter(betas_meas, e_meas, marker = "s", facecolors = "none", color = colors[i], label = r'$L =$' + str(L))
plt.legend()
plt.show()


fig, ax = plt.subplots()
ax.plot(betas_1, mean_abs_mag_analytical(betas_1),color = "r", label = "analytical solution")
ax.plot(betas_2, mean_abs_mag_analytical(betas_2), color = "r")
ax.axhline(0, 0, 0.4408, color = "r")
ax.axvline(0.4406, 0, 0.4, color = "r")
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_ylabel(r'absolut mean magnetisation per particle $\langle |m| \rangle$')
ax.set_xlabel(r'temperature $\beta$')
plt.title("absolut mean magnetisation per lattice point of small planar lattices")
plt.grid()
#ax.scatter(betas_meas, mea_mag_ana)

LS = [2,3,4]
colors = ["b", "g", "orange"]
for i,L in enumerate(LS):
    df = pd.read_table("output_files/explicit_vals_L=" + str(L) + ".txt")
    betas_meas = np.array(df.iloc[:,0])
    abs_mag_meas = np.array(df.iloc[:,3])
    e_meas = np.array(df.iloc[:,1])
    ax.scatter(betas_meas, abs_mag_meas, marker = "s", facecolors = "none", color = colors[i], label = r'$L =$' + str(L))
plt.legend()
plt.show()
