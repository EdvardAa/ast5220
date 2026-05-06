import numpy as np
import matplotlib.pyplot as plt

x_decoupling_recombination = -6.9853
x_decoupling_recombination_Saha = -7.36938
x_decoupling_scattering = -6.98843
x_decoupling_Saha_scattering = -5.18139

data = np.loadtxt("recombination.txt")
x = data[:,0]
index = np.argmin(np.abs(x - (-12.0)))
x = x[index:]
Xe = data[index:,1]
ne = data[index:,2]
tau = data[index:,3]
dtau_dx = data[index:,4]
ddtau_ddx = data[index:,5]
g_tilde = data[index:,6]
dg_tildedx = data[index:,7]
ddg_tildeddx = data[index:,8]

Xe_reion = data[index:,9]
ne_reion = data[index:,10]
tau_reion = data[index:,11]
dtau_dx_reion = data[index:,12]
ddtau_ddx_reion = data[index:,13]
g_tilde_reion = data[index:,14]
dg_tildedx_reion = data[index:,15]
ddg_tildeddx_reion = data[index:,16]
Xe_saha = data[index:,17]
s = data[index:,18]

plt.plot(x, ne, label=r"$n_e(x)$")
plt.plot(x, ne_reion, label=r"$n_e$ with reionization", linestyle="dashed")
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$n_e$ [$1 / m^3$]", fontsize=18)
plt.title(r"$n_e$ as a function of $x$", fontsize=20)
plt.legend(fontsize=12)
plt.yscale("log")
plt.savefig("figures2/ne_vs_x.pdf")
plt.show()

plt.plot(x, Xe, label=r"$X_e(x)$")
plt.plot(x, Xe_saha, label=r"$X_e$ from Saha equation")
plt.axvline(x_decoupling_recombination, color='blue', linestyle=':', 
            label=f'Rec (Peebles) x={x_decoupling_recombination:.2f}')
plt.axvline(x_decoupling_recombination_Saha, color='red', linestyle=':', 
            label=f'Rec (Saha) x={x_decoupling_recombination_Saha:.2f}')
plt.yscale("log")
plt.ylim(1e-4, 2)
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$X_e$ [1]", fontsize=18)
plt.title(r"X_e with Saha and Peebles vs only Saha", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/Xe_Saha_vs_Peebles.pdf")
plt.show()

plt.plot(x, Xe, label=r"$X_e(x)$")
plt.plot(x, Xe_reion, label=r"$X_e$ with reionization", linestyle="dashed")
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$X_e$ [1]", fontsize=18)
plt.title(r"Electron fraction $X_e$ as a function of $x$", fontsize=20)
plt.legend(fontsize=12)
plt.yscale("log")
plt.savefig("figures2/Xe_vs_x_reion.pdf")
plt.show()

plt.plot(x[:-3],tau[:-3], label=r"$\tau(x)$")
plt.plot(x[:-3], -dtau_dx[:-3], label=r"$-d\tau/dx$")
plt.plot(x[:-3], ddtau_ddx[:-3], label=r"$d^2\tau/dx^2$")
plt.axvline(x_decoupling_scattering, color='green', linestyle=':', 
            label=f'Scattering x={x_decoupling_scattering:.2f}')
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\tau$ [1]", fontsize=18)
plt.yscale("log")
plt.title(r"Optical depth $\tau$ and derivatives", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/tau_and_derivatives_vs_x.pdf")
plt.show()

plt.plot(x[:-3], tau_reion[:-3], label=r"$\tau$ with reionization")
plt.plot(x[:-3], -dtau_dx_reion[:-3], label=r"$-d\tau/dx$ with reionization")
plt.plot(x[:-3], ddtau_ddx_reion[:-3], label=r"$d^2\tau/dx^2$ with reionization")
plt.axvline(x_decoupling_scattering, color='green', linestyle=':', 
            label=f'Scattering x={x_decoupling_scattering:.2f}')
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\tau$ [1]", fontsize=18)
plt.yscale("log")
plt.title(r"Optical depth $\tau$ with reionization", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/tau_and_derivatives_vs_x_reion.pdf")
plt.show()

plt.plot(x, g_tilde / np.sum(np.abs(g_tilde)), label=r"$\tilde{g}(x)$")
plt.plot(x, dg_tildedx / np.sum(np.abs(dg_tildedx)), label=r"$d\tilde{g}/dx$", linestyle="--")
plt.plot(x, ddg_tildeddx / np.sum(np.abs(ddg_tildeddx)), label=r"$d^2\tilde{g}/dx^2$", linestyle=":")
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\tilde{g}$ [1]", fontsize=18)
plt.xlim(-8,0)
plt.title(r"Scaled visibility function $\tilde{g}$ and derivatives", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/g_tilde_and_derivatives_vs_x.pdf")
plt.show()

plt.plot(x, g_tilde_reion / np.sum(np.abs(g_tilde_reion)), label=r"$\tilde{g}$ with reionization")
plt.plot(x, dg_tildedx_reion / np.sum(np.abs(dg_tildedx_reion)), label=r"$d\tilde{g}/dx$ with reionization", linestyle="--")
plt.plot(x, ddg_tildeddx_reion / np.sum(np.abs(ddg_tildeddx_reion)), label=r"$d^2\tilde{g}/dx^2$ with reionization", linestyle=":")
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\tilde{g}$ [1]", fontsize=18)
plt.xlim(-8,0)
plt.title(r"Scaled visibility function $\tilde{g}$ with reionization", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/g_tilde_and_derivatives_vs_x_reion.pdf")
plt.show()

plt.plot(x, s, label=r"$s(x)$")
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$s$ [Mpc]", fontsize=18)
plt.title(r"Sound horizon $s$ as a function of $x$", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures2/s_vs_x.pdf")
plt.show()