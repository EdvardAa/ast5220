import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, mu, sigma):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5 * ((x - mu)/sigma)**2)

data = np.loadtxt("cosmology.txt")
x = data[:,0]
eta = data[:,1]
Hp = data[:,2]
dHpdx = data[:,3]
OmegaB = data[:,4]
OmegaCDM = data[:,5]
OmegaLambda = data[:,6]
OmegaR = data[:,7]
OmegaNu = data[:,8]
OmegaK = data[:,9]
t = data[:,10]
ddHpddx = data[:,11]

omega_matter = 0
omega_radiation = 1/3
omega_dark_energy = -1
omegas = [omega_matter, omega_radiation]
omega_strings = ["Matter domination", "Radiation domination"]
colors = ["orange", "red"]

a = np.exp(x)
a0 = 1.0
H = Hp * a
z = a0/a - 1

Omega_relativistic = OmegaR + OmegaNu
Omega_matter = OmegaB + OmegaCDM
c = 2.99792458e8
Gyr = 1e9*365*24*3600
Mpc = 3.08567758e22




plt.plot(x, t/Gyr, label="t(x) cosmic time")
plt.plot(x, eta/(c*Gyr), label=r"$\eta(x)$ conformal time")
plt.xlabel("x", fontsize=18)
plt.ylabel("t (Gyr)", fontsize=18)
plt.title("Time as a function of x", fontsize=20)
plt.yscale("log")
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures/time_vs_x.pdf")
plt.show()

plt.plot(x, dHpdx/Hp, label="dHp/dx / Hp")
for i in range(2):
    plt.axhline(-(1+3*omegas[i])/2, color=colors[i], linestyle='--', label=omega_strings[i])
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\frac{dH_p/dx}{H_p}$", fontsize=18)
plt.title(r"$\frac{dH_p/dx}{H_p}$ as a function of $x$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures/dHpdx_over_Hp_vs_x.pdf")
plt.show()

plt.plot(x, ddHpddx/Hp, label="ddHp/ddx / Hp")
for i in range(2):
    plt.axhline(((1+3*omegas[i])/2)**2, color=colors[i], linestyle='--', label=omega_strings[i])
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\frac{d^2H_p/dx^2}{dH_p/dx}$", fontsize=18)
plt.title(r"$\frac{d^2H_p/dx^2}{dH_p/dx}$ as a function of $x$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures/ddHpddx_over_Hp_vs_x.pdf")
plt.show()

plt.plot(x, Hp*Mpc/(1000*100))
plt.xlabel("x", fontsize=18)
plt.ylabel("Hp", fontsize=18)
plt.title(r"$H_p(x)$ ($\frac{100km}{s \cdot Mpc}$)", fontsize=18)
plt.yscale("log")
plt.tight_layout()
plt.savefig("figures/Hp_vs_x.pdf")
plt.show()

plt.plot(x, eta*Hp/c)
plt.xlabel("x", fontsize=18)
plt.ylabel(r"$\eta$*Hp", fontsize=18)
plt.title(r"$\frac{\eta H_p}{c}$ as a function of $x$", fontsize=20)
plt.tight_layout()
plt.savefig("figures/etaHp_over_c_vs_x.pdf")
plt.show()

plt.plot(x, Omega_relativistic, label="Radiation")
plt.plot(x, Omega_matter, label="Matter")
plt.plot(x, OmegaLambda, label="Dark energy")
plt.xlabel("x", fontsize=18)
plt.ylabel("Omega", fontsize=18)
plt.title("Density parameters as a function of x", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures/Omega_vs_x.pdf")
plt.show()

plt.plot(t/Gyr, np.exp(x))
plt.xlabel("t (Gyr)", fontsize=18)
plt.ylabel("a(t)", fontsize=18)
plt.title("Scale factor a(t) as a function of time t", fontsize=20)
plt.savefig("figures/a_vs_t.pdf")
plt.show()

data_supernova = np.loadtxt("results_supernovafitting.txt", skiprows=200)
chi2 = data_supernova[:,0]
chi2min = np.min(chi2)
h = data_supernova[:,1]
OmegaM_1sigma = data_supernova[chi2 < chi2min + 3.53,2]
OmegaK_1sigma = data_supernova[chi2 < chi2min + 3.53,3]
OmegaLambda_1sigma = 1 - OmegaM_1sigma - OmegaK_1sigma

OmegaM_2sigma = data_supernova[chi2 < chi2min + 8.02,2]
OmegaK_2sigma = data_supernova[chi2 < chi2min + 8.02,3]
OmegaLambda_2sigma = 1 - OmegaM_2sigma - OmegaK_2sigma

OmegaM = data_supernova[:,2]
OmegaK = data_supernova[:,3]
OmegaLambda = 1 - OmegaM - OmegaK

OmegaM_analytical = np.linspace(0.0, 1.0, 100)
OmegaLambda_analytical = 1.0 - OmegaM_analytical

plt.plot(OmegaM_analytical, OmegaLambda_analytical, label="Flat universe", linestyle="dashed")
plt.scatter(OmegaM_2sigma, OmegaLambda_2sigma, c='r', label=r"2$\sigma$ constraint")
plt.scatter(OmegaM_1sigma, OmegaLambda_1sigma, c='y', label=r"1$\sigma$ constraint")
plt.xlabel("OmegaM", fontsize=18)
plt.ylabel("OmegaLambda", fontsize=18)
plt.legend(fontsize=12)
plt.title("Density parameters as a function of OmegaM", fontsize=20)
plt.savefig("figures/OmegaLambda_vs_OmegaM.pdf")
plt.show()

std_M = np.std(OmegaM_1sigma)
std_K = np.std(OmegaK_1sigma)
std_Lambda = np.std(OmegaLambda_1sigma)

mean_M = np.mean(OmegaM)
mean_K = np.mean(OmegaK)
mean_Lambda = np.mean(OmegaLambda)

min_chi2_arg = np.argmin(chi2)
best_fit_M = data_supernova[min_chi2_arg,2]
best_fit_K = data_supernova[min_chi2_arg,3]
best_fit_Lambda = 1 - best_fit_M - best_fit_K

x = np.linspace(0.0, 1.0, 100)
plt.hist(OmegaLambda, bins=40, density=True, label="Samples")
plt.plot(x, gaussian(x, mean_Lambda, std_Lambda), label="Best fit")
plt.axvline(best_fit_Lambda, color='black', linestyle='--', label="Best fit Lambda")
plt.xlabel("OmegaLambda", fontsize=18)
plt.ylabel("Number of samples", fontsize=18)
plt.title(r"Distribution of $\Omega_\Lambda$ samples", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures/OmegaLambda_histogram.pdf")
plt.show()

data_supernova = np.loadtxt("data/supernovadata.txt", skiprows=1)
z_data = data_supernova[:,0]
dL_data = data_supernova[:,1]
error_data = data_supernova[:,2]
plt.errorbar(z_data, dL_data, yerr=error_data, fmt='o', label="Data with error bars", color='red')
plt.xlabel("Redshift z", fontsize=18)
plt.ylabel("Luminosity distance dL (Mpc)", fontsize=18)
plt.title("Luminosity distance as a function of redshift", fontsize=20)
plt.legend(fontsize=12)
plt.savefig("figures/dL_vs_z.pdf")
plt.show()