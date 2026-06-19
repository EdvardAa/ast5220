import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("data/perturbations_k0.100000.txt")
data2 = np.loadtxt("data/perturbations_k0.010000.txt")
data3 = np.loadtxt("data/perturbations_k0.001000.txt")

x_1 = data1[:,0]
theta0_1 = data1[:,1]
theta1_1 = data1[:,2]
theta2_1 = data1[:,3]
Phi_1 = data1[:,4]
Psi_1 = data1[:,5]
Pi_1 = data1[:,6]
delta_cdm_1 = data1[:,7]
delta_b_1 = data1[:,8]
v_cdm_1 = data1[:,9]
v_b_1 = data1[:,10]
nu0_1 = data1[:,11]
nu1_1 = data1[:,12]
nu2_1 = data1[:,13]
theta0p_1 = data1[:,14]
theta1p_1 = data1[:,15]
theta2p_1 = data1[:,16]

x_2 = data2[:,0]
theta0_2 = data2[:,1]
theta1_2 = data2[:,2]
theta2_2 = data2[:,3]
Phi_2 = data2[:,4]
Psi_2 = data2[:,5]
Pi_2 = data2[:,6]
delta_cdm_2 = data2[:,7]
delta_b_2 = data2[:,8]
v_cdm_2 = data2[:,9]
v_b_2 = data2[:,10]
nu0_2 = data2[:,11]
nu1_2 = data2[:,12]
nu2_2 = data2[:,13]
theta0p_2 = data2[:,14]
theta1p_2 = data2[:,15]
theta2p_2 = data2[:,16]


x_3 = data3[:,0]
theta0_3 = data3[:,1]
theta1_3 = data3[:,2]
theta2_3 = data3[:,3]
Phi_3 = data3[:,4]
Psi_3 = data3[:,5]
Pi_3 = data3[:,6]
delta_cdm_3 = data3[:,7]
delta_b_3 = data3[:,8]
v_cdm_3 = data3[:,9]
v_b_3 = data3[:,10]
nu0_3 = data3[:,11]
nu1_3 = data3[:,12]
nu2_3 = data3[:,13]
theta0p_3 = data3[:,14]
theta1p_3 = data3[:,15]
theta2p_3 = data3[:,16]

plt.plot(x_1, Phi_1, label=r"$k = 0.1$")
plt.plot(x_1, Psi_1, label=r"$k = 0.1$", linestyle='--')
plt.plot(x_1, Phi_1 + Psi_1, label=r"$k = 0.1$", linestyle=':')
plt.plot(x_2, Phi_2, label=r"$k = 0.01$")
plt.plot(x_2, Psi_2, label=r"$k = 0.01$", linestyle='--')
plt.plot(x_2, Phi_2 + Psi_2, label=r"$k = 0.01$", linestyle=':')
plt.plot(x_3, Phi_3, label=r"$k = 0.001$")
plt.plot(x_3, Psi_3, label=r"$k = 0.001$", linestyle='--')
plt.plot(x_3, Phi_3 + Psi_3, label=r"$k = 0.001$", linestyle=':')
plt.xlabel("x", fontsize=18)
plt.ylabel("Perturbations", fontsize=18)
plt.title("Metric perturbations", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Gravitational_Perturbations.pdf")
plt.show()

plt.plot(x_1, theta2_1, label=r"$\Theta_2$")
plt.plot(x_1, nu2_1, label=r"$\nu_2$", linestyle='--')
plt.xlabel("x", fontsize=18)
plt.ylabel("Quadropoles", fontsize=18)
plt.title("Photon and neutrino quadropoles, $k = 0.1$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Quadropoles_k0.1.pdf")
plt.show()

plt.plot(x_2, theta2_2, label=r"$\Theta_2$")
plt.plot(x_2, nu2_2, label=r"$\nu_2$", linestyle='--')
plt.xlabel("x", fontsize=18)
plt.ylabel("Quadropoles", fontsize=18)
plt.title("Photon and neutrino quadropoles, $k = 0.01$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Quadropoles_k0.01.pdf")
plt.show()

plt.plot(x_3, theta2_3, label=r"$\Theta_2$")
plt.plot(x_3, nu2_3, label=r"$\nu_2$", linestyle='--')
plt.xlabel("x", fontsize=18)
plt.ylabel("Quadropoles", fontsize=18)
plt.title("Photon and neutrino quadropoles, $k = 0.001$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Quadropoles_k0.001.pdf")
plt.show()

plt.plot(x_1, delta_cdm_1, label=r"$\delta_{cdm}$ for $k = 0.1$")
plt.plot(x_1, delta_b_1, label=r"$\delta_b$ for $k = 0.1$", linestyle='--')
plt.plot(x_2, delta_cdm_2, label=r"$\delta_{cdm}$ for $k = 0.01$")
plt.plot(x_2, delta_b_2, label=r"$\delta_b$ for $k = 0.01$", linestyle='--')
plt.plot(x_3, delta_cdm_3, label=r"$\delta_{cdm}$ for $k = 0.001$")
plt.plot(x_3, delta_b_3, label=r"$\delta_b$ for $k = 0.001$", linestyle='--')
plt.yscale("log")
plt.xlabel("x", fontsize=18)
plt.ylabel("Density Perturbations", fontsize=18)
plt.title("$\delta$ for baryons and dark matter", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Density_Perturbations_b_cdm.pdf")
plt.show()

plt.plot(x_1, v_cdm_1, label=r"$v_{cdm}$ for $k = 0.1$")
plt.plot(x_1, v_b_1, label=r"$v_b$ for $k = 0.1$", linestyle='--')
plt.plot(x_2, v_cdm_2, label=r"$v_{cdm}$ for $k = 0.01$")
plt.plot(x_2, v_b_2, label=r"$v_b$ for $k = 0.01$", linestyle='--')
plt.plot(x_3, v_cdm_3, label=r"$v_{cdm}$ for $k = 0.001$")
plt.plot(x_3, v_b_3, label=r"$v_b$ for $k = 0.001$", linestyle='--')
plt.yscale("log")
plt.xlabel("x", fontsize=18)
plt.ylabel("Velocity Perturbations", fontsize=18)
plt.title("$v$ for baryons and dark matter", fontsize=20)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("figures3/Velocity_Perturbations_b_cdm.pdf")
plt.show()

plt.plot(x_1, 4*theta0_1, label=r"$\delta_\gamma$, $k = 0.1$")
plt.plot(x_1, 4*nu0_1, label=r"$\delta_{\nu}$, $k = 0.1$", linestyle='--')
plt.plot(x_2, 4*theta0_2, label=r"$\delta_\gamma$, $k = 0.01$")
plt.plot(x_2, 4*nu0_2, label=r"$\delta_{\nu}$, $k = 0.01$", linestyle='--')
plt.plot(x_3, 4*theta0_3, label=r"$\delta_\gamma$, $k = 0.001$")
plt.plot(x_3, 4*nu0_3, label=r"$\delta_{\nu}$, $k = 0.001$", linestyle='--')
plt.title("$\delta$ for photons and neutrinos", fontsize=20)
plt.xlabel("x", fontsize=18)
plt.ylabel("Density Perturbations", fontsize=18)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Density_Perturbations_gamma_nu.pdf")
plt.show()

plt.plot(x_1, -3*theta1_1, label=r"$v_\gamma$")
plt.plot(x_1, -3*nu1_1, label=r"$v_{\nu}$", linestyle='--')
plt.title("$v$ for photons and neutrinos, $k = 0.1$", fontsize=20)
plt.xlabel("x", fontsize=18)
plt.ylabel("Velocity Perturbations", fontsize=18)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Velocity_Perturbations_gamma_nu_k0.1.pdf")
plt.show()

plt.plot(x_2, -3*theta1_2, label=r"$v_\gamma$")
plt.plot(x_2, -3*nu1_2, label=r"$v_{\nu}$", linestyle='--')
plt.title("$v$ for photons and neutrinos, $k = 0.01$", fontsize=20)
plt.xlabel("x", fontsize=18)
plt.ylabel("Velocity Perturbations", fontsize=18)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Velocity_Perturbations_gamma_nu_k0.01.pdf")
plt.show()

plt.plot(x_3, -3*theta1_3, label=r"$v_\gamma$")
plt.plot(x_3, -3*nu1_3, label=r"$v_{\nu}$", linestyle='--')
plt.title("$v$ for photons and neutrinos, $k = 0.001$", fontsize=20)
plt.xlabel("x", fontsize=18)
plt.ylabel("Velocity Perturbations", fontsize=18)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Velocity_Perturbations_gamma_nu_k0.001.pdf")
plt.show()

plt.plot(x_1, theta0p_1, label=r"$\Theta_0^P$")
plt.plot(x_1, theta1p_1, label=r"$\Theta_1^P$")
plt.plot(x_1, theta2p_1, label=r"$\Theta_2^P$")
plt.xlabel("x", fontsize=18)
plt.ylabel("Multipoles", fontsize=18)
plt.title("Polarization multipoles, $k = 0.1$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Polarization_Multipoles_k0.1.pdf")
plt.show()

plt.plot(x_2, theta0p_2, label=r"$\Theta_0^P$")
plt.plot(x_2, theta1p_2, label=r"$\Theta_1^P$")
plt.plot(x_2, theta2p_2, label=r"$\Theta_2^P$")
plt.xlabel("x", fontsize=18)
plt.ylabel("Multipoles", fontsize=18)
plt.title("Polarization multipoles, $k = 0.01$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Polarization_Multipoles_k0.01.pdf")
plt.show()

plt.plot(x_3, theta0p_3, label=r"$\Theta_0^P$")
plt.plot(x_3, theta1p_3, label=r"$\Theta_1^P$")
plt.plot(x_3, theta2p_3, label=r"$\Theta_2^P$")
plt.xlabel("x", fontsize=18)
plt.ylabel("Multipoles", fontsize=18)
plt.title("Polarization multipoles, $k = 0.001$", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures3/Polarization_Multipoles_k0.001.pdf")
plt.show()