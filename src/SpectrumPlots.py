import numpy as np
import matplotlib.pyplot as plt

cell_data = np.loadtxt("data/cells.txt")
ells = cell_data[:, 0]
cell_TT = cell_data[:, 1]
cells_EE = cell_data[:, 2]
cells_TE = cell_data[:, 3]

k_data = np.loadtxt("data/k_spectra.txt")
k = k_data[:,0]
matter_PS = k_data[:,1]
ThetaT15 = k_data[:,2]
ThetaT500 = k_data[:,3]
ThetaT1000 = k_data[:,4]
ThetaE15 = k_data[:,5]
ThetaE500 = k_data[:,6]
ThetaE1000 = k_data[:,7]

low_l_TT_data = np.loadtxt("data/low_l_TT_data.txt")
low_l_TT_ell = low_l_TT_data[:, 0]
low_l_TT_Cell = low_l_TT_data[:, 1]
low_l_TT_err_UP = low_l_TT_data[:, 2]
low_l_TT_err_DOWN = low_l_TT_data[:, 3]

high_l_TT_data = np.loadtxt("data/high_l_TT_data.txt")
high_l_TT_ell = high_l_TT_data[:, 0]
high_l_TT_Cell = high_l_TT_data[:, 1]
high_l_TT_err_UP = high_l_TT_data[:, 2]
high_l_TT_err_DOWN = high_l_TT_data[:, 3]

high_l_EE_data = np.loadtxt("data/high_l_EE_data.txt")
high_l_EE_ell = high_l_EE_data[:, 0]
high_l_EE_Cell = high_l_EE_data[:, 1]
high_l_EE_err_UP = high_l_EE_data[:, 2]
high_l_EE_err_DOWN = high_l_EE_data[:, 3]

high_l_TE_data = np.loadtxt("data/high_l_TE_data.txt")
high_l_TE_ell = high_l_TE_data[:, 0]
high_l_TE_Cell = high_l_TE_data[:, 1]
high_l_TE_err_UP = high_l_TE_data[:, 2]
high_l_TE_err_DOWN = high_l_TE_data[:, 3]

low_l_errors = [low_l_TT_err_DOWN, low_l_TT_err_UP]
high_l_errors = [high_l_TT_err_DOWN, high_l_TT_err_UP]
ee_errors = [high_l_EE_err_DOWN, high_l_EE_err_UP]
te_errors = [high_l_TE_err_DOWN, high_l_TE_err_UP]

plt.figure(figsize=(10, 6))

plt.plot(ells, cell_TT, label="Theory (TT)", color="black", linewidth=2, zorder=1)

# Plot Low-l observations with asymmetric errors
plt.errorbar(
    low_l_TT_ell, 
    low_l_TT_Cell, 
    yerr=low_l_errors, 
    fmt='o', 
    markersize=4,
    color='crimson', 
    ecolor='crimson', 
    elinewidth=1.2, 
    capsize=2, 
    label="Low-$\ell$ Data",
    zorder=2
)

# Plot High-l observations with asymmetric errors
plt.errorbar(
    high_l_TT_ell, 
    high_l_TT_Cell, 
    yerr=high_l_errors, 
    fmt='s', 
    markersize=3,
    color='royalblue', 
    ecolor='royalblue', 
    elinewidth=1.0, 
    capsize=1.5, 
    label="High-$\ell$ Data",
    zorder=3
)

plt.xlabel(r"Multipole moment ($\ell$)", fontsize=18)
plt.ylabel(r"$D_\ell \equiv \ell(\ell+1) C_\ell / 2\pi$ [$\mu\text{K}^2$]", fontsize=18)
plt.title("CMB Temperature Anisotropy Power Spectrum ($TT$)", fontsize=20)

plt.xscale("log")

plt.xlim(2, max(ells))
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(loc="upper right", fontsize=11)

# Save and render configuration
plt.tight_layout()
plt.savefig("figures4/TT_spectrum.pdf", dpi=300)
plt.show()

plt.figure(figsize=(10, 6))

plt.plot(ells, cells_EE, label="Theory (EE)", color="black", linewidth=2, zorder=1)

plt.errorbar(
    high_l_EE_ell, 
    high_l_EE_Cell, 
    yerr=ee_errors, 
    fmt='o', 
    markersize=3.5,
    color='darkorchid', 
    ecolor='darkorchid', 
    elinewidth=1.0, 
    capsize=1.5, 
    label="High-$\ell$ EE Data",
    zorder=2
)

plt.xscale("log")

plt.xlabel(r"Multipole moment ($\ell$)", fontsize=18)
plt.ylabel(r"$D_\ell^{EE} \equiv \ell(\ell+1) C_\ell^{EE} / 2\pi$ [$\mu\text{K}^2$]", fontsize=18)
plt.title("CMB E-Mode Polarization Power Spectrum ($EE$)", fontsize=20)

plt.xlim(2, max(ells))
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(loc="upper right", fontsize=11)
plt.tight_layout()
plt.savefig("figures4/EE_spectrum.pdf", dpi=300)
plt.show()

plt.figure(figsize=(10, 6))

plt.plot(ells, cells_TE, label="Theory (TE)", color="black", linewidth=2, zorder=1)

plt.errorbar(
    high_l_TE_ell, 
    high_l_TE_Cell, 
    yerr=te_errors, 
    fmt='o', 
    markersize=3.5,
    color='forestgreen', 
    ecolor='forestgreen', 
    elinewidth=1.0, 
    capsize=1.5, 
    label="High-$\ell$ TE Data",
    zorder=2
)

plt.xscale("log")

plt.xlabel(r"Multipole moment ($\ell$)", fontsize=18)
plt.ylabel(r"$D_\ell^{TE} \equiv \ell(\ell+1) C_\ell^{TE} / 2\pi$ [$\mu\text{K}^2$]", fontsize=18)
plt.title("CMB Temperature-Polarization Cross Spectrum ($TE$)", fontsize=20)

plt.xlim(2, max(ells))
plt.axhline(0, color='gray', linestyle='-', alpha=0.5, zorder=0) # baseline reference
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(loc="upper right", fontsize=11)
plt.tight_layout()
plt.savefig("figures4/TE_spectrum.pdf", dpi=300)
plt.show()

dr7_data = np.loadtxt("data/dr7_data.txt")
k_dr7 = dr7_data[:,0]
P_dr7 = dr7_data[:,1]
P_error_dry = dr7_data[:,2]

wmap_act_data = np.loadtxt("data/dr7_data.txt")
k_wmap_act = wmap_act_data[:,0]
P_wmap_act = wmap_act_data[:,1]
P_error_wmap_act = np.abs(wmap_act_data[:,2] - P_wmap_act)

dr7_data = np.loadtxt("data/dr7_data.txt")
k_dr7, P_dr7, P_error_dr7 = dr7_data[:, 0], dr7_data[:, 1], dr7_data[:, 2]

# Ensure wmap data is correctly loaded (if it's a different file, adjust path)
wmap_act_data = np.loadtxt("data/wmap_act_data.txt") 
k_wmap, P_wmap, P_error_wmap = wmap_act_data[:, 0], wmap_act_data[:, 1], np.abs(wmap_act_data[:, 2] - wmap_act_data[:, 1])

plt.figure(figsize=(10, 6))

# 1. Plot the theoretical solution
plt.plot(k, matter_PS, label="Theory", color="black", linewidth=2, zorder=1)

# 2. Plot DR7 observations
plt.errorbar(
    k_dr7, P_dr7, 
    yerr=P_error_dr7, 
    fmt='o', markersize=4,
    color='crimson', ecolor='crimson',
    elinewidth=1.2, capsize=2,
    label="SDSS Galaxies (DR7 LRG)", zorder=2
)

# 3. Plot WMAP/ACT observations
plt.errorbar(
    k_wmap, P_wmap, 
    yerr=P_error_wmap, 
    fmt='s', markersize=3,
    color='royalblue', ecolor='royalblue',
    elinewidth=1.0, capsize=1.5,
    label="Cosmic Microwave Background (WMAP+ACT)", zorder=3
)

# Formatting
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Wavenumber $k$ [h/Mpc]", fontsize=18)
plt.ylabel(r"Matter Power Spectrum $P(k)$ [(Mpc/h)$^3$]", fontsize=18)
plt.title("Linear Matter Power Spectrum", fontsize=20)

plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(loc="lower left", fontsize=11)

plt.tight_layout()
plt.savefig("figures4/Matter_PS.pdf", dpi=300)
plt.show()

plt.plot(k, ThetaT15, label=r"$\ell = 15$")
plt.plot(k, ThetaT500, label=r"$\ell = 500$")
plt.plot(k, ThetaT1000, label=r"$\ell = 1000$")
plt.xscale("log")
plt.xlabel("Wavenumber k [h/Mpc]", fontsize=18)
plt.ylabel("Multipoles $\Theta_\ell$", fontsize=18)
plt.title("Photon temperature transfer functions", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures4/Photon_Multipoles_x0.pdf")
plt.show()

plt.plot(k, ThetaE15, label=r"$\ell = 15$")
plt.plot(k, ThetaE500, label=r"$\ell = 500$")
plt.plot(k, ThetaE1000, label=r"$\ell = 1000$")
plt.xscale("log")
plt.xlabel("Wavenumber k [h/Mpc]", fontsize=18)
plt.ylabel("Multipoles $\Theta^E_\ell$", fontsize=18)
plt.title("Photon polarization transfer functions", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures4/Polarization_Multipoles_x0.pdf")
plt.show()

fig, ax=plt.subplots(3)
ax[0].plot(k, ThetaT15**2 / k, label="$\ell = 15$")
ax[1].plot(k, ThetaT500**2 / k, label="$\ell = 500$")
ax[2].plot(k, ThetaT1000**2 / k, label="$\ell = 1000$")
ax[0].set_xscale("log")
ax[1].set_xscale("log")
ax[2].set_xscale("log")
ax[0].set_title("Photon temperature $C_\ell$ integrand", fontsize=20)
ax[1].set_ylabel("Integrand $|\Theta_\ell|^2 / k$", fontsize=18)
ax[2].set_xlabel("Wavenumber k [h/Mpc]", fontsize=18)
ax[0].legend(fontsize=12)
ax[1].legend(fontsize=12)
ax[2].legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures4/Photon_integrand_x0.pdf")
plt.show()

plt.plot(k, ThetaE15**2 / k, label=r"$\ell = 15$")
plt.plot(k, ThetaE500**2 / k, label=r"$\ell = 500$")
plt.plot(k, ThetaE1000**2 / k, label=r"$\ell = 1000$")
plt.xscale("log")
plt.xlabel("Wavenumber k [h/Mpc]", fontsize=18)
plt.ylabel("Integrand $|\Theta^E_\ell|^2 / k$", fontsize=18)
plt.title("Photon polarization $C_\ell$ integrand", fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("figures4/Polarization_integrand_x0.pdf")
plt.show()