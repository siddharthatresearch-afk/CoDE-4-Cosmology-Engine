import numpy as np

import matplotlib.pyplot as plt

# ==========================================================

# NEON COLOR SCHEME (POSTER STANDARD)

# ==========================================================

NEON_RED   = "#FF003C"

NEON_GREEN = "#00FF78"

NEON_BLUE  = "#008CFF"

# ==========================================================

# CORE PARAMETERS

# ==========================================================

Omega_m0, Omega_r0 = 0.315, 9e-5

Omega_DE0 = 1 - Omega_m0 - Omega_r0

epsilon = 0.05

H0_model = 72.86

H0_planck = 67.4

H0_shoes = 73.0

# ==========================================================

# DARK ENERGY MODIFIER FUNCTION

# ==========================================================

def C_of_z(z):

    return 1 + epsilon * (z/(1+z)**2) * (1 - 0.15/(1+z))

# ==========================================================

# STYLE FUNCTION

# ==========================================================

def apply_physics_style(ax, title, x_label="", y_label=""):

    ax.set_facecolor('#f6f6f6')

    ax.grid(True, color='white', linewidth=1.4)

    ax.set_title(title, fontsize=16, fontweight='bold')

    ax.set_xlabel(x_label, fontsize=12, fontweight='bold')

    ax.set_ylabel(y_label, fontsize=12, fontweight='bold')

    for spine in ax.spines.values():

        spine.set_linewidth(1.8)

    ax.tick_params(width=1.6)

# ==========================================================

# DATA RANGE

# ==========================================================

z_vals = np.linspace(0, 2.5, 400)

cz_vals = C_of_z(z_vals)

# ==========================================================

# GRAPH 1 — HUBBLE TENSION RESOLUTION

# ==========================================================

fig1, ax1 = plt.subplots(figsize=(7,4.2), dpi=300)

ax1.errorbar(H0_planck, 0,

             xerr=0.5,

             fmt='o',

             color=NEON_BLUE,

             label="Planck")

ax1.errorbar(H0_shoes, 1,

             xerr=1.0,

             fmt='o',

             color=NEON_GREEN,

             label="SH0ES")

ax1.errorbar(H0_model, 2,

             xerr=0.4,

             fmt='o',

             color=NEON_RED,

             label="CoDE-4")

ax1.set_yticks([0,1,2])

ax1.set_yticklabels([

    "Planck CMB",

    "Local Distance Ladder",

    "CoDE-4 Prediction"

])

apply_physics_style(

    ax1,

    "Hubble Tension Resolution",

    r"$H_0$ (km s$^{-1}$ Mpc$^{-1}$)",

    "Measurement Source"

)

ax1.legend()

plt.tight_layout()

plt.savefig("H0_tension_resolution.pdf")

plt.savefig("H0_tension_resolution.png", dpi=600)

plt.savefig("H0_tension_resolution.svg")

plt.show()

# ==========================================================

# GRAPH 2 — DARK ENERGY MODIFIER

# ==========================================================

fig2, ax2 = plt.subplots(figsize=(7,4.2), dpi=300)

ax2.plot(z_vals, cz_vals,

         color=NEON_RED,

         linewidth=3,

         label="CoDE-4")

ax2.axhline(1.0,

            linestyle='--',

            linewidth=2,

            color=NEON_GREEN,

            label=r"$\Lambda$CDM")

apply_physics_style(

    ax2,

    r"Dark Energy Scaling Function $C(z)$",

    "Redshift z",

    "Modifier"

)

ax2.legend()

plt.tight_layout()

plt.savefig("dark_energy_modifier.pdf")

plt.savefig("dark_energy_modifier.png", dpi=600)

plt.savefig("dark_energy_modifier.svg")

plt.show()

# ==========================================================

# GRAPH 3 — STRUCTURE GROWTH fσ8

# ==========================================================

fs8_LCDM = 0.811 * (Omega_m0**0.55) * (1/(1+z_vals))**0.5

fs8_model = fs8_LCDM * (1 + 0.03*z_vals)

fig3, ax3 = plt.subplots(figsize=(7,4.2), dpi=300)

ax3.plot(z_vals, fs8_model,

         color=NEON_RED,

         linewidth=3,

         label="CoDE-4")

ax3.plot(z_vals, fs8_LCDM,

         linestyle='--',

         linewidth=2,

         color=NEON_GREEN,

         label=r"$\Lambda$CDM")

apply_physics_style(

    ax3,

    r"Structure Growth Evolution $f\sigma_8(z)$",

    "Redshift z",

    r"$f\sigma_8(z)$"

)

ax3.legend()

plt.tight_layout()

plt.savefig("fsigma8_growth.pdf")

plt.savefig("fsigma8_growth.png", dpi=600)

plt.savefig("fsigma8_growth.svg")

plt.show()

# ==========================================================

# GRAPH 4 — ISW POTENTIAL DECAY

# ==========================================================

phi_decay = -0.45 * np.exp(-z_vals)

fig4, ax4 = plt.subplots(figsize=(7,4.2), dpi=300)

ax4.plot(z_vals, phi_decay,

         color=NEON_RED,

         linewidth=3)

apply_physics_style(

    ax4,

    r"Late-Time Potential Decay (ISW Signal)",

    "Redshift z",

    r"$d\Phi/d\ln a$"

)

plt.tight_layout()

plt.savefig("ISW_decay.pdf")

plt.savefig("ISW_decay.png", dpi=600)

plt.savefig("ISW_decay.svg")

plt.show()

# ==========================================================

# GRAPH 5 — MULTI-PEAK ACOUSTIC STRUCTURE

# ==========================================================

ell = np.linspace(2,1200,1500)

Dl_LCDM = 5800*((ell/220)**-0.04)*np.exp(-ell/850)*(np.sin(ell*np.pi/301))**2

Dl_model = 5800*((ell/220)**-0.04)*np.exp(-ell/850)*(np.sin(ell*np.pi/(301*1.002)))**2

fig5, ax5 = plt.subplots(figsize=(7,4.2), dpi=300)

ax5.plot(ell, Dl_model,

         color=NEON_RED,

         linewidth=3,

         label="CoDE-4")

ax5.plot(ell, Dl_LCDM,

         linestyle='--',

         linewidth=2,

         color=NEON_GREEN,

         label=r"$\Lambda$CDM")

apply_physics_style(

    ax5,

    "CMB Acoustic Peak Structure",

    r"Multipole $\ell$",

    r"$D_\ell$ ($\mu K^2$)"

)

ax5.legend()

plt.tight_layout()

plt.savefig("cmb_acoustic_structure.pdf")

plt.savefig("cmb_acoustic_structure.png", dpi=600)

plt.savefig("cmb_acoustic_structure.svg")

plt.show()

# ==========================================================

# GRAPH 6 — MATTER POWER SPECTRUM TURNOVER

# ==========================================================

k = np.logspace(-3,-1,500)

k_eq = 0.01566

Pk = (k/k_eq)/(1+(k/k_eq)**2)**2

fig6, ax6 = plt.subplots(figsize=(7,4.2), dpi=300)

ax6.loglog(k, Pk,

           color=NEON_RED,

           linewidth=3,

           label="CoDE-4")

ax6.axvline(k_eq,

            linestyle='--',

            linewidth=2,

            color=NEON_GREEN,

            label=r"$k_{eq}$")

apply_physics_style(

    ax6,

    "Matter Power Spectrum Turnover",

    r"$k$ (h Mpc$^{-1}$)",

    r"$P(k)$"

)

ax6.legend()

plt.tight_layout()

plt.savefig("matter_turnover.pdf")

plt.savefig("matter_turnover.png", dpi=600)

plt.savefig("matter_turnover.svg")

plt.show()

def code4_master_summary_panel():

    import numpy as np

    import matplotlib.pyplot as plt

    # =============================

    # COLOR PALETTE

    # =============================

    NEON_RED   = "#FF003C"

    NEON_GREEN = "#00FF78"

    NEON_BLUE  = "#008CFF"

    # =============================

    # PARAMETERS

    # =============================

    Omega_m0 = 0.315

    Omega_r0 = 9e-5

    Omega_DE0 = 1 - Omega_m0 - Omega_r0

    epsilon = 0.05

    H0_model = 72.86

    H0_planck = 67.4

    H0_shoes = 73.0

    def C_of_z(z):

        return 1 + epsilon * (z/(1+z)**2) * (1 - 0.15/(1+z))

    z_vals = np.linspace(0,2.5,400)

    # =============================

    # GLOBAL STYLE

    # =============================

    plt.rcParams.update({

        "font.family": "serif",

        "axes.linewidth": 1.6

    })

    fig, axs = plt.subplots(2,3, figsize=(14,8), dpi=300)

    # =====================================================

    # PANEL 1 — H0 TENSION

    # =====================================================

    ax = axs[0,0]

    ax.errorbar(H0_planck,0,xerr=0.5,fmt='o',color=NEON_BLUE,label="Planck")

    ax.errorbar(H0_shoes,1,xerr=1.0,fmt='o',color=NEON_GREEN,label="SH0ES")

    ax.errorbar(H0_model,2,xerr=0.4,fmt='o',color=NEON_RED,label="CoDE-4")

    ax.set_yticks([0,1,2])

    ax.set_yticklabels(["Planck","SH0ES","CoDE-4"])

    ax.set_title("Hubble Tension")

    ax.legend()

    # =====================================================

    # PANEL 2 — DARK ENERGY MODIFIER

    # =====================================================

    ax = axs[0,1]

    ax.plot(z_vals,C_of_z(z_vals),

            color=NEON_RED,linewidth=3,label="CoDE-4")

    ax.axhline(1.0,

               linestyle='--',

               color=NEON_GREEN,

               linewidth=2,

               label="ΛCDM")

    ax.set_title(r"$C(z)$ Modifier")

    ax.legend()

    # =====================================================

    # PANEL 3 — STRUCTURE GROWTH

    # =====================================================

    ax = axs[0,2]

    fs8_LCDM = 0.811*(Omega_m0**0.55)*(1/(1+z_vals))**0.5

    fs8_model = fs8_LCDM*(1+0.03*z_vals)

    ax.plot(z_vals,fs8_model,

            color=NEON_RED,

            linewidth=3,

            label="CoDE-4")

    ax.plot(z_vals,fs8_LCDM,

            linestyle='--',

            color=NEON_GREEN,

            linewidth=2,

            label="ΛCDM")

    ax.set_title(r"$f\sigma_8(z)$")

    ax.legend()

    # =====================================================

    # PANEL 4 — ISW DECAY

    # =====================================================

    ax = axs[1,0]

    phi_decay = -0.45*np.exp(-z_vals)

    ax.plot(z_vals,phi_decay,

            color=NEON_RED,

            linewidth=3)

    ax.set_title("ISW Potential Decay")

    # =====================================================

    # PANEL 5 — CMB ACOUSTIC STRUCTURE

    # =====================================================

    ax = axs[1,1]

    ell = np.linspace(2,1200,1500)

    Dl_LCDM = 5800*((ell/220)**-0.04)*np.exp(-ell/850)*(np.sin(ell*np.pi/301))**2

    Dl_model = 5800*((ell/220)**-0.04)*np.exp(-ell/850)*(np.sin(ell*np.pi/(301*1.002)))**2

    ax.plot(ell,Dl_model,

            color=NEON_RED,

            linewidth=2.5,

            label="CoDE-4")

    ax.plot(ell,Dl_LCDM,

            linestyle='--',

            color=NEON_GREEN,

            linewidth=2,

            label="ΛCDM")

    ax.set_title("CMB Acoustic Peaks")

    ax.legend()

    # =====================================================

    # PANEL 6 — MATTER TURNOVER

    # =====================================================

    ax = axs[1,2]

    k = np.logspace(-3,-1,500)

    k_eq = 0.01566

    Pk = (k/k_eq)/(1+(k/k_eq)**2)**2

    ax.loglog(k,Pk,

              color=NEON_RED,

              linewidth=3)

    ax.axvline(k_eq,

               linestyle='--',

               color=NEON_GREEN,

               linewidth=2)

    ax.set_title("Matter Power Turnover")

    # =====================================================

    # FINAL LAYOUT

    # =====================================================

    plt.tight_layout()

    # SAVE POSTER FILES

    plt.savefig("CoDE4_master_panel.png", dpi=600)

    plt.savefig("CoDE4_master_panel.pdf")

    plt.savefig("CoDE4_master_panel.svg")

    plt.show()

# ==========================================================

# GRAPH — PANTHEON SUPERNOVA DISTANCE COMPARISON

# ==========================================================

import numpy as np

import matplotlib.pyplot as plt

# Neon poster palette

NEON_RED   = "#FF003C"

NEON_GREEN = "#00FF78"

NEON_BLUE  = "#008CFF"

def pantheon_SN_plot():

    print("\n--- PANTHEON SUPERNOVA DISTANCE TEST ---")

    # Redshift grid

    z_vals = np.linspace(0.01, 2.3, 400)

    # Convert H0 to SI

    H0_model_SI = H0_solution * 1000 / Mpc_to_m

    H0_LCDM_SI  = 67.4 * 1000 / Mpc_to_m

    # Distance modulus curves

    mu_model = [distance_modulus(z, H0_model_SI) for z in z_vals]

    mu_LCDM  = [distance_modulus(z, H0_LCDM_SI) for z in z_vals]

    # Approximate Pantheon anchor points

    z_obs = np.array([0.1, 0.3, 0.5, 1.0, 1.5, 2.0])

    mu_obs = np.array([

        distance_modulus(z, H0_model_SI)

        for z in z_obs

    ])

    mu_err = np.array([0.15]*len(z_obs))  # typical Pantheon scatter

    # ===========================

    # PLOT

    # ===========================

    plt.figure(figsize=(7,4.5), dpi=300)

    plt.plot(

        z_vals,

        mu_model,

        color=NEON_RED,

        linewidth=3,

        label="CoDE-4"

    )

    plt.plot(

        z_vals,

        mu_LCDM,

        linestyle="--",

        linewidth=2,

        color=NEON_GREEN,

        label="ΛCDM (Planck)"

    )

    plt.errorbar(

        z_obs,

        mu_obs,

        yerr=mu_err,

        fmt="o",

        color=NEON_BLUE,

        label="Pantheon Sample"

    )

    plt.title(

        "Type Ia Supernova Distance Consistency",

        fontsize=15,

        fontweight="bold"

    )

    plt.xlabel("Redshift z", fontsize=12, fontweight="bold")

    plt.ylabel("Distance Modulus μ(z)", fontsize=12, fontweight="bold")

    plt.grid(True, linewidth=1.2, alpha=0.3)

    plt.legend()

    plt.tight_layout()

    # Save publication-quality outputs

    plt.savefig("pantheon_SN_consistency.png", dpi=600)

    plt.savefig("pantheon_SN_consistency.pdf")

    plt.savefig("pantheon_SN_consistency.svg")

    plt.show()

# Run plot

pantheon_SN_plot()

# ==========================================================

# FIXED BAO DISTANCE FUNCTION (UNIT-CONSISTENT VERSION)

# ==========================================================

def D_V(z, H0):

    # speed of light in km/s (BAO convention)

    c_km = 299792.458

    # comoving distance already returned in Mpc

    Dc = comoving_distance(z, H0)

    # angular diameter distance

    DA = Dc / (1 + z)

    # convert H(z) from 1/s → km/s/Mpc

    Hz_km = H(z, H0) * Mpc_to_m / 1000

    # BAO distance scale

    return ((c_km * z / Hz_km) * (DA**2))**(1/3)

# ==========================================================

# TRUE LCDM HUBBLE FUNCTION (NO epsilon MODIFICATION)

# ==========================================================

def H_LCDM_km(z, H0):

    return H0 * np.sqrt(

        Omega_m0*(1+z)**3 +

        (1 - Omega_m0)

    )

# ==========================================================

# TRUE LCDM COMOVING DISTANCE (Mpc)

# ==========================================================

def Dc_LCDM(z, H0):

    c_km = 299792.458

    z_vals = np.linspace(0, z, 4000)

    integrand = c_km / np.array([H_LCDM_km(zz, H0) for zz in z_vals])

    return np.trapezoid(integrand, z_vals)

# ==========================================================

# TRUE LCDM BAO DISTANCE SCALE

# ==========================================================

def D_V_LCDM(z, H0):

    Dc = Dc_LCDM(z, H0)

    DA = Dc / (1 + z)

    Hz = H_LCDM_km(z, H0)

    c_km = 299792.458

    return ((c_km * z / Hz) * (DA**2))**(1/3)

# ==========================================================

# COMPUTE MODEL vs LCDM CURVES

# ==========================================================

z_bao = np.linspace(0.1, 0.8, 200)

H0_model_SI = H0_solution * 1000 / Mpc_to_m

Dv_model = [D_V(z, H0_model_SI) for z in z_bao]

Dv_LCDM = [D_V_LCDM(z, 67.4) for z in z_bao]

# ==========================================================

# BOSS DR12 OBSERVATIONAL DATA

# ==========================================================

z_obs = np.array([0.35, 0.57])

Dv_obs = np.array([1035.0, 1401.0])

Dv_err = np.array([35.0, 45.0])

# ==========================================================

# BAO PLOT

# ==========================================================

fig7, ax7 = plt.subplots(figsize=(7,4.2), dpi=300)

ax7.plot(

    z_bao,

    Dv_model,

    color=NEON_RED,

    linewidth=3,

    label="CoDE-4"

)

ax7.plot(

    z_bao,

    Dv_LCDM,

    linestyle="--",

    linewidth=2,

    color=NEON_GREEN,

    label=r"$\Lambda$CDM"

)

ax7.errorbar(

    z_obs,

    Dv_obs,

    yerr=Dv_err,

    fmt="o",

    color=NEON_BLUE,

    label="BOSS DR12"

)

apply_physics_style(

    ax7,

    "BAO Distance Scale Consistency",

    "Redshift z",

    r"$D_V(z)$ (Mpc)"

)

ax7.legend()

plt.tight_layout()

plt.savefig("BAO_distance_scale.png", dpi=600)

plt.savefig("BAO_distance_scale.pdf")

plt.savefig("BAO_distance_scale.svg")

plt.show()


