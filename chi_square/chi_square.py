# ============================================================

# CoDE-4 COMPRESSED OBSERVABLE χ² ANALYSIS

# ============================================================

import numpy as np

print("\n" + "="*70)

print("CoDE-4 COMPRESSED OBSERVABLE χ² ANALYSIS")

print("="*70)

# ============================================================

# OBSERVATIONAL VALUES

# ============================================================

H0_obs = 73.04

sigma_H0 = 1.04

R_obs = 1.7502

sigma_R = 0.0046

l1_obs = 220.6

sigma_l1 = 0.5

S8_obs = 0.776

sigma_S8 = 0.017

z_eq_obs = 3400

sigma_z_eq = 100

# ============================================================

# CoDE-4 MODEL VALUES

# ============================================================

H0_model = 72.86

R_model = 1.75002

l1_model = 220.72

S8_model = 0.8396

z_eq_model = 3499

# ============================================================

# CHI-SQUARE FUNCTION

# ============================================================

def chi2(model, obs, sigma):

    return ((model - obs) ** 2) / (sigma ** 2)

# ============================================================

# INDIVIDUAL χ² VALUES

# ============================================================

chi2_H0 = chi2(H0_model, H0_obs, sigma_H0)

chi2_R = chi2(R_model, R_obs, sigma_R)

chi2_l1 = chi2(l1_model, l1_obs, sigma_l1)

chi2_S8 = chi2(S8_model, S8_obs, sigma_S8)

chi2_eq = chi2(z_eq_model, z_eq_obs, sigma_z_eq)

# ============================================================

# TOTAL χ²

# ============================================================

chi2_total = (

    chi2_H0 +

    chi2_R +

    chi2_l1 +

    chi2_S8 +

    chi2_eq

)

# ============================================================

# OUTPUT

# ============================================================

print("\n--- INDIVIDUAL χ² CONTRIBUTIONS ---")

print(f"H0 χ²              = {chi2_H0:.4f}")

print(f"Shift Parameter χ² = {chi2_R:.4f}")

print(f"Acoustic Peak χ²   = {chi2_l1:.4f}")

print(f"S8 χ²              = {chi2_S8:.4f}")

print(f"z_eq χ²            = {chi2_eq:.4f}")

print("\n--- TOTAL χ² ---")

print(f"Total χ² = {chi2_total:.4f}")

print("\n--- INTERPRETATION ---")

if chi2_total < 5:

    print("Fit Status: EXCELLENT compressed-observable consistency")

elif chi2_total < 15:

    print("Fit Status: GOOD compressed-observable consistency")

elif chi2_total < 30:

    print("Fit Status: MODERATE consistency")

else:

    print("Fit Status: Significant observational tension")

print("="*70)
