# CoDE-4 Theoretical Framework

## Overview

CoDE-4 (Cosmological Dynamical Evolution Framework – Version 4) is a phenomenological late-time cosmological framework designed to introduce a bounded perturbative deformation to the standard ΛCDM expansion history while preserving near-standard early-universe behavior.

The framework was constructed with the following primary objectives:

- preserve recombination-era geometric consistency,
- maintain near-ΛCDM early-universe recovery,
- allow mild late-time expansion deformation,
- preserve acoustic structure,
- and explore weakly dynamical effective dark-energy evolution.

Rather than replacing standard cosmology at the fundamental level, CoDE-4 currently functions as an exploratory phenomenological framework for studying perturbative late-time modifications to cosmological expansion dynamics.

---

# Core Postulates

## 1. Near-ΛCDM Early-Universe Recovery

Any late-time deformation must asymptotically vanish at high redshift in order to preserve:
- recombination consistency,
- acoustic peak structure,
- matter-radiation equality,
- and CMB geometric constraints.

---

## 2. Bounded Perturbative Expansion Deformation

The modification must remain bounded across cosmic evolution and avoid:
- singularities,
- runaway growth,
- or unstable ultraviolet behavior.

---

## 3. Smooth Late-Time Enhancement

The framework allows mild late-time enhancement of the effective expansion sector in order to explore potential alleviation of the Hubble tension.

---

## 4. Recursive Acoustic Consistency

The deformation is constrained such that the acoustic structure of the universe remains approximately preserved despite late-time vacuum evolution.

---

# Modified Expansion Structure

The CoDE-4 expansion framework modifies the effective vacuum sector through a bounded redshift-dependent correction function:

$$
C(z)=1+\epsilon\left(\frac{z}{(1+z)^2}\right)\left(1-\frac{0.15}{1+z}\right)
$$

This structure was constructed phenomenologically to satisfy several desired asymptotic conditions:

- low-redshift enhancement,
- high-redshift suppression,
- smooth bounded evolution,
- and near-standard early-universe recovery.

The modified Friedmann structure becomes:

$$
H^2(z)=H_0^2\left[
\Omega_m(1+z)^3
+\Omega_r(1+z)^4
+\Omega_\Lambda C(z)
\right]
$$

where:
- $\Omega_m$ = matter density parameter,
- $\Omega_r$ = radiation density parameter,
- $\Omega_\Lambda$ = effective vacuum density parameter,
- $C(z)$ = perturbative deformation sector.

---

# Effective Dark-Energy Interpretation

The modification may be interpreted phenomenologically as a weakly dynamical effective dark-energy sector:

$$
\rho_X(z)=\rho_\Lambda C(z)
$$

leading to an effective equation-of-state evolution:

$$
w_{\mathrm{eff}}(z)
=
-1+\frac{1+z}{3C(z)}\frac{dC}{dz}
$$

The resulting evolution remains close to ΛCDM while allowing mild late-time deviations from a pure cosmological constant.

---

# Asymptotic Behavior

The framework was designed to recover near-standard cosmological evolution at high redshift.

As $z \to \infty$:
- $C(z) \to 1$
- $w_{\mathrm{eff}}(z) \to -1$
- the expansion history approaches ΛCDM behavior.

This asymptotic recovery is essential for preserving:
- recombination physics,
- early structure formation,
- and CMB geometric consistency.

---

# Stability Considerations

The perturbative structure of CoDE-4 remains bounded throughout cosmic evolution and avoids:
- finite-redshift singularities,
- divergent vacuum behavior,
- and unstable asymptotic growth.

The deformation remains smooth across both low-redshift and high-redshift limits.

---

# Perturbation Interpretation

The framework also admits a perturbative growth interpretation through modified expansion evolution.

The linear matter perturbation equation is written as:

$$
\delta''
+
\left(
\frac{3}{a}
+
\frac{H'}{H}
\right)\delta'
-
\frac{3\Omega_mH_0^2}{2a^5H^2}\delta
=
0
$$

This allows exploration of:
- structure growth,
- clustering evolution,
- growth-index behavior,
- and weak-lensing consistency.

---

# CPL Mapping

The effective equation-of-state evolution may be approximately mapped into CPL form:

$$
w(a)=w_0+w_a(1-a)
$$

with current CoDE-4 estimates yielding approximately:
- $w_0 \approx -0.986$
- $w_a \approx -0.052$

This places the framework within weakly dynamical dark-energy territory while remaining close to ΛCDM behavior.

---

# Phenomenological Interpretation

Overall, CoDE-4 behaves as:
- a weakly dynamical late-time cosmological framework,
- with near-ΛCDM early-universe recovery,
- perturbative late-time vacuum evolution,
- preserved acoustic geometry,
- and bounded asymptotic structure.

The framework currently remains exploratory and phenomenological in nature.

---

# Current Limitations

The current implementation does not yet include:
- full covariance likelihood analyses,
- Bayesian evidence comparison,
- Cobaya MCMC parameter estimation,
- or CLASS/CAMB perturbation integration.

These remain important future directions for further validation and refinement.

---

# Future Directions

Planned future extensions include:
- covariance-aware χ² analysis,
- full MCMC pipelines,
- CLASS/CAMB integration,
- expanded perturbation evolution,
- Bayesian comparison against ΛCDM,
- and larger observational consistency studies.
