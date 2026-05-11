# Modifier Function Derivation

## Overview

The central mathematical structure of the CoDE-4 framework is the bounded deformation function:

$$
C(z) = 1 + \epsilon \left(\frac{z}{(1+z)^2}\right)\left(1-\frac{0.15}{1+z}\right)
$$

This function acts as the primary perturbative modification to the effective vacuum sector of the cosmological expansion framework.

The modifier was intentionally constructed to satisfy several stability and phenomenological conditions simultaneously:
- bounded late-time deformation,
- asymptotic recovery of standard cosmology,
- finite perturbative behavior,
- weak low-redshift modification,
- and stable high-redshift suppression.

The function therefore acts as:
- a phenomenological expansion deformation,
rather than:
- a strongly divergent cosmological correction.

---

# Construction Philosophy

The modifier was not derived from a fundamental field-theoretic Lagrangian.

Instead, the function was phenomenologically constructed to satisfy:
- cosmological boundedness conditions,
- asymptotic recovery requirements,
- observational smoothness,
- and perturbative stability constraints.

The construction philosophy was guided by:
- preservation of early-universe consistency,
- weak late-time expansion modification,
- and numerical stability throughout cosmological evolution.

---

# Primary Deformation Structure

The dominant perturbative component is:

$$
\frac{z}{(1+z)^2}
$$

This term was selected because it naturally satisfies several desirable asymptotic properties.

At low redshift:
- the correction becomes non-negligible,
- allowing weak late-time expansion deformation.

At high redshift:

$$
z \to \infty
$$

the term satisfies:

$$
\frac{z}{(1+z)^2} \to 0
$$

This produces asymptotic suppression of the deformation sector.

---

# Low-Redshift Expansion

For small redshift fields near the current cosmic epoch ($z \ll 1$), the full joint modifier function $C(z)$ is analytic and can be modeled via Taylor-series expansion to describe the immediate dark energy trajectory:

$$
C(z) \approx 1 + \epsilon \left( 0.85z - 1.55z^2 + 2.20z^3 \right) + \mathcal{O}(z^4)
$$

This demonstrates:
- smooth perturbative low-redshift behavior,
- finite correction amplitudes,
- and analytic Taylor-expandable structure.

The modifier therefore behaves as:
- a weak perturbative correction near the present cosmological epoch.

---

# High-Redshift Suppression

At sufficiently large redshift:

$$
z \to \infty
$$

the dominant deformation satisfies:

$$
\frac{z}{(1+z)^2} \sim \frac{1}{z} \to 0
$$

This suppression is one of the central stability properties of the framework because it ensures:
- recovery of radiation domination,
- preservation of recombination-era physics,
- stable Big Bang Nucleosynthesis evolution,
- and near-standard early-universe expansion.

The deformation therefore does not dominate asymptotic cosmological evolution.

---

# Secondary Suppression Structure

The modifier additionally contains the correction factor:

$$
\left( 1-\frac{0.15}{1+z} \right)
$$

This term acts as:
- a bounded amplitude regulator,
- a late-time suppression stabilizer,
- and a smooth correction-shaping component.

At high redshift:

$$
\frac{0.15}{1+z} \to 0
$$

which preserves asymptotic recovery.

At low redshift:
- the correction weakly reshapes the deformation profile,
- while preserving bounded evolution.

The secondary suppression structure therefore helps maintain:
- smooth perturbative behavior,
- finite late-time evolution,
- and numerical stability.

---

# Analytical Smooth Differentiability

To avoid numerical finite-difference calculation errors inside the cosmological continuity solvers, the framework utilizes the exact analytical derivative of the modifier profile with respect to redshift:

$$
\frac{dC}{dz} = \epsilon \left[ \frac{1-z}{(1+z)^3} \left( 1 - \frac{0.15}{1+z} \right) + \left( \frac{z}{(1+z)^2} \right) \left( \frac{0.15}{(1+z)^2} \right) \right]
$$

This derivative directly feeds into the reconstruction of the effective dark energy equation of state:

$$
w_{\mathrm{eff}}(z) = -1 + \frac{1+z}{3C(z)} \frac{dC}{dz}
$$

---

# Effective Vacuum Modification

The deformation function modifies the effective vacuum sector through:

$$
\rho_X(z) = \rho_\Lambda C(z)
$$

leading to the modified expansion equation:

$$
H^2(z) = H_0^2 \left[ \Omega_m(1+z)^3 + \Omega_r(1+z)^4 \left( 1+\frac{\epsilon z}{(1+z)^3} \right) + \Omega_\Lambda C(z) \right]
$$

The deformation therefore acts primarily through:
- bounded late-time vacuum evolution,
while preserving:
- near-standard matter scaling,
- stable radiation behavior,
- and asymptotic early-universe recovery.

---

# Boundedness Properties

An important construction requirement was prevention of divergent cosmological evolution.

The modifier satisfies:
- finite low-redshift behavior,
- finite intermediate-redshift behavior,
- and asymptotic high-redshift suppression.

The framework therefore does not currently exhibit:
- runaway vacuum growth,
- asymptotic divergence,
- or singular perturbative behavior.

The modifier remains bounded throughout the tested cosmological range.

---

# Numerical Stability Role

The modifier additionally plays an important role in numerical stability.

Because the correction:
- remains smooth,
- remains differentiable,
- and remains asymptotically suppressed,

the resulting numerical evolution remains:
- stable,
- finite,
- and non-divergent across the tested cosmological domain.

This helps preserve:
- perturbation integration stability,
- recombination consistency,
- and smooth background evolution.

---

# Phenomenological Interpretation

Overall, the CoDE-4 modifier function behaves as:
- a bounded phenomenological deformation of late-time cosmological expansion,
- with asymptotic suppression,
- stable perturbative structure,
- and preserved early-universe recovery.

The modifier was intentionally constructed to produce:
- weak late-time cosmological deviation,
while avoiding:
- strong high-redshift disruption,
- divergent perturbative behavior,
- and unstable cosmological evolution.

The framework therefore currently behaves as:
- a weakly dynamical phenomenological expansion-deformation cosmology framework.
