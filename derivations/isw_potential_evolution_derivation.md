# ISW Potential Evolution Derivation

## Overview

One of the most important late-time cosmological consistency diagnostics is the Integrated Sachs-Wolfe (ISW) effect.

The ISW sector probes:
- late-time gravitational potential evolution,
- dark-energy-driven expansion dynamics,
- perturbative stability,
- and large-scale cosmological structure evolution.

Because the CoDE-4 framework modifies the late-time cosmological expansion sector through a bounded perturbative deformation, the framework naturally alters:
- gravitational potential decay,
- late-time metric evolution,
- and ISW anisotropy behavior.

The purpose of this derivation is therefore to:
- reconstruct the effective ISW potential evolution,
- analyze late-time gravitational decay behavior,
- evaluate perturbative stability,
- and connect the framework directly to observational ISW consistency diagnostics.

---

# Physical Origin of the ISW Effect

The Integrated Sachs-Wolfe effect arises because:
- photons propagating through evolving gravitational potentials experience net energy shifts,
- and time-dependent metric evolution produces additional CMB anisotropies.

In a perfectly matter-dominated universe:
- gravitational potentials remain approximately constant,
- and the ISW contribution becomes negligible.

However:
- accelerated expansion causes gravitational potentials to decay,
- producing measurable late-time ISW signatures.

The ISW effect therefore acts as:
- a direct probe of dark-energy-driven expansion,
- and a sensitive stability diagnostic for modified cosmological models.

---

# Gravitational Potential Reconstruction

Within the CoDE-4 framework, the normalized linear gravitational potential on sub-horizon scales is reconstructed phenomenologically through:

$$
\Phi(a) = \frac{\delta(a)}{a}
$$

where:
- $\delta(a)$ is the linear matter perturbation amplitude,
- and $a$ is the cosmological scale factor.

This reconstruction corresponds directly to the normalized potential-evolution structure used inside the CoDE-4 numerical engine. The framework therefore tracks relative gravitational well evolution rather than an absolute gauge-fixed relativistic metric potential.

---

# Perturbation Evolution Background

The perturbation sector evolves with respect to scale factor $a$ (where $' \equiv d/da$) through:

$$
\delta'' + \left( \frac{3}{a} + \frac{H'}{H} \right)\delta' - \frac{3 H_0^2 \Omega_m}{2 a^5 H^2} \left( 1 + \beta f(z) \right) \delta = 0
$$

where:
- the background expansion history modifies perturbation friction,
- and the $\beta$ sector modifies effective clustering behavior.

The resulting perturbation evolution directly feeds into the gravitational potential reconstruction, and therefore into ISW evolution.

---

# Potential Evolution Structure

The ISW effect depends on:
- the time evolution of the gravitational potential,
- rather than the potential amplitude itself.

The CoDE-4 numerical engine evaluates this behavior through:

$$
\frac{d\Phi}{d\ln a}
$$

The sign of this quantity determines whether gravitational potentials decay, remain constant, or grow dynamically.

---

# Standard Cosmological Expectation

In standard late-time dark-energy-dominated cosmology:
- accelerated expansion suppresses clustering growth,
- gravitational wells become shallower,
- and gravitational potentials decay.

Therefore:

$$
\frac{d\Phi}{d\ln a} < 0
$$

acts as the expected stable late-time behavior. A positive value would instead indicate growing gravitational potentials, unstable clustering amplification, or non-standard perturbative dynamics.

---

# Numerical ISW Reconstruction

The CoDE-4 numerical engine reconstructs normalized perturbation growth, gravitational potential evolution, and ISW decay behavior directly. Current numerical reconstruction produces approximately:
- $\frac{d\Phi}{d\ln a}(z=0) \approx -0.4705$
- $\frac{d\Phi}{d\ln a}(z=0.5) \approx -0.2719$
- $\frac{d\Phi}{d\ln a}(z=1.0) \approx -0.1446$
- $\frac{d\Phi}{d\ln a}(z=2.0) \approx -0.0465$

These results indicate:
- stable late-time potential decay,
- bounded perturbative behavior,
- and standard dark-energy-dominated ISW evolution.

---

# Physical Interpretation

The reconstructed negative potential evolution implies:
- gravitational wells become progressively shallower,
- structure growth becomes increasingly suppressed,
- and accelerated expansion dominates late-time dynamics.

The framework therefore preserves:
- standard qualitative ISW behavior,
- stable metric evolution,
- and bounded perturbative cosmological dynamics.

---

# High-Redshift Recovery

At sufficiently large redshift:

$$
z \to \infty
$$

the deformation satisfies:

$$
C(z) \to 1
$$

while the perturbative correction satisfies:

$$
\frac{\epsilon z}{(1+z)^3} \to 0
$$

The framework therefore asymptotically approaches:
- near-standard matter domination,
- stable early perturbation evolution,
- and approximately constant gravitational potentials.

This recovery is important because the ISW effect should become suppressed during matter domination, preserving standard early-universe cosmological behavior.

---

# Stability Interpretation

The ISW sector acts as one of the strongest dynamical stability diagnostics of the framework. Current numerical evolution demonstrates:
- no growing late-time potentials,
- no runaway perturbative amplification,
- and no unstable gravitational clustering behavior.

The framework therefore remains:
- dynamically stable,
- perturbatively bounded,
- and numerically well-behaved throughout tested cosmological evolution.

---

# Connection to Observations

The ISW sector connects directly to:
- large-angle CMB anisotropies,
- late-time structure evolution,
- weak gravitational lensing,
- and galaxy-CMB cross-correlation measurements.

Because the CoDE-4 framework preserves negative potential decay, stable perturbative evolution, and near-standard late-time gravitational dynamics, the framework remains phenomenologically compatible with standard late-time ISW behavior and stable cosmological metric evolution.

---

# Phenomenological Interpretation

Overall, the CoDE-4 framework reconstructs:
- a stable late-time ISW sector,
- bounded gravitational potential evolution,
- weak perturbative deformation,
- and standard dark-energy-dominated potential decay behavior.

The framework therefore currently behaves as:
- a phenomenologically stable modified-expansion cosmology framework,
- with weak late-time perturbative deformation,
- while preserving stable gravitational potential evolution across cosmological history.
