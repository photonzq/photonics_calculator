# Photonics Calculator

Live page: **https://photonzq.github.io/photonics_calculator/**

`index.html` — a single self-contained page of interactive back-of-the-envelope
calculators for laser and integrated-photonics work. No build step and no server:
open it from the link above, or open the file straight from disk in a browser.
Tailwind and KaTeX load from a CDN, so keep a network connection for the styling
and equations.

## Sections

| # | Section | Calculators |
|---|---------|-------------|
| 1 | Spectral & Frequency Domain | wavelength ↔ frequency ↔ wavenumber ↔ photon energy (eV) with photon flux, spectral-band and telecom-band identification, wavelength presets applied to every calculator; differential bandwidth Δλ ↔ Δf ↔ Δṽ with the exact (non-linearised) band width |
| 2 | Ultrafast Pulses & Dispersion | transform-limited pulses; GDD ↔ total dispersion D; pulse elongation through a dispersive length (exact for Gaussian, exact numerical curve for sech²) |
| 3 | Fresnel Reflection & Interfaces | s- and p-polarized reflection and transmission vs. angle, with Brewster and critical angles; TIR phase shifts and evanescent depth |
| 4 | Cavity Resonators & Guided Waves | microresonator loss, Q and linewidth, bus coupling (κ_ex, Q_ex, loaded Q, regime, extinction, build-up, photon lifetime, ring radius) and Q₀/Q_ex extraction from a measured resonance; Fabry–Pérot cavity; waveguide Bragg grating with duty cycle and order, grating strength κ (from index contrast or entered), peak reflectivity, stop band, FWHM / null bandwidth, effective length and reflection spectrum; grating coupler (period ↔ fiber angle, radiating orders, angular dispersion, 1-/3-dB bandwidth, uniform-grating overlap, optimal scattering strength and CE estimate with spectrum) |
| 5 | Nonlinear Optics & Soliton Dynamics | soliton parameters and the characteristic lengths L_D and L_NL; Kerr microcomb / dissipative Kerr soliton (D₁, D₂, γ, parametric threshold, soliton existence range, duration, bandwidth, comb power and spectrum) |
| 6 | Gaussian Beam Propagation & Focusing | w(z), R(z), Gouy phase and on-axis intensity; thin-lens focusing; Gaussian mode / fiber coupling (size, axial, lateral and tilt mismatch; MFD or step-index fiber); draggable caustic diagram, transverse profile and lens views |

A global control sets the number of significant figures shown everywhere, and
every calculator has a collapsible "Physics notes & assumptions" panel that
states its conventions and range of validity.

## Repository history

This repo starts from a snapshot of the three HTML files as they stood before
any consolidation, so the two that were retired afterwards remain recoverable:

```bash
# earlier single-purpose nm/THz converter
git show archive/wavelength-converter-old:wavelength_frequency_converter_old.html > old_converter.html

# standalone Gaussian beam tool (now Section 6)
git show archive/gaussian-beam-standalone:gaussian_beam_calculator.html > gaussian_beam_calculator.html
```

`git log --follow index.html` follows the page across its renames from
`wavelength_frequency_converter.html`.

## Hosting

GitHub Pages serves this from the `main` branch, root folder — the page is
`index.html`, so the repository URL loads the calculator directly. `.nojekyll`
is present so the files are published verbatim rather than run through Jekyll.
Every asset is either inline or an absolute `https://` CDN URL and no path is
root-relative, so the page works unchanged at a project subpath, at a user page,
or behind a custom domain.

## Note on Google Drive

The working copy lives inside a synced Drive folder, and Drive syncs `.git/`
like any other directory, so avoid running git here from two machines at the same
time. The canonical copy is the GitHub remote:

```bash
git clone https://github.com/photonzq/photonics_calculator.git
```
