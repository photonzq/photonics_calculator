# Photonics Calculator

`photonics_calculator.html` — a single self-contained page of interactive
back-of-the-envelope calculators for laser and integrated-photonics work.
No build step and no server: open the file in a browser (Tailwind and KaTeX
load from a CDN, so keep a network connection for the styling and equations).

## Sections

| # | Section | Calculators |
|---|---------|-------------|
| 1 | Spectral & Frequency Domain | wavelength ↔ frequency ↔ wavenumber with spectral-band and telecom-band identification; differential bandwidth Δλ ↔ Δf ↔ Δṽ |
| 2 | Ultrafast Pulses & Dispersion | transform-limited pulses; GDD ↔ total dispersion D; pulse elongation through a dispersive length |
| 3 | Fresnel Reflection & Interfaces | s- and p-polarized reflection and transmission vs. angle, with Brewster and critical angles |
| 4 | Cavity Resonators & Guided Waves | microresonator loss, Q and linewidth; Fabry–Pérot cavity; waveguide Bragg grating |
| 5 | Nonlinear Optics & Soliton Dynamics | soliton parameters and the characteristic lengths L_D and L_NL |
| 6 | Gaussian Beam Propagation & Focusing | w(z), R(z), Gouy phase and on-axis intensity; thin-lens focusing; draggable caustic diagram, transverse profile and lens views |

A global precision control sets the displayed digits everywhere, and the
Gaussian beam section can pull its wavelength from Section 1.

## Repository history

This repo starts from a snapshot of the three HTML files as they stood before
any consolidation, so the two that were retired afterwards remain recoverable:

```bash
# earlier single-purpose nm/THz converter
git show archive/wavelength-converter-old:wavelength_frequency_converter_old.html > old_converter.html

# standalone Gaussian beam tool (now Section 6)
git show archive/gaussian-beam-standalone:gaussian_beam_calculator.html > gaussian_beam_calculator.html
```

`git log --follow photonics_calculator.html` follows the file across its rename
from `wavelength_frequency_converter.html`.

## Note on Google Drive

The working copy lives inside a synced Drive folder, and Drive syncs `.git/`
like any other directory, so avoid running git here from two machines at the same
time. The canonical copy is the GitHub remote:

```bash
git clone https://github.com/photonzq/photonics_calculator.git
```
