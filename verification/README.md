# Verification

`reference_checks.py` recomputes quantities the page shows, by a different
numerical route where possible, for fixed inputs that can be typed into the page.

```bash
python reference_checks.py      # needs numpy and scipy; runs in a few seconds
```

| # | Check | Independent route | Page value | Reference |
|---|-------|-------------------|-----------|-----------|
| 1 | Depleted SHG (1550 nm, defaults, 1 W) | adaptive ODE (DOP853) vs. $P_1\tanh^2(\sqrt{\eta P_1}L)$ | 961.7 mW | 961.715 mW |
| 2 | SFG 1550 + 1310 nm, $\Delta k_r L/2 = 1.67$, 500 + 300 mW | three-wave ODE, Manley–Rowe check | 108.2 µW | 108.171 µW |
| 3 | DFG/OPA 1550 / 780 nm, $\Delta k_r L/2 = 0.79$, 800 mW pump, 5 mW signal | three-wave ODE; textbook $\Gamma^2$ | 55.27 mW idler, 10.9 dB (11.2 dB small-signal) | 55.266 mW, 10.862 dB (11.168 dB) |
| 4 | OPA gain, phase-matched, 1560 / 775 nm, 500 mW | ODE vs. $\cosh^2(\Gamma L)$ | 8.52 dB | 8.519 dB |
| 5 | QPM acceptance with GVD ($n_{g1} = n_{g2}$, $\beta_2 = 100$ fs²/mm) | brute-force FWHM of $\operatorname{sinc}^2$ on a fine grid | 23.7 THz | 23.748 THz |
| 6 | Ring OPO threshold (defaults) | time-domain three-mode coupled-mode run at 0.9 / 1.1 / 1.5 $P_{\text{th}}$ | 490.5 nW, $g_0/2\pi$ = 2.562 MHz | 490.55 nW, 2.5624 MHz; seed decays below, grows above |
| 7 | Ring spectrum with drop port and backscatter (20/20/15 MHz, $g/2\pi$ = 40 MHz) | 2×2 linear solve; energy balance; limits | T 58.76 %, D 4.086 % + 8.646 % bwd, R 11.53 % | identical; balance closes to $4\times10^{-16}$ |
| 8 | Resonant EO comb slope ($\mathcal{F}_L$ = 2500, $\beta$ = 1) | FFT of $1/(1 - r e^{i\beta\sin\theta})$ | 0.0109 dB/line | 0.01091 dB/line |

Inputs not listed are the page defaults. For the QPM rows, set the process,
wavelengths, effective indices, powers and (for rows 2 and 3) the poling period
exactly as given in the script's section headers.
