"""Independent reference calculations for the calculator page (index.html).

Each section re-derives a quantity the page computes, using a different numerical
route where possible (adaptive ODE solvers, linear solves, FFTs), and prints the
value the page should show for the stated inputs. Run:

    python reference_checks.py

Requires numpy and scipy. Runtime is a few seconds (the ring-OPO time-domain
runs dominate).
"""
import numpy as np
from scipy.integrate import solve_ivp

c = 299792458.0
eps0 = 8.8541878128e-12
hbar = 1.054571817e-34


def header(title):
    print("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78)


# ---------------------------------------------------------------------------
# QPM: normalized efficiencies (same definitions as the page)
#   d_eff = (2 / (m pi)) d33 |sin(m pi D)|
#   eta_SHG = 8 pi^2 d_eff^2 / (eps0 c A n1^2 n3 lambda1^2)
#   eta_SFG(lambda3) = 8 pi^2 d_eff^2 / (eps0 c A n1 n2 n3 lambda3^2)
#   eta_DFG = 8 pi^2 d_eff^2 / (eps0 c A n1 n2 n3 lambda2^2)   (idler generated)
# ---------------------------------------------------------------------------
def d_eff(d33, m=1, duty=0.5):
    return d33 * 2 / (m * np.pi) * abs(np.sin(m * np.pi * duty))


def eta_shg(deff, A, n1, n3, l1):
    return 8 * np.pi ** 2 * deff ** 2 / (eps0 * c * A * n1 ** 2 * n3 * l1 ** 2)


def eta_three(deff, A, n1, n2, n3, lam):
    return 8 * np.pi ** 2 * deff ** 2 / (eps0 * c * A * n1 * n2 * n3 * lam ** 2)


def shg_ode(P1, eta, L, dk):
    """dA2/dz = -i k A1^2 e^{-i dk z}, dA1/dz = -i k A2 A1* e^{+i dk z}, k^2 = eta, |A|^2 = P."""
    k = np.sqrt(eta)

    def f(z, y):
        A1 = y[0] + 1j * y[1]
        A2 = y[2] + 1j * y[3]
        d2 = -1j * k * A1 ** 2 * np.exp(-1j * dk * z)
        d1 = -1j * k * A2 * np.conj(A1) * np.exp(1j * dk * z)
        return [d1.real, d1.imag, d2.real, d2.imag]

    s = solve_ivp(f, [0, L], [np.sqrt(P1), 0, 0, 0], rtol=1e-11, atol=1e-14, method="DOP853")
    y = s.y[:, -1]
    return abs(y[2] + 1j * y[3]) ** 2, abs(y[0] + 1j * y[1]) ** 2


def three_wave_ode(P, lams, eta_sfg3, L, dkr):
    """Three-wave mixing in b_j = A_j / sqrt(omega_j) (|A_j|^2 = P_j), G^2 = eta_SFG w1 w2 / w3.
    db1/dz = -i G b3 b2* e^{i dk z}, db2/dz = -i G b3 b1* e^{i dk z}, db3/dz = -i G b1 b2 e^{-i dk z}."""
    w = [2 * np.pi * c / lam for lam in lams]
    G = np.sqrt(eta_sfg3 * w[0] * w[1] / w[2])

    def f(z, y):
        b1 = y[0] + 1j * y[1]
        b2 = y[2] + 1j * y[3]
        b3 = y[4] + 1j * y[5]
        e = np.exp(1j * dkr * z)
        d1 = -1j * G * b3 * np.conj(b2) * e
        d2 = -1j * G * b3 * np.conj(b1) * e
        d3 = -1j * G * b1 * b2 * np.conj(e)
        return [d1.real, d1.imag, d2.real, d2.imag, d3.real, d3.imag]

    y0 = []
    for Pj, wj in zip(P, w):
        y0 += [np.sqrt(Pj / wj), 0.0]
    s = solve_ivp(f, [0, L], y0, rtol=1e-12, atol=1e-22, method="DOP853")
    y = s.y[:, -1]
    return [abs(y[2 * i] + 1j * y[2 * i + 1]) ** 2 * w[i] for i in range(3)], G, w


def residual_mismatch(lams, ns, period, m=1):
    l1, l2, l3 = lams
    n1, n2, n3 = ns
    dk = 2 * np.pi * (n3 / l3 - n1 / l1 - n2 / l2)
    return dk - np.sign(dk) * 2 * np.pi * m / period, dk


def sinc2(x):
    return 1.0 if abs(x) < 1e-12 else (np.sin(x) / x) ** 2


# Page defaults shared by the QPM cases
D33, AEFF, LQ = 25e-12, 2e-12, 5e-3
DEFF = d_eff(D33)

header("1. Depleted SHG versus P2 = P1 tanh^2(sqrt(eta P1) L)  (page: SHG, 1550 nm, n1 = 1.88, n3 = 2.07)")
eta = eta_shg(DEFF, AEFF, 1.88, 2.07, 1550e-9)
print(f"eta_SHG = {eta:.6e} 1/(W m^2) = {eta * 1e-4 * 100:.2f} %/W/cm^2   (page: 2143 %/W/cm^2)")
for P1 in (0.01, 0.2, 1.0):
    p2, p1 = shg_ode(P1, eta, LQ, 0.0)
    law = P1 * np.tanh(np.sqrt(eta * P1) * LQ) ** 2
    print(f"  P1 = {P1:5.2f} W: ODE P2 = {p2:.6e} W, tanh^2 law = {law:.6e} W, P1+P2 = {p1 + p2:.9f} W")
print("  page at P1 = 10 mW: 517.2 uW; at 1 W: 961.7 mW")

header("2. SFG with residual mismatch  (page: SFG 1550 + 1310 nm, n = 1.88/1.90/2.08, Lambda = 3.7547312919999998 um, 500 + 300 mW)")
l1, l2 = 1550e-9, 1310e-9
l3 = 1 / (1 / l1 + 1 / l2)
ns = (1.88, 1.90, 2.08)
dkr, dk = residual_mismatch((l1, l2, l3), ns, 3.7547312919999998e-6)
e3 = eta_three(DEFF, AEFF, *ns, l3)
P, _, w = three_wave_ode([0.5, 0.3, 0.0], (l1, l2, l3), e3, LQ, dkr)
und = e3 * 0.5 * 0.3 * LQ ** 2 * sinc2(dkr * LQ / 2)
print(f"  dk_r L/2 = {dkr * LQ / 2:.4f}, undepleted P3 = {und:.4e} W (page 1.343 W)")
print(f"  coupled-mode P3 = {P[2]:.5e} W (page 108.2 uW); Manley-Rowe P1/w1 + P3/w3 = {P[0] / w[0] + P[2] / w[2]:.6e} (initial {0.5 / w[0]:.6e})")

header("3. DFG / OPA with residual mismatch  (page: DFG signal 1550 nm, pump 780 nm, Lambda = 4.10472078 um, pump 800 mW, signal 5 mW)")
l1, l3 = 1550e-9, 780e-9
l2 = 1 / (1 / l3 - 1 / l1)
dkr, dk = residual_mismatch((l1, l2, l3), ns, 4.10472078e-6)
e3 = eta_three(DEFF, AEFF, *ns, l3)
e2 = eta_three(DEFF, AEFF, *ns, l2)
P, G, w = three_wave_ode([0.005, 0.0, 0.8], (l1, l2, l3), e3, LQ, dkr)
und = e2 * 0.8 * 0.005 * LQ ** 2 * sinc2(dkr * LQ / 2)
Gam2 = G ** 2 * 0.8 / w[2]
g2 = Gam2 - (dkr / 2) ** 2
small = 1 + (Gam2 / g2 * np.sinh(np.sqrt(g2) * LQ) ** 2 if g2 > 0 else Gam2 / -g2 * np.sin(np.sqrt(-g2) * LQ) ** 2)
print(f"  dk_r L/2 = {dkr * LQ / 2:.4f}, undepleted idler = {und:.4e} W (page 16.6 mW)")
print(f"  coupled-mode idler = {P[1]:.5e} W (page 55.27 mW), signal gain = {10 * np.log10(P[0] / 0.005):.3f} dB (page 10.9 dB)")
print(f"  small-signal gain 1 + (Gamma^2/g^2) sinh^2(gL) = {10 * np.log10(small):.3f} dB (page 11.2 dB)")
Gam2_std = 8 * np.pi ** 2 * DEFF ** 2 * 0.8 / (eps0 * c * ns[0] * ns[1] * ns[2] * AEFF * l1 * l2)
print(f"  Gamma^2 from G^2 P3/w3 = {Gam2:.6e}; textbook 8 pi^2 d^2 P3 / (eps0 c n1 n2 n3 A l1 l2) = {Gam2_std:.6e}")

header("4. Small-signal OPA gain, phase-matched, near degeneracy  (page: DFG 1560 / 775 nm, n = 1.88/1.88/2.07, 500 mW pump, 1 uW signal)")
l1, l3 = 1560e-9, 775e-9
l2 = 1 / (1 / l3 - 1 / l1)
ns2 = (1.88, 1.88, 2.07)
e3 = eta_three(DEFF, AEFF, *ns2, l3)
P, G, w = three_wave_ode([1e-6, 0.0, 0.5], (l1, l2, l3), e3, LQ, 0.0)
Gam = np.sqrt(G ** 2 * 0.5 / w[2])
print(f"  eta_SFG(lambda3) = {e3 * 1e-4 * 100:.1f} %/W/cm^2 (page 8573)")
print(f"  ODE gain = {10 * np.log10(P[0] / 1e-6):.4f} dB, cosh^2(Gamma L) = {10 * np.log10(np.cosh(Gam * LQ) ** 2):.4f} dB (page 8.52 dB)")

header("5. QPM acceptance with group-velocity dispersion: FWHM of sinc^2[(a W + b W^2) L / 2]")
L = 5e-3
for a, b, label in ((0.0, -1e-25, "a = 0, b = -1e-25 s^2/m (GV-matched DFG, beta2 = 100 fs^2/mm each)"),
                    (0.15 / c, 0.0, "a = 0.15/c, b = 0 (first order)"),
                    (0.15 / c, -1e-25, "a = 0.15/c, b = -1e-25 s^2/m")):
    W = np.linspace(-3e14, 3e14, 2_000_001)
    eff = np.sinc((a * W + b * W ** 2) * L / 2 / np.pi) ** 2   # numpy sinc(x) = sin(pi x)/(pi x)
    i0 = np.argmin(np.abs(W))
    lo = i0
    while lo > 0 and eff[lo] >= 0.5:
        lo -= 1
    hi = i0
    while hi < len(W) - 1 and eff[hi] >= 0.5:
        hi += 1
    fwhm = (W[hi] - W[lo]) / (2 * np.pi)
    D = 2 * 1.391557 / L
    if b == 0:
        closed = 2 * D / abs(a) / (2 * np.pi)
    elif a == 0:
        closed = 2 * np.sqrt(D / abs(b)) / (2 * np.pi)
    else:
        roots = []
        for rhs in (D, -D):
            disc = a * a + 4 * b * rhs
            if disc >= 0:
                roots += [(-a + np.sqrt(disc)) / (2 * b), (-a - np.sqrt(disc)) / (2 * b)]
        closed = (min(r for r in roots if r > 0) - max(r for r in roots if r < 0)) / (2 * np.pi)
    print(f"  {label}: numerical FWHM = {fwhm / 1e12:.4f} THz, |dk L/2| <= 1.3916 rule = {closed / 1e12:.4f} THz")
print("  page (GV-matched DFG, 1560 nm): 23.7 THz")

header("6. Triply resonant ring OPO: threshold formula versus time-domain coupled-mode run")
# Page defaults: SHG process at 1550 nm (signal = idler = 1550 nm, pump 775 nm), n = 1.88/1.88/2.07,
# group indices 2.30/2.30/2.45 (v_g from their mean), L_rt = 440 um, kappa_p/2pi = 400 MHz,
# kappa_ex,p/2pi = 200 MHz, kappa_s/2pi = kappa_i/2pi = 200 MHz.
ls_, li_, lp_ = 1550e-9, 1550e-9, 775e-9
ws, wi, wp = [2 * np.pi * c / lam for lam in (ls_, li_, lp_)]
vg = c / np.mean([2.30, 2.30, 2.45])
Lr = 440e-6
kp, kexp, ks, ki = [2 * np.pi * f for f in (400e6, 200e6, 200e6, 200e6)]
etap = kexp / kp
e3 = eta_three(DEFF, AEFF, 1.88, 1.88, 2.07, lp_)
g0 = np.sqrt(e3 * hbar * ws * wi * vg ** 3 / (wp * Lr))
Pth = kp * ks * ki * wp ** 2 * Lr / (16 * etap * vg ** 3 * e3 * ws * wi)
print(f"  eta_SFG(775 nm) = {e3 * 1e-4 * 100:.1f} %/W/cm^2, g0/2pi = {g0 / 2 / np.pi / 1e6:.4f} MHz (page 2.562 MHz)")
print(f"  P_th = {Pth * 1e9:.2f} nW (page 490.5 nW)")


def opo_run(Pin, T=1.5e-6):
    """da_s/dt = -(ks/2) a_s - i g0 a_p a_i*, da_i/dt = -(ki/2) a_i - i g0 a_p a_s*,
    da_p/dt = -(kp/2) a_p - i g0 a_s a_i + sqrt(kex_p) s;  |a|^2 photons, |s|^2 photons/s."""
    s = np.sqrt(Pin / (hbar * wp))

    def f(t, y):
        ap = y[0] + 1j * y[1]
        a_s = y[2] + 1j * y[3]
        ai = y[4] + 1j * y[5]
        dp = -kp / 2 * ap - 1j * g0 * a_s * ai + np.sqrt(kexp) * s
        ds = -ks / 2 * a_s - 1j * g0 * ap * np.conj(ai)
        di = -ki / 2 * ai - 1j * g0 * ap * np.conj(a_s)
        return [dp.real, dp.imag, ds.real, ds.imag, di.real, di.imag]

    sol = solve_ivp(f, [0, T], [0, 0, 1e-3, 0, 1e-3, 0], rtol=1e-8, atol=1e-10, max_step=5e-10)
    return abs(sol.y[2, -1] + 1j * sol.y[3, -1]) ** 2


for fac in (0.9, 1.1, 1.5):
    print(f"  P = {fac:.1f} P_th: signal photons after 1.5 us (seed 1e-6) = {opo_run(fac * Pth):.3e}")
print("  below threshold the seed decays, above it grows: the formula marks the transition")

header("7. Ring spectrum with drop port and backscatter: linear solve, energy balance, limits")


def ring_linear(delta, k0, kex, kex2, g):
    """Solve da+/dt = z a+ + i g a- + sqrt(kex) s, da-/dt = z a- + i g a+ = 0 in steady state (s = 1)."""
    z = 1j * delta - (k0 + kex + kex2) / 2
    M = np.array([[z, 1j * g], [1j * g, z]])
    ap, am = np.linalg.solve(M, np.array([-np.sqrt(kex), 0.0]))
    T = abs(1 - np.sqrt(kex) * ap) ** 2
    return dict(T=T, D=kex2 * abs(ap) ** 2, Db=kex2 * abs(am) ** 2, R=kex * abs(am) ** 2,
                A=k0 * (abs(ap) ** 2 + abs(am) ** 2))


def ring_closed(delta, k0, kex, kex2, g):
    """The page's closed form: a+ = -sqrt(kex) z / (z^2 + g^2), a- = -i g a+ / z."""
    z = 1j * delta - (k0 + kex + kex2) / 2
    ap = -np.sqrt(kex) * z / (z * z + g * g)
    am = -1j * g * ap / z
    return dict(T=abs(1 - np.sqrt(kex) * ap) ** 2, D=kex2 * abs(ap) ** 2, Db=kex2 * abs(am) ** 2, R=kex * abs(am) ** 2)


err_bal, err_form = 0.0, 0.0
for k0, kex, kex2, g in ((0, 1, 1, 0.7), (0, 1, 0.3, 2.5), (0.4, 1, 0.3, 2.5), (0.2, 2, 0, 1.3), (1, 1, 0, 0)):
    for d in np.linspace(-5, 5, 41):
        a = ring_linear(d, k0, kex, kex2, g)
        b = ring_closed(d, k0, kex, kex2, g)
        err_bal = max(err_bal, abs(a["T"] + a["D"] + a["Db"] + a["R"] + a["A"] - 1))
        err_form = max(err_form, max(abs(a[k] - b[k]) for k in ("T", "D", "Db", "R")))
print(f"  max |T + D_fwd + D_bwd + R + absorbed - 1| = {err_bal:.2e}")
print(f"  max |linear solve - closed form| = {err_form:.2e}")
print(f"  T0 undercoupled (k0 = 3, kex = 1): {ring_linear(0, 3, 1, 0, 0)['T']:.6f} vs ((k0-kex)/k)^2 = {((3 - 1) / 4) ** 2:.6f}")
print(f"  T0 critical (k0 = kex): {ring_linear(0, 1, 1, 0, 0)['T']:.2e}")
print(f"  lossless symmetric add-drop D0: {ring_linear(0, 0, 1, 1, 0)['D']:.12f}")
# Page case: FSR 100 GHz, kappa0/2pi = kex/2pi = 20 MHz, kex2/2pi = 15 MHz, g/2pi = 40 MHz
tw = 2 * np.pi
r0 = ring_linear(0, tw * 20e6, tw * 20e6, tw * 15e6, tw * 40e6)
print(f"  page case (20/20/15 MHz, g/2pi = 40 MHz) at zero detuning: T = {100 * r0['T']:.3f} %, "
      f"D = {100 * r0['D']:.3f} % + {100 * r0['Db']:.3f} % bwd, R = {100 * r0['R']:.3f} %")
print("  page: T = 58.76 %, D = 4.086 % + 8.646 % bwd, R = 11.53 %")

header("8. Resonant EO comb: FFT of E(theta) = 1 / (1 - r e^{i beta sin theta}) versus the decay law")
for r, beta in ((0.97, 1.0), (0.99, 2.0), (0.9987442, 1.0)):
    F = np.pi * np.sqrt(r) / (1 - r)
    N = 1 << 18
    th = 2 * np.pi * np.arange(N) / N
    E = 1 / (1 - r * np.exp(1j * beta * np.sin(th)))
    Pn = np.abs(np.fft.fft(E) / N) ** 2
    n = np.arange(200, 600) if r > 0.995 else np.arange(20, 60)
    slope = np.polyfit(n, 10 * np.log10(Pn[n]), 1)[0]
    law = 20 * np.log10(np.e) * np.arcsinh(-np.log(r) / beta)
    print(f"  r = {r}, beta = {beta}: F = {F:.1f}, FFT slope = {-slope:.5f} dB/line, "
          f"20 log10(e) asinh(-ln r / beta) = {law:.5f} dB/line")
print("  page (F_L = 2500, beta = 1): 0.0109 dB/line")
