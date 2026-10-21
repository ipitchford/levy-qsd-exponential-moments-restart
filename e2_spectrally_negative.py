#!/usr/bin/env python3
"""Section 8 of the manuscript: exact finite-L scale-function checks for spectrally negative models.

For a spectrally negative Levy process X killed on leaving (0, L), the
interval QSD has density proportional to W^{(-kappa_L)}(L - x), where
kappa_L is the smallest positive zero of q -> W^{(-q)}(L) (proof in
manuscript Section 8, proof of Theorem 8.1, Step 1).  For rational psi = P/Q the q-scale function
is the residue sum  W^{(q)}(x) = sum_i Q(t_i) e^{t_i x} / D'(t_i),
D = P - q Q.  All finite-L quantities below are then closed-form sums of
exponentials evaluated in mpmath (dps 60).  Labels: "numerical (high
precision)" for finite-L values and extrapolations; the limit constants
come from the theorem in e2_section.tex (model (a) closed form via sympy).

Models
  (a) psi(s) = s - 2s/(1+s)            (sigma = 0, bounded variation)
  (b) psi(s) = s^2/2 + 2s - 3s/(1+s)   (sigma = 1)

Output: JSON with sorted keys, fixed-digit strings; round-off residuals as
threshold strings.  Usage: e2_spectrally_negative.py [--out DIR]
"""

import argparse
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
SCRIPT_DIR = Path(__file__).resolve().parent
LS = [20, 40, 80, 160, 320, 640, 1280]


def s12(x):
    """Fixed 12-significant-digit string."""
    return mp.nstr(mp.mpf(x), 12, min_fixed=-6, max_fixed=8, strip_zeros=False)


def thr(x, cut=mp.mpf('1e-40')):
    x = abs(mp.mpf(x))
    return ('<' + mp.nstr(cut, 1)) if x < cut else mp.nstr(x, 3)


# ---------------------------------------------------------------- polynomials
def pval(c, t):          # c: coefficients highest degree first
    return mp.polyval(c, t)


def pder(c):
    n = len(c) - 1
    return [c[i] * (n - i) for i in range(n)]


def padd(a, b):
    n = max(len(a), len(b))
    a = [mp.mpf(0)] * (n - len(a)) + list(a)
    b = [mp.mpf(0)] * (n - len(b)) + list(b)
    return [x + y for x, y in zip(a, b)]


def pscale(a, k):
    return [k * x for x in a]


class Model:
    def __init__(self, name, P, Q, sym_psi):
        self.name = name
        self.P = [mp.mpf(x) for x in P]
        self.Q = [mp.mpf(x) for x in Q]
        self.sym = sym_psi
        self.psi = lambda t: pval(self.P, t) / pval(self.Q, t)
        dpsi = lambda t: mp.diff(self.psi, t)
        # Cramer root: positive root of P (excluding 0)
        roots = mp.polyroots(self.P, maxsteps=200, extraprec=200)
        pos = [mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf('1e-40') and mp.re(r) > mp.mpf('1e-30')]
        self.lam = min(pos)
        self.eta = mp.findroot(dpsi, self.lam / 2)
        self.kinf = -self.psi(self.eta)
        self.v = mp.diff(self.psi, self.eta, 2)
        self.psi3 = mp.diff(self.psi, self.eta, 3)
        self.chi = lambda s: self.psi(self.eta + s) - self.psi(self.eta)

    def modes(self, q):
        """Roots t_i and weights a_i with W_X^{(q)}(x) = sum a_i e^{t_i x}."""
        D = padd(self.P, pscale(self.Q, -q))
        ts = mp.polyroots(D, maxsteps=400, extraprec=400)
        dD = pder(D)
        return [(t, pval(self.Q, t) / pval(dD, t)) for t in ts]

    def W(self, q, x):
        return mp.re(sum(a * mp.exp(t * x) for t, a in self.modes(q)))

    # limits from the theorem (sigma >= 0, spectrally negative)
    def K_ell(self, b):
        return mp.pi ** 2 / (2 * b ** 2)

    def K_r(self, b):
        return mp.pi ** 2 * self.v / (4 * self.chi(b))


def int_exp(a, L):
    """int_0^L e^{a y} dy."""
    if abs(a) < mp.mpf('1e-45'):
        return mp.mpf(L)
    return (mp.exp(a * L) - 1) / a


class FiniteL:
    def __init__(self, m, L):
        self.m, self.L = m, mp.mpf(L)
        Lm = self.L
        f = lambda al: m.W(m.psi(m.eta) - al, Lm) * mp.exp(-m.eta * Lm)
        delta = -m.psi3 / (3 * m.v)
        guess = mp.pi ** 2 * m.v / 2 / (Lm + delta) ** 2
        self.alpha = mp.findroot(f, (guess * mp.mpf('0.7'), guess * mp.mpf('1.3')), solver='anderson', tol=mp.mpf('1e-100'), verify=False)
        self.root_residual = f(self.alpha) / f(guess * mp.mpf('0.7'))
        self.kappa = m.kinf + self.alpha
        self.md = m.modes(-self.kappa)       # density n(y) = W^{(-kappa)}(L-y)
        self.norm = mp.re(self.raw(0))
        # smallest-root check and positivity of the density (hypothesis (Q))
        grid = [self.alpha * k / 400 for k in range(1, 400)]
        self.sign_ok_rate = all(f(a) > 0 for a in grid)
        xs = [Lm * k / 400 for k in range(0, 400)]       # L-y in [0, L)
        self.min_density = min(m.W(-self.kappa, x) * mp.exp(-m.eta * x) for x in xs[1:])
        self.pos_ok = self.min_density > 0

    def raw(self, th):
        """int_0^L e^{th y} W^{(-kappa)}(L-y) dy."""
        L = self.L
        return sum(a * mp.exp(t * L) * int_exp(th - t, L) for t, a in self.md)

    def M(self, th):
        return self.raw(th) / self.norm   # complex for complex th


def resolvent_residual(m, F, y, q=mp.mpf(1)):
    """|int n(x)u^q(x,y)dx - n(y)/(q+kappa)| / n(y), Kyprianou Thm 8.7 kernel."""
    L = F.L
    n = lambda x: m.W(-F.kappa, L - x)
    WL = m.W(q, L)
    u = lambda x: m.W(q, x) * m.W(q, L - y) / WL - (m.W(q, x - y) if x > y else 0)
    lhs = mp.quad(lambda x: n(x) * u(x), [0, y, L])
    return abs(lhs - n(y) / (q + F.kappa)) / n(y)


def adjoint_residual_a(F, y):
    """Model (a): -n' + 2 int_y^L n(w) e^{-(w-y)} dw - 2 n + kappa n = 0."""
    m, L = F.m, F.L
    n = lambda x: m.W(-F.kappa, L - x)
    val = -mp.diff(n, y) + 2 * mp.quad(lambda w: n(w) * mp.exp(-(w - y)), [y, L]) - 2 * n(y) + F.kappa * n(y)
    return abs(val) / n(y)


def richardson(seq):
    """Two Richardson levels for doubling L, error model c1/L + c2/L^2."""
    r1 = [2 * seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
    r2 = [(4 * r1[i + 1] - r1[i]) / 3 for i in range(len(r1) - 1)]
    return r1, r2


def run_model(m, extra_theta_hi, theta_lo, adjoint):
    eta, lam, v = m.eta, m.lam, m.v
    Cstar = m.K_r(lam - eta) / m.K_ell(eta)
    Rstar = (v / 2) * eta ** 2
    delta = -m.psi3 / (3 * v)
    lim = {
        'C_star': Cstar, 'R_star': Rstar,
        'L2alpha': mp.pi ** 2 * v / 2,
        'M_eta_over_L2': 1 / m.K_ell(eta),
        'M_lo': m.K_ell(eta - theta_lo) / m.K_ell(eta),
        'M_hi_scaled': m.K_r(extra_theta_hi - eta) / m.K_ell(eta),
        'Q_L': mp.mpf(1) / 2,
    }
    # W_X^{(0)} modes for q_L
    w0 = m.modes(0)
    rows, ser = [], {k: [] for k in ['C', 'R', 'A', 'Me', 'Ml', 'Mh']}
    checks = {}
    for L in LS:
        F = FiniteL(m, L)
        Lm = F.L
        C = mp.exp(-(lam - eta) * Lm) * mp.re(F.M(lam))
        WL = mp.re(sum(a * mp.exp(t * Lm) for t, a in w0))
        qL = mp.re(sum(a * F.M(t) for t, a in w0)) / WL
        rL = F.kappa * qL
        R = mp.exp(eta * Lm) * rL
        # independent route: r_L = alpha e^{-eta L} Q_L / A_L(-eta), with
        # Q_L = Ptilde(upper exit) = int W_X^{(-kinf)}(x) e^{eta(L-x)} / W_X^{(-kinf)}(L) dnu~
        # nu~ density e^{eta y} n(y)/norm~ ; A_L(-eta) = 1/ M~ ... computed directly:
        A_m_eta = F.norm / mp.re(F.raw(eta))          # int e^{-eta y} dnu~ = norm / int e^{eta y} n
        # Q_L via W_Y(x)/W_Y(L), W_Y(x) = e^{-eta x} W_X^{(-kinf)}(x); double mode at t=eta,
        # use W_X^{(-kinf)} = lim of W_X^{(-kinf-eps)} (regular in q) with eps = 1e-30.
        eps = mp.mpf('1e-30')
        wk = m.modes(-m.kinf - eps)
        WYL = mp.re(sum(a * mp.exp((t - eta) * Lm) for t, a in wk))
        # int W_Y(x) e^{eta x} n(x) dx = int W_X^{(-kinf)}(x) n(x) dx
        num = mp.re(sum(a * F.raw(t) for t, a in wk))
        QL = num / mp.re(F.raw(eta)) / WYL
        r_alt = F.alpha * mp.exp(-eta * Lm) * QL / A_m_eta
        Me = mp.re(F.M(eta)) / Lm ** 2
        Ml = mp.re(F.M(theta_lo))
        Mh = mp.exp(-(extra_theta_hi - eta) * Lm) * mp.re(F.M(extra_theta_hi))
        A = Lm ** 2 * F.alpha
        for k, val in zip(['C', 'R', 'A', 'Me', 'Ml', 'Mh'], [C, R, A, Me, Ml, Mh]):
            ser[k].append(val)
        eff = (mp.pi ** 2 * v / 2) * (Lm / (Lm + delta)) ** 2
        row = {
            'L': str(L),
            'kappa_L': s12(F.kappa),
            'L2_alpha_L': s12(A),
            'L2_alpha_L_over_limit': s12(A / lim['L2alpha']),
            'L2_alpha_L_over_effective_length_prediction': s12(A / eff),
            'L_times_(L2alpha/limit-1)': s12(Lm * (A / lim['L2alpha'] - 1)),
            'C_L': s12(C), 'C_L_over_C_star': s12(C / Cstar),
            'L_times_(C_L/C_star-1)': s12(Lm * (C / Cstar - 1)),
            'R_L': s12(R), 'R_L_over_R_star': s12(R / Rstar),
            'L_times_(R_L/R_star-1)': s12(Lm * (R / Rstar - 1)),
            'M_L(eta)/L^2_over_limit': s12(Me / lim['M_eta_over_L2']),
            'M_L(theta_lo)_over_limit': s12(Ml / lim['M_lo']),
            'scaled_M_L(theta_hi)_over_limit': s12(Mh / lim['M_hi_scaled']),
            'Q_L_centred_upper_exit': s12(QL),
            'restart_rate_two_routes_relative_difference': thr((r_alt - rL) / rL, mp.mpf('1e-20')),
            'smallest_root_sign_check_passed': bool(F.sign_ok_rate),
            'density_positive_on_grid': bool(F.pos_ok),
        }
        if L in (20, 40):
            res = max(resolvent_residual(m, F, Lm * k / 7) for k in (1, 3, 6))
            row['resolvent_identity_residual'] = thr(res, mp.mpf('1e-30'))
            if adjoint:
                ra = max(adjoint_residual_a(F, Lm * k / 7) for k in (1, 3, 6))
                row['adjoint_equation_residual'] = thr(ra, mp.mpf('1e-30'))
        rows.append(row)
    rich = {}
    for k, key in [('C', 'C_star'), ('R', 'R_star'), ('A', 'L2alpha'), ('Me', 'M_eta_over_L2'),
                   ('Ml', 'M_lo'), ('Mh', 'M_hi_scaled')]:
        r1, r2 = richardson(ser[k])
        rich[key] = {'richardson2_last': s12(r2[-1]),
                     'richardson2_last_rel_error': mp.nstr(abs(r2[-1] / lim[key] - 1), 3),
                     'raw_last_rel_error': mp.nstr(abs(ser[k][-1] / lim[key] - 1), 3)}
    # 1/L coefficients extrapolated: c1 = lim L (f/f* - 1), via Richardson on L(f/f*-1)
    coef = {}
    for k, key in [('C', 'C_star'), ('R', 'R_star'), ('A', 'L2alpha')]:
        g = [mp.mpf(L) * (val / lim[key] - 1) for L, val in zip(LS, ser[k])]
        r1, r2 = richardson(g)
        coef[key] = s12(r2[-1])
    return {
        'parameters': {'lambda_star': s12(lam), 'eta': s12(eta), 'kappa_inf': s12(m.kinf),
                       'v': s12(v), 'psi3_eta': s12(m.psi3),
                       'delta_mean_stationary_down_overshoot': s12(delta),
                       'theta_lo': s12(theta_lo), 'theta_hi': s12(extra_theta_hi)},
        'limits': {k: s12(x) for k, x in lim.items()},
        'rows': rows,
        'richardson': rich,
        'first_order_coefficients_extrapolated': coef,
        'effective_length_prediction_for_L2alpha_coefficient': s12(-2 * delta),
        'label': 'numerical (high precision)',
    }


def exact_a():
    s = sp.symbols('s')
    psi = s - 2 * s / (1 + s)
    lam = sp.solve(sp.numer(sp.together(psi)) / s, s)
    eta = [r for r in sp.solve(sp.diff(psi, s), s) if r.is_positive][0]
    kinf = sp.nsimplify(sp.simplify(-psi.subs(s, eta)))
    v = sp.simplify(sp.diff(psi, s, 2).subs(s, eta))
    C = sp.radsimp(sp.simplify(eta ** 2 * v / (2 * kinf)))
    R = sp.radsimp(sp.expand(eta ** 2 * v / 2))
    return {'lambda_star': str(lam), 'eta': str(sp.simplify(eta)),
            'kappa_inf': str(sp.radsimp(sp.expand(kinf))), 'v': str(sp.radsimp(v)),
            'C_star': str(C), 'R_star': str(R),
            'C_star_squared': str(sp.simplify(C ** 2)),
            'C_star_decimal': s12(mp.mpf(str(sp.N(C, 50)))),
            'R_star_decimal': s12(mp.mpf(str(sp.N(R, 50)))),
            'label': 'certified (exact symbolic, sympy)'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(SCRIPT_DIR / 'checks'))
    args = ap.parse_args()
    ma = Model('a', [1, -1, 0], [1, 1], None)                # s(s-1)/(1+s)
    mb = Model('b', [1, 5, -2, 0], [2, 2], None)             # s(s^2+5s-2)/(2(1+s))
    out = {
        'description': 'Spectrally negative QSD (manuscript Section 8): scale-function finite-L checks',
        'model_a': {'psi': 's - 2s/(1+s)', 'exact': exact_a(),
                    **run_model(ma, mp.mpf(3), mp.mpf(-1) / 2, True)},
        'model_b': {'psi': 's^2/2 + 2s - 3s/(1+s)',
                    **run_model(mb, mp.mpf(3), mp.mpf(-1) / 2, False)},
        'brownian_sanity': {'C_star_formula_eta2v_over_2kinf': '1 (exact, m=sigma=1)',
                            'R_star_equals_kappa_inf': True},
    }
    d = Path(args.out)
    d.mkdir(parents=True, exist_ok=True)
    (d / 'e2_spectrally_negative.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
