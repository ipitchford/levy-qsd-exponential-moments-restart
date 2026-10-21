#!/usr/bin/env python3
"""Section 9 of the manuscript: effective-length (O(1/L)) corrections for the two-sided exponential model.

High-precision (mpmath) solver of the exact four-mode finite-interval QSD
problem for the centred (Esscher-tilted at eta) process, plus the closed-form
constants delta, g = -pi^2 chi'''(0)/(6 v) = (pi^2/2)(delta_- - delta_+),
c_4 = (5 chi'''(0)^2 - 3 v chi''''(0))/(72 v).

Usage:  e1_effective_length.py [--out DIR]   (default DIR: checks/ next to this file)
Output is deterministic: fixed significant digits, sorted keys, residuals as
threshold strings.
"""
import argparse
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 70
DIG = 20          # significant digits for finite-L values
CDIG = 30         # significant digits for constants
LS = (40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)


def S(x, n=DIG):
    return mp.nstr(mp.mpf(x), n, min_fixed=-5, max_fixed=8)


def thr(x, k):
    """Threshold string for a round-off level residual."""
    x = abs(mp.mpf(x))
    assert x < mp.mpf(10) ** (-k), (x, k)
    return '<1e-%d' % k


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pmul(a, b):
    o = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            o[i + j] += x * y
    return o


def psc(a, k):
    return [k * x for x in a]


class Model:
    """psi(s) = mu s + s^2/2 + lu s/(bu-s) - ld s/(bd+s).

    Upward Levy density lu*bu*e^{-bu y} (y>0), downward ld*bd*e^{bd y} (y<0).
    """

    def __init__(self, mu, lu, bu, ld, bd, eta0, lam0):
        mu, lu, bu, ld, bd = (mp.mpf(str(x)) for x in (mu, lu, bu, ld, bd))
        self.lu, self.ld, self.bd = lu, ld, bd
        self.psi = lambda z: mu * z + z ** 2 / 2 + lu * z / (bu - z) - ld * z / (bd + z)
        d1 = lambda z: mu + z + lu * bu / (bu - z) ** 2 - ld * bd / (bd + z) ** 2
        self.eta = e = mp.findroot(d1, mp.mpf(eta0))
        self.lam = mp.findroot(self.psi, mp.mpf(lam0))
        assert 0 < e < self.lam < bu
        self.kinf = -self.psi(e)
        self.v = 1 + 2 * lu * bu / (bu - e) ** 3 + 2 * ld * bd / (bd + e) ** 3
        self.k3 = 6 * lu * bu / (bu - e) ** 4 - 6 * ld * bd / (bd + e) ** 4
        self.k4 = 24 * lu * bu / (bu - e) ** 5 + 24 * ld * bd / (bd + e) ** 5
        self.bp, self.bm = bp, bm = bu - e, bd + e
        cu, cd = lu * bu, ld * bd
        mut = cd / bm ** 2 - cu / bp ** 2          # centred drift: chi'(0)=0
        # N(s) = 2(bp+s)(bm-s)(chi(-s)+alpha) = P0(s) + alpha*base2(s)
        base = pmul([bp, 1], [bm, -1])
        P = pmul(base, [0, -2 * mut, 1])
        P = padd(P, psc(pmul([0, -1], [bm, -1]), 2 * cu / bp))
        P = padd(P, psc(pmul([0, 1], [bp, 1]), 2 * cd / bm))
        self.P0, self.base2 = P, psc(base, 2)
        self.chk = {'chi_double_zero': max(abs(P[0]), abs(P[1]))}
        r = sorted(mp.re(x) for x in mp.polyroots([P[4], P[3], P[2]]))
        self.ap, self.am = -r[0], r[1]
        assert self.ap > bp and self.am > bm
        self.dp = 1 / bp - 1 / self.ap
        self.dm = 1 / bm - 1 / self.am
        self.delta = self.dp + self.dm
        self.g = -mp.pi ** 2 * self.k3 / (6 * self.v)
        self.c4 = (5 * self.k3 ** 2 - 3 * self.v * self.k4) / (72 * self.v)
        self.Gp = lambda z: bp * (self.ap - z) / (self.ap * (bp - z))
        self.Gm = lambda z: bm * (self.am - z) / (self.am * (bm - z))
        b = self.lam - e
        self.Cstar = self.Gp(b) / self.Gm(e)
        self.Rstar = self.kinf * self.Gp(-e) / self.Gm(e)
        # consistency of closed forms
        chi = lambda u: self.psi(e + u) - self.psi(e)
        u0 = mp.mpf('0.3')
        self.chk['factorisation'] = abs(chi(u0) - u0 ** 2 * (self.ap - u0) * (self.am + u0)
                                        / (2 * (bp - u0) * (bm + u0)))
        self.chk['v_equals_apam_over_bpbm'] = abs(self.v - self.ap * self.am / (bp * bm))
        self.chk['k3_equals_3v(dp-dm)'] = abs(self.k3 - 3 * self.v * (self.dp - self.dm))
        self.chk['k3_numeric_diff'] = abs(self.k3 - mp.diff(self.psi, e, 3))
        self.chk['Cstar_endpoint_form'] = abs(e ** 2 * self.Gp(-e) / (b ** 2 * self.Gm(-b)) - self.Cstar)
        self.chk['Rstar_endpoint_form'] = abs(e ** 2 * self.Gp(-e) * (self.am / (2 * bm)
                                              + lu / (bp ** 2 * self.Gm(-bp))) - self.Rstar)

    def roots(self, al):
        p = padd(self.P0, psc(self.base2, al))
        rr = mp.polyroots(p[::-1], maxsteps=200, extraprec=200)
        tiny = mp.mpf(10) ** (-40)
        real = sorted(mp.re(z) for z in rr if abs(mp.im(z)) < tiny)
        cx = [z for z in rr if mp.im(z) > tiny]
        assert len(real) == 2 and len(cx) == 1, (al, rr)
        return real[0], real[1], cx[0]

    def matrix(self, al, L):
        sn, sp, s = self.roots(al)
        bp, bm = self.bp, self.bm
        e1, e2, ec = mp.exp(sn * L), mp.exp(-sp * L), mp.exp(s * L)
        # columns: e^{sn y}, e^{sp(y-L)}, Im e^{sy}, Re e^{sy} (slow part Im(K e^{sy}))
        M = mp.matrix([[1, e2, 0, 1],
                       [1 / (sn + bp), e2 / (sp + bp), mp.im(1 / (s + bp)), mp.re(1 / (s + bp))],
                       [e1, 1, mp.im(ec), mp.re(ec)],
                       [e1 / (bm - sn), 1 / (bm - sp), mp.im(ec / (bm - s)), mp.re(ec / (bm - s))]])
        return sn, sp, s, M

    def solve(self, L):
        L = mp.mpf(L)
        a0 = mp.pi ** 2 * self.v / (2 * (L + self.delta) ** 2) + self.c4 * mp.pi ** 4 / L ** 4
        D = lambda a: mp.det(self.matrix(a, L)[3])
        al = mp.findroot(D, (a0, a0 * (1 + mp.mpf(10) ** -8)), solver='secant',
                         tol=mp.mpf(10) ** (-120))
        sn, sp, s, M = self.matrix(al, L)
        x = []
        for j in range(4):
            sub = mp.matrix(3, 3)
            cols = [k for k in range(4) if k != j]
            for i in range(3):
                for jj, k in enumerate(cols):
                    sub[i, jj] = M[i, k]
            x.append((-1) ** j * mp.det(sub))
        x = mp.matrix(x)
        bres = mp.norm(M * x, mp.inf) / mp.norm(x, mp.inf)
        c1, c2, K = x[0], x[1], mp.mpc(x[2], x[3])
        p = lambda y: c1 * mp.exp(sn * y) + c2 * mp.exp(sp * (y - L)) + mp.im(K * mp.exp(s * y))
        if p(L / 2) < 0:
            c1, c2, K = -c1, -c2, -K
        p = lambda y: c1 * mp.exp(sn * y) + c2 * mp.exp(sp * (y - L)) + mp.im(K * mp.exp(s * y))
        dp = lambda y: (c1 * sn * mp.exp(sn * y) + c2 * sp * mp.exp(sp * (y - L))
                        + mp.im(K * s * mp.exp(s * y)))

        def I(th):   # int_0^L e^{-th y} p(y) dy
            return (c1 * (mp.exp((sn - th) * L) - 1) / (sn - th)
                    + c2 * (mp.exp(-th * L) - mp.exp(-sp * L)) / (sp - th)
                    + mp.im(K * (mp.exp((s - th) * L) - 1) / (s - th)))

        def J(b):    # int_0^L e^{-b(L-y)} p(y) dy
            return (c1 * (mp.exp(sn * L) - mp.exp(-b * L)) / (sn + b)
                    + c2 * (1 - mp.exp(-(sp + b) * L)) / (sp + b)
                    + mp.im(K * (mp.exp(s * L) - mp.exp(-b * L)) / (s + b)))

        e, b = self.eta, self.lam - self.eta
        Z = I(e)
        C = J(b) / Z
        up = -dp(L) / 2 + self.lu * J(self.bp)        # e^{eta L} x upper exit flux
        R = up / Z
        low = dp(0) / 2 + self.ld * I(e + self.bd)    # lower exit flux
        flux = abs((low + mp.exp(-e * L) * up) / Z - (self.kinf + al)) / (self.kinf + al)
        # direct residuals of the four boundary conditions, relative to |p|_max scale
        sc = max(abs(c1), abs(c2), abs(K))
        r_jl = abs(c1 / (sn + self.bp) + c2 * mp.exp(-sp * L) / (sp + self.bp)
                   + mp.im(K / (s + self.bp))) / sc
        r_jr = abs(c1 * mp.exp(sn * L) / (self.bm - sn) + c2 / (self.bm - sp)
                   + mp.im(K * mp.exp(s * L) / (self.bm - s))) / sc
        w, xi = mp.im(s), mp.re(s)
        thl = mp.atan(w / (self.bp + xi)) - mp.atan(w / (-sn + xi))
        thr_ = mp.atan(w / (self.bm - xi)) - mp.atan(w / (sp - xi))
        phase = w * L + thl + thr_ - mp.pi
        vals = [p(L * k / 128) for k in range(1, 128)]
        assert min(vals) > 0 and Z > 0
        return dict(alpha=al, C=C, R=R, bres=bres, flux=flux, rjl=r_jl, rjr=r_jr,
                    pbc=max(abs(p(0)), abs(p(L))) / sc, phase=phase, xi=xi, w=w,
                    thl=thl, thr=thr_)

    def second_order(self):
        """omega^2 coefficients of the even end-functional ratios Q_C, Q_R (L = infinity)."""
        w = mp.mpf(10) ** -15
        a = mp.findroot(lambda al: mp.im(self.roots(al)[2]) - w, self.v * w ** 2 / 2)
        sn, sp, s = self.roots(a)
        xi = mp.re(s)
        thl = mp.atan(w / (self.bp + xi)) - mp.atan(w / (-sn + xi))
        thr_ = mp.atan(w / (self.bm - xi)) - mp.atan(w / (sp - xi))
        e, b = self.eta, self.lam - self.eta
        Ll = -mp.sin(thl) / (e - sn) + mp.im(mp.expj(thl) / (e - s))
        Lr = lambda bb: -mp.sin(thr_) / (bb + sp) + mp.im(mp.expj(thr_) / (bb + xi - 1j * w))
        Rf = (sp * mp.sin(thr_) - xi * mp.sin(thr_) + w * mp.cos(thr_)) / 2 + self.lu * Lr(self.bp)
        hC = (Lr(b) / Ll / self.Cstar - 1) / w ** 2
        hR = (Rf / Ll / self.Rstar - 1) / w ** 2
        return hC, hR


def richardson(Ls, f, orders):
    """Repeated Richardson on doubling sequence, f = f_inf + sum c_k L^{-k}."""
    vals = [f[L] for L in Ls]
    for k in orders:
        vals = [(2 ** k * vals[i + 1] - vals[i]) / (2 ** k - 1) for i in range(len(vals) - 1)]
    return vals[-1]


def run_model(m, Ls, full=True):
    out = {}
    hC, hR = m.second_order()
    g, d = m.g, m.delta
    sC = g ** 2 / 2 - 2 * d * g + mp.pi ** 2 * hC
    sR = g ** 2 / 2 - 2 * d * g + mp.pi ** 2 * hR
    out['constants'] = {
        'eta': S(m.eta, CDIG), 'lambda_star': S(m.lam, CDIG), 'kappa_inf': S(m.kinf, CDIG),
        'v': S(m.v, CDIG), 'chi3': S(m.k3, CDIG), 'chi4': S(m.k4, CDIG),
        'beta_plus': S(m.bp, CDIG), 'beta_minus': S(m.bm, CDIG),
        'a_plus': S(m.ap, CDIG), 'a_minus': S(m.am, CDIG),
        'delta_plus': S(m.dp, CDIG), 'delta_minus': S(m.dm, CDIG), 'delta': S(d, CDIG),
        'g_C_equals_g_R': S(g, CDIG), 'c4': S(m.c4, CDIG), 'c4_pi4': S(m.c4 * mp.pi ** 4, CDIG),
        'C_star': S(m.Cstar, CDIG), 'R_star': S(m.Rstar, CDIG),
        'pi2v_over_2': S(mp.pi ** 2 * m.v / 2, CDIG),
        'second_order_C_coefficient_computed': S(sC, 12),
        'second_order_R_coefficient_computed': S(sR, 12),
    }
    out['closed_form_checks'] = {k: thr(v, 60) for k, v in sorted(m.chk.items())}
    rows, fC, fR, fA = [], {}, {}, {}
    for L in Ls:
        r = m.solve(L)
        Lm = mp.mpf(L)
        a2 = Lm ** 2 * r['alpha']
        pred = mp.pi ** 2 * m.v / 2 * (Lm / (Lm + d)) ** 2
        eC, eR = r['C'] / m.Cstar - 1, r['R'] / m.Rstar - 1
        fC[L], fR[L], fA[L] = Lm * eC, Lm * eR, Lm ** 2 * (a2 - pred)
        row = {
            'L': L, 'L2_alpha_L': S(a2), 'C_L': S(r['C']), 'R_L': S(r['R']),
            'L2alpha_minus_effective_length_prediction': S(a2 - pred, 12),
            'L2_times_(L2alpha_minus_prediction)': S(fA[L], 12),
            'L_times_(C_L/C_star-1)': S(fC[L], 15),
            'L_times_(R_L/R_star-1)': S(fR[L], 15),
            'L2_times_(C_L/C_star-1-g/L)': S(Lm ** 2 * (eC - g / Lm), 12),
            'L2_times_(R_L/R_star-1-g/L)': S(Lm ** 2 * (eR - g / Lm), 12),
            'theta_left_over_omega': S(r['thl'] / r['w'], 15),
            'theta_right_over_omega': S(r['thr'] / r['w'], 15),
            'xi_L': S(r['xi'] * Lm, 15),
            'boundary_matrix_residual': thr(r['bres'], 60),
            'jump_kernel_residual_left': thr(r['rjl'], 60),
            'jump_kernel_residual_right': thr(r['rjr'], 60),
            'dirichlet_residual': thr(r['pbc'], 60),
            'total_exit_flux_identity_relative_error': thr(r['flux'], 60),
            'phase_equation_residual': thr(r['phase'], 50),
        }
        rows.append(row)
    out['finite_L'] = rows
    tail = list(Ls[-4:])
    out['richardson'] = {
        'from_L': tail,
        'L(C_L/C*-1)_limit': S(richardson(tail, fC, (1, 2, 3)), 12),
        'L(R_L/R*-1)_limit': S(richardson(tail, fR, (1, 2, 3)), 12),
        'L2(L2alpha-pred)_limit': S(richardson(tail, fA, (1, 2, 3)), 12),
        'C_L_limit': S(richardson(tail, {L: m.solve(L)['C'] for L in tail}, (1, 2, 3)), 15),
    }
    out['max_abs_deviation_L>=640'] = {
        'L(C_L/C*-1)-g': S(max(abs(fC[L] - g) for L in Ls if L >= 640), 6),
        'L(R_L/R*-1)-g': S(max(abs(fR[L] - g) for L in Ls if L >= 640), 6),
        'L2alpha-(pi^2v/2)(L/(L+delta))^2': S(max(abs(fA[L]) / mp.mpf(L) ** 2 for L in Ls if L >= 640), 6),
    }
    return out


def brownian(a='0.8'):
    """BM with drift -a: eta=a, delta=0, chi'''=0 => g=0. Exact closed forms."""
    a = mp.mpf(a)
    rows = []
    for L in LS:
        k = mp.pi / L
        C = mp.mpf(1)   # exact: e^{-aL} M_L(2a) = 1
        R = (a ** 2 + k ** 2) / (2 * (1 + mp.exp(-a * L)))
        rows.append({'L': L, 'C_L': S(C), 'L_times_(R_L/R_star-1)': S(L * (R / (a ** 2 / 2) - 1), 12),
                     'L2_alpha_L_over_(pi2v/2)': '1'})
    return {'drift': S(a, 4), 'delta': '0', 'g': '0', 'rows': rows,
            'note': 'C_L=1 exactly, L^2 alpha_L=pi^2/2 exactly, R_L/R_*-1=O(L^-2): consistent with g=0, delta=0'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(Path(__file__).resolve().parent / 'checks'))
    args = ap.parse_args()
    m1 = Model(-1, 2, 4, 3, 1, '0.94', '2')
    e = m1.eta
    # closed forms specific to the source model
    v = m1.v
    delta_cf = (5 - mp.sqrt((5 - 4 * e) ** 2 + 4 * v * (4 - e) * (1 + e)) / v) / ((4 - e) * (1 + e))
    g_cf = (mp.pi ** 2 / 2) * (3 - 2 * e - (5 - 4 * e) / v) / ((4 - e) * (1 + e))
    g_cf2 = -mp.pi ** 2 * (-18 / (1 + e) ** 4 + 48 / (4 - e) ** 4) / (6 * (1 + 6 / (1 + e) ** 3 + 16 / (4 - e) ** 3))
    eta_poly = (e - 1) * (1 + e) ** 2 * (4 - e) ** 2 - 3 * (4 - e) ** 2 + 8 * (1 + e) ** 2
    res = {
        'description': 'Effective-length check (manuscript Section 9); exact four-mode solver at mp.dps=70',
        'source_model': {'psi': '-s+s^2/2-3s/(1+s)+2s/(4-s)', **run_model(m1, LS)},
        'source_model_closed_forms': {
            'eta_quintic_residual': thr(eta_poly, 60),
            'delta_closed_form_minus_numeric': thr(delta_cf - m1.delta, 60),
            'g_closed_form_minus_numeric': thr(g_cf - m1.g, 60),
            'g_cumulant_form_minus_numeric': thr(g_cf2 - m1.g, 60),
            'old_float_script_L80': {'L2_alpha': '11.65138333', 'C_L': '0.85538487', 'R_L': '0.91294774'},
        },
        'exploratory_second_model': {
            'psi': '-2s+s^2/2+s/(2-s)-s/(5+s)',
            **run_model(Model(-2, 1, 2, 1, 5, '0.5', '1.5'), (80, 160, 320, 640, 1280, 2560, 5120))},
        'exploratory_brownian_drift': brownian(),
    }
    old = res['source_model']['finite_L'][1]
    assert abs(mp.mpf(old['L2_alpha_L']) - mp.mpf('11.65138333')) < 1e-7
    assert abs(mp.mpf(old['C_L']) - mp.mpf('0.85538487')) < 1e-7
    assert abs(mp.mpf(old['R_L']) - mp.mpf('0.91294774')) < 1e-7
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    f = out / 'e1_effective_length.json'
    f.write_text(json.dumps(res, indent=2, sort_keys=True) + '\n')
    print('wrote', f)


if __name__ == '__main__':
    main()
