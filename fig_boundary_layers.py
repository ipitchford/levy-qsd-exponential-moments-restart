#!/usr/bin/env python3
"""Figure 1: the centred QSD of the source model is a sine bulk with two
ladder-renewal boundary layers.  Exact finite-interval solution (four-mode
boundary-value problem, as in manuscript_constant_check.py), float64.
Label: numerical illustration (no claim rests on this figure).
Usage: python3 fig_boundary_layers.py [--out DIR]
"""
import argparse, json, os, tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from manuscript_constant_check import decimal_constants

ROOT = Path(__file__).resolve().parent


def centred_density(L, eta, kinf, v):
    bp, bm = 4 - eta, 1 + eta
    def modes(alpha):
        roots = np.roots([1, 5 - 4 * eta, 2 * alpha - v * bp * bm, 2 * alpha * (bp - bm), -2 * alpha * bp * bm])
        real = sorted([z.real for z in roots if abs(z.imag) < 1e-9])
        s = [z for z in roots if z.imag > 1e-9][0]
        sn, sp = real
        en, ep, ec = np.exp(sn * L), np.exp(-sp * L), np.exp(s * L)
        M = np.array([[1, ep, 1, 0], [en, 1, ec.real, ec.imag],
                      [1 / (sn + bp), ep / (sp + bp), (1 / (s + bp)).real, (1 / (s + bp)).imag],
                      [en / (bm - sn), 1 / (bm - sp), (ec / (bm - s)).real, (ec / (bm - s)).imag]])
        return sn, sp, s, M
    sc = np.pi ** 2 * v / (2 * L ** 2)
    grid = sc * np.geomspace(.12, 2.0, 180)
    dets = [np.linalg.det(modes(x)[3]) for x in grid]
    br = [(x, y) for x, y, fx, fy in zip(grid[:-1], grid[1:], dets[:-1], dets[1:]) if fx * fy < 0][0]
    alpha = brentq(lambda a: np.linalg.det(modes(a)[3]), *br, xtol=1e-16, rtol=1e-13)
    sn, sp, s, M = modes(alpha)
    coeff = np.linalg.svd(M)[2][-1]
    def f(y):
        z = np.exp(s * y)
        return coeff @ np.array([np.exp(sn * y), np.exp(sp * (y - L)), z.real, z.imag])
    if f(L / 2) < 0:
        coeff = -coeff
    y = np.linspace(0, L, 4001)
    p = np.array([f(t) for t in y])
    p /= np.trapezoid(p, y)
    return y, p, alpha


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, default=ROOT / 'checks')
    out = ap.parse_args().out; out.mkdir(parents=True, exist_ok=True)
    k = decimal_constants()
    eta, kinf, v = (float(k[x]) for x in ('eta', 'c', 'v'))
    a_p, a_m = float(k['A']), float(k['D'])
    bp, bm = 4 - eta, 1 + eta
    dp, dm = 1 / bp - 1 / a_p, 1 / bm - 1 / a_m
    os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='mplcfg_'))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(1, 3, figsize=(12.8, 4.0), layout='constrained')
    cols = {12: '#b97726', 48: '#26746d', 192: '#73365b'}
    u = np.linspace(0, 1, 401)
    for L, c in cols.items():
        y, p, _ = centred_density(L, eta, kinf, v)
        axs[0].plot(y / L, L * p, color=c, label=f'L = {L}')
        z = np.linspace(0, 1.5, 300)
        axs[1].plot(y[y <= 1.5], L ** 2 * p[y <= 1.5], color=c)
        axs[2].plot((L - y)[(L - y) <= 1.5], L ** 2 * p[(L - y) <= 1.5], color=c)
    axs[0].plot(u, np.pi / 2 * np.sin(np.pi * u), 'k--', lw=1, label=r'$(\pi/2)\sin(\pi u)$')
    axs[0].set(xlabel=r'scaled position $u=y/L$', ylabel=r'$L\times$ centred QSD density', title='(a) bulk: sine profile')
    axs[0].legend(frameon=False, fontsize=8)
    z = np.linspace(0, 1.5, 300)
    for ax, d, a, name, lab in [(axs[1], dp, a_p, 'left', r'y'), (axs[2], dm, a_m, 'right', r'L-y')]:
        ax.plot(z, np.pi ** 2 / 2 * (z + d * (1 - np.exp(-a * z))), 'k--', lw=1, label=r'limit $\frac{\pi^2}{2}m V(\cdot)$')
        ax.plot(z, np.pi ** 2 / 2 * (z + d), ':', color='grey', lw=1, label='linear asymptote')
        ax.set(xlabel=f'distance to {name} end, ${lab}$', ylabel=r'$L^2\times$ density', ylim=(0, None))
        ax.legend(frameon=False, fontsize=8)
    axs[1].set_title(r'(b) left layer: ascending ladder $V_+$')
    axs[2].set_title(r'(c) right layer: descending ladder $V_-$')
    fig.savefig(out / 'boundary_layers.pdf', metadata={'CreationDate': None, 'Producer': None, 'Creator': None})
    plt.close(fig)
    print('wrote', out / 'boundary_layers.pdf')


if __name__ == '__main__':
    main()
