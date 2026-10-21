#!/usr/bin/env python3
"""Semantic negative controls: plausible wrong versions of the claims must be rejected.

Standard library only. Exit status 0 means every control behaved as required.
1. Swapping the ladder orientation at the two ends (manuscript remark "Which ladder sits at which end")
   gives C* = 1.2926..., which must lie outside the certified enclosure of C*.
2. The ratio-one reading of the conjecture (C* = 1) must lie outside the enclosure.
3. A one-character change in a saved check output must be detected by SHA256SUMS.
4. The corrected budget scale: along B = exp((eta + eps) L) the restart success
   probability tends to one, so the linear relation p ~ r B must fail there.
"""
import hashlib, json, math, shutil, tempfile
from decimal import Decimal as D, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def constants():
    with localcontext() as ctx:
        ctx.prec = 60
        psi = lambda s: -s + s * s / 2 - 3 * s / (s + 1) + 2 * s / (4 - s)
        dpsi = lambda s: -1 + s - 3 / (s + 1) ** 2 + 8 / (4 - s) ** 2
        lo, hi = D(0), D(2)
        for _ in range(220):
            mid = (lo + hi) / 2
            if dpsi(mid) > 0:
                hi = mid
            else:
                lo = mid
        eta = (lo + hi) / 2
        rho0 = -psi(eta)
        d, e, b = 4 - eta, 1 + eta, 2 - eta
        p = 5 - 2 * eta
        disc = (p * p + 32 * rho0 / (eta * eta)).sqrt()
        A, B = -eta - (-p - disc) / 2, eta + (-p + disc) / 2
        Gp = lambda s: d * (A - s) / (A * (d - s))
        Gm = lambda s: e * (B - s) / (B * (e - s))
        return eta, Gp(b) / Gm(eta), Gm(b) / Gp(eta)


def main():
    enc = json.loads((ROOT / 'checks' / 'e4_ratio_one.json').read_text())['C_star_outward_interval_40d']
    lo, hi = D(enc[0]), D(enc[1])
    eta, correct, swapped = constants()
    assert lo <= correct <= hi, 'control baseline: correct orientation must lie in the enclosure'
    assert not (lo <= swapped <= hi) and abs(swapped - D('1.2926')) < D('1e-4'), 'swapped orientation not rejected'
    print(f'OK  swapped ladder orientation gives {swapped:.10f}, outside the certified enclosure')
    assert not (lo <= 1 <= hi), 'ratio-one reading not rejected'
    print('OK  ratio-one reading C* = 1 lies outside the certified enclosure')

    sums = {}
    for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
        h, name = line.split(maxsplit=1)
        sums[name.lstrip('*').removeprefix('./')] = h
    target = 'checks/e4_ratio_one.json'
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / 'tampered.json'
        shutil.copy(ROOT / target, f)
        raw = bytearray(f.read_bytes())
        i = raw.index(b'0.85014864899')
        raw[i + 12] = ord('8') if raw[i + 12] != ord('8') else ord('7')
        f.write_bytes(bytes(raw))
        assert hashlib.sha256((ROOT / target).read_bytes()).hexdigest() == sums[target], 'baseline hash mismatch'
        assert hashlib.sha256(f.read_bytes()).hexdigest() != sums[target], 'tampering not detected'
    print('OK  a one-digit change in checks/e4_ratio_one.json is detected by SHA256SUMS')

    R, L, eps = 0.9061044256808939, 30.0, 0.3
    B = math.exp((float(eta) + eps) * L)
    r = R * math.exp(-float(eta) * L)
    p = 1 - math.exp(-r * B)
    assert p > 0.999 and r * B > 1000, 'beyond e^{eta L} the linear regime p ~ r B must fail'
    print(f'OK  at B = exp((eta+0.3)L), L = 30: rB = {r * B:.3g} but p = {p:.6f}; linear regime rejected')
    print('ALL NEGATIVE CONTROLS PASSED')


if __name__ == '__main__':
    main()
