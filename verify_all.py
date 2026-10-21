#!/usr/bin/env python3
"""Replay every check and compare its outputs byte for byte with checks/.

Usage: python3 verify_all.py [--figures]
Each script is rerun with --out pointing at a fresh temporary directory.
Data outputs (.json, .csv, .txt) must be byte-identical to the archived copies.
With --figures, PDF/PNG figures are regenerated and compared too (they are
byte-identical with the pinned library versions listed in README.md, but may
legitimately differ under other Matplotlib versions).
Exit status 0 means every compared file matched.
"""
import argparse, hashlib, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ('verify.py', ['verification.json', 'lattice_diagnostics.csv', 'saturation_counterexample.csv'],
     ['moment_and_saturation.pdf', 'moment_and_saturation.png', 'lattice_scaling.pdf', 'lattice_scaling.png']),
    ('double_exponential_prefactor.py', ['prefactor_decimal.txt'], []),
    ('manuscript_constant_check.py', ['manuscript_constant_check.json'], []),
    ('e4_ratio_one.py', ['e4_ratio_one.json'], []),
    ('e1_effective_length.py', ['e1_effective_length.json'], []),
    ('e2_spectrally_negative.py', ['e2_spectrally_negative.json'], []),
    ('fig_boundary_layers.py', [], ['boundary_layers.pdf']),
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--figures', action='store_true')
    figures = ap.parse_args().figures
    ok = True
    with tempfile.TemporaryDirectory(prefix='lq_replay_') as tmp:
        tmp = Path(tmp)
        for script, data, figs in SCRIPTS:
            if not data and not figures:
                print(f'skip  {script} (figure only; use --figures)')
                continue
            args = [sys.executable, str(ROOT / script), '--out', str(tmp)]
            if script == 'verify.py' and not figures:
                args.append('--no-plots')
            r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
            if r.returncode != 0:
                ok = False
                print(f'FAIL  {script} exited with {r.returncode}\n{r.stderr[-2000:]}')
                continue
            for name in data + (figs if figures else []):
                a, b = ROOT / 'checks' / name, tmp / name
                same = a.exists() and b.exists() and sha(a) == sha(b)
                ok &= same
                print(f"{'OK   ' if same else 'DIFF '} {name}  {sha(b)[:16] if b.exists() else 'missing'}")
    print('ALL OUTPUTS REPRODUCED BYTE-IDENTICALLY' if ok else 'MISMATCH: see lines above')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
