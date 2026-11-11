"""``python -m altharm <data-file-schema.xml>`` -- harmonics from a finished run.

The whole point of the entry point is that the only thing it needs told is
where the run is.  Everything else -- the cell, the magnetic sublattices, the
circle's plane, the usable ``symprec``, which ``l`` symmetry allows -- is read
or measured from the run itself, because each of those is something that has
silently produced a wrong answer in this project when it was supplied by hand.
"""

from __future__ import annotations

import argparse
import sys

from .analyze import analyze_qe
from .symmetry import SymmetryToleranceError, format_symprec_scan


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="altharm",
        description="Angular harmonics of an altermagnet's spin splitting, "
                    "straight from a Quantum ESPRESSO bands run.",
    )
    p.add_argument("xml", help="path to <prefix>.save/data-file-schema.xml")
    p.add_argument("--band", default="fermi",
                   help="0-based band index, or 'fermi' (default) for the band "
                        "nearest E_F. Check the reported E - E_F: a band index "
                        "is not portable between runs.")
    p.add_argument("--lmax", type=int, default=None,
                   help="highest harmonic fitted (default: n_kpoints // 4)")
    p.add_argument("--symprec", type=float, default=None,
                   help="spglib tolerance (default: scan and pick the tightest "
                        "one that still finds a sublattice-exchanging operation)")
    p.add_argument("--soc", action="store_true",
                   help="the run included spin-orbit coupling, so E(k) = E(-k) "
                        "per spin channel does not hold and odd l are not forbidden")
    p.add_argument("--scan", action="store_true",
                   help="also print the full symprec scan")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    band = args.band if args.band == "fermi" else int(args.band)
    try:
        a = analyze_qe(
            args.xml,
            band=band,
            lmax=args.lmax,
            symprec=args.symprec,
            include_time_reversal=not args.soc,
        )
    except SymmetryToleranceError as exc:
        print(f"altharm: {exc}", file=sys.stderr)
        return 2
    if args.scan:
        print(format_symprec_scan(a.symprec_rows))
        print()
    print(a.summary())
    return 0


if __name__ == "__main__":
    sys.exit(main())
