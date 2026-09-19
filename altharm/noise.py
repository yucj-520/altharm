"""The symmetry-forbidden harmonics, used as a measured noise floor.

This is the organising idea of the package.  Fit every harmonic up to ``lmax``,
not only the ones symmetry allows.  A coefficient at a forbidden ``l`` cannot be
physical, so whatever amplitude it carries is the numerical noise of the
calculation -- basis-set discontinuity as the plane-wave count changes around
the circle, incomplete eigenvalue convergence, residual density error.

That turns the question "is this high harmonic real?" into arithmetic: an
allowed harmonic is believable when it stands above the forbidden floor, and the
ratio is the number to quote.  It is not a substitute for convergence tests, but
it is the one diagnostic that costs nothing extra, because the forbidden
coefficients come out of the same fit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class NoiseReport:
    allowed: dict          # l -> amplitude, in the units of the fitted values
    forbidden: dict        # l -> amplitude
    floor: float           # RMS over all forbidden harmonics
    floor_max: float       # largest single forbidden amplitude
    residual_rms: float    # what the basis could not fit at all

    def snr(self):
        """Allowed amplitude divided by the forbidden floor, per ``l``."""
        if self.floor <= 0.0:
            return {l: np.inf for l in self.allowed}
        return {l: a / self.floor for l, a in self.allowed.items()}

    def credible(self, threshold: float = 5.0):
        """Allowed harmonics standing at least ``threshold`` times above the floor."""
        return {l: a for l, a in self.allowed.items() if self.snr()[l] >= threshold}

    def summary(self, unit_name="eV", unit_scale=1.0):
        lines = [
            f"noise floor (RMS over forbidden l) : {self.floor / unit_scale:12.6g} {unit_name}",
            f"worst single forbidden harmonic    : {self.floor_max / unit_scale:12.6g} {unit_name}",
            f"unfitted residual (RMS)            : {self.residual_rms / unit_scale:12.6g} {unit_name}",
            "",
            f"{'l':>4}  {'status':>9}  {'amplitude':>14}  {'/floor':>9}",
        ]
        snr = self.snr()
        for l in sorted(set(self.allowed) | set(self.forbidden)):
            if l in self.allowed:
                lines.append(
                    f"{l:4d}  {'allowed':>9}  {self.allowed[l] / unit_scale:14.6g}"
                    f"  {snr[l]:9.1f}"
                )
            else:
                lines.append(
                    f"{l:4d}  {'forbidden':>9}  {self.forbidden[l] / unit_scale:14.6g}"
                    f"  {'--':>9}"
                )
        return "\n".join(lines)


def noise_report(fit, allowed_l_set, skip_l0: bool = True) -> NoiseReport:
    """Split a fit into symmetry-allowed signal and symmetry-forbidden noise.

    ``l = 0`` is excluded from the floor by default.  It is the angular average
    of the splitting, which for a clean altermagnet is zero but which also
    absorbs any rigid offset between the two spin channels -- a Fermi-level or
    band-alignment shift, not angular noise, and including it would inflate the
    floor for a reason that has nothing to do with the harmonics.
    """
    amps = fit.amplitude_by_l()
    allowed, forbidden = {}, {}
    for l, a in amps.items():
        if skip_l0 and l == 0:
            continue
        (allowed if l in allowed_l_set else forbidden)[l] = a

    if forbidden:
        vals = np.array(list(forbidden.values()), dtype=float)
        floor = float(np.sqrt(np.mean(vals ** 2)))
        floor_max = float(np.max(vals))
    else:
        floor, floor_max = 0.0, 0.0

    return NoiseReport(
        allowed=allowed,
        forbidden=forbidden,
        floor=floor,
        floor_max=floor_max,
        residual_rms=fit.residual_rms,
    )
