"""Top-hat wavelength bands shared by the still and video simulators."""

import math

MAX_SPECTRAL_SAMPLES = 128


def short_rel(bandwidth):
    """Shortest wavelength divided by the arithmetic band centre."""
    if not 0 <= bandwidth < 2:
        raise ValueError("bandwidth must be in [0, 2) for a positive lower wavelength")
    return 1 - bandwidth / 2


def from_um(lower, upper):
    """Convert GUI endpoints in micrometres to the existing nm/fraction pair."""
    try:
        lower, upper = float(lower), float(upper)
    except (TypeError, ValueError) as exc:
        raise ValueError("Enter numeric lower and upper wavelengths (µm)") from exc
    if not (math.isfinite(lower) and math.isfinite(upper) and 0 < lower <= upper):
        raise ValueError("Wavelengths must be finite, with 0 < lower <= upper (µm)")
    centre = lower + (upper - lower) / 2
    bandwidth = (upper - lower) / centre
    if not math.isfinite(centre * 1000) or bandwidth >= 2:
        raise ValueError("Wavelength range is too large to represent")
    return centre * 1000, bandwidth


def edges_um(lambda_nm, bandwidth):
    centre = lambda_nm * 1e-3
    return centre * short_rel(bandwidth), centre * (1 + bandwidth / 2)


def n_lambda(w020, bandwidth):
    short = short_rel(bandwidth)
    if bandwidth == 0:
        return 1
    # Midpoint bins: <=10% of the shortest wavelength and <= half a
    # wave of rim defocus change per bin, even at the short edge.
    bins = max(10 * bandwidth / short, 2 * abs(w020) * bandwidth / short**2)
    if not math.isfinite(bins) or bins > MAX_SPECTRAL_SAMPLES:
        raise ValueError(f"Band needs more than {MAX_SPECTRAL_SAMPLES} wavelength "
                         "samples; narrow the band or reduce defocus")
    return max(int(3 + 20 * bandwidth + 0.7 * bandwidth * abs(w020)),
               math.ceil(bins))
