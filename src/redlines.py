"""
redlines.py - draw failure-limit ("redline") warnings on a matplotlib axis.

A redline is a limit that, once crossed, means the engine is at risk of failing (wall past its
melting point, coolant past its boiling limit, stress past yield, ...). draw_redline() draws the limit
in red, marks every station that crosses it, shades the offending span, and writes a banner on the plot.
"""

import numpy as np

RED = "#d62728"


def exceeded_mask(y, limit, above=True):
    """Boolean array: where y is past the limit (above it if above=True, below it otherwise)."""
    y = np.asarray(y, dtype=float)
    limit = np.broadcast_to(np.asarray(limit, dtype=float), y.shape)
    return (y > limit) if above else (y < limit)


def span_text(x, mask, unit="in", name="z"):
    """'z = a to b in' for the first..last offending station (name/unit configurable for non-axial x axes)."""
    xs = np.asarray(x, dtype=float)[mask]
    lo, hi = float(xs.min()), float(xs.max())
    u = f" {unit}" if unit else ""
    return f"{name} = {lo:.2f}{u}" if np.isclose(lo, hi) else f"{name} = {lo:.2f} to {hi:.2f}{u}"


def draw_redline(ax, x, y, limit, label, above=True, draw_limit=True, mask=None, banner_slot=0, x_name="z", x_unit="in"):
    """
    Draw one redline on `ax` and return a violation message (str) or None if the limit is respected.

    Parameters
    ----------
    x, y : arrays          plotted quantity along the axis
    limit : float | array  the redline value (ignored for the check when `mask` is given)
    label : str            what the redline is, e.g. "Hot wall melting point"
    above : bool           True if crossing means y > limit, False if y < limit
    draw_limit : bool      draw the limit line itself (False if the axis already shows it)
    mask : bool array      override the crossing test (e.g. a ratio computed elsewhere)
    banner_slot : int      stack index for the warning banner so several do not overlap
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if mask is None:
        mask = exceeded_mask(y, limit, above)
    mask = np.asarray(mask, dtype=bool)

    if draw_limit and limit is not None:
        lim_arr = np.broadcast_to(np.asarray(limit, dtype=float), x.shape)
        ax.plot(x, lim_arr, color=RED, linestyle="--", linewidth=1.6, label=f"REDLINE: {label}", zorder=4)

    if not mask.any():
        return None

    # Shade the offending span and mark each crossing station.
    ax.fill_between(x, 0, 1, where=mask, transform=ax.get_xaxis_transform(), color=RED, alpha=0.12, zorder=1)
    ax.scatter(x[mask], y[mask], marker="x", color=RED, s=36, zorder=5, label=f"{label} EXCEEDED")

    msg = f"{label} exceeded at {span_text(x, mask, x_unit, x_name)}"
    ax.text(
        0.99, 0.97 - 0.09 * banner_slot, f"⚠ REDLINE: {msg}",
        transform=ax.transAxes, va="top", ha="right", fontsize=9, fontweight="bold", color="white",
        bbox=dict(boxstyle="round,pad=0.25", facecolor=RED, edgecolor="none", alpha=0.9), zorder=6,
    )
    return msg
