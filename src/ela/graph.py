from itertools import chain
from typing import Iterable, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb

from .slice import Slice

_default_theme = "whitegrid"
sb.set_theme(style=_default_theme)


def show(dem_data: np.ndarray, filename: str, title: Optional[str] = None):
    sb.set_theme(style="white")

    fig, axes = plt.subplots()
    ramp = "viridis"
    axes.imshow(dem_data, cmap=ramp)
    fig.colorbar(axes.pcolor(dem_data, cmap=ramp))
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename)

    sb.set_theme(style=_default_theme)


def prep_slices_for_histogram(slices: Iterable[Slice]) -> np.ndarray:
    # must be large so int() doesn't "quantize" the reconstruction
    NUM_SAMPLES = 1000000
    reconstruction = chain(*[[s.mid]*int(s.weight*NUM_SAMPLES) for s in slices])
    return np.fromiter(reconstruction, dtype=float)


def plot_elevation_histogram(data: Union[np.ndarray, Iterable[Slice]], ela: list[float], interval: float, filename: str,
                             title: str = "Vertical Distribution of Surface Area") -> None:
    """
    Create and save to disk a PNG histogram of the glacier surface information in `data`.
    A traditional and cumulative histogram will be produced and saved to `filename`. A vertical
    red line will be plotted at each elevation in `ela`.
    """
    values: np.ndarray
    if isinstance(data, np.ndarray):
        values = data
    else:
        try:
            values = prep_slices_for_histogram(data)
        except AttributeError:
            raise ValueError("parameter 'data' is not a numpy array or iterable of Slice")

    fig, axes = plt.subplots(nrows=2, sharex=True)
    fig.set_size_inches(w=8, h=8)

    # traditional histogram
    ax: plt.Axes = sb.histplot(values,
                               binwidth=interval, stat="percent",
                               ax=axes[0])
    ax.set_ylabel("% of Total")
    ax.set_xlabel("Elevation")
    for elev in ela:
        ax.axvline(elev, color='red')

    # cumulative histogram
    ax: plt.Axes = sb.histplot(values,
                               binwidth=interval, stat="percent", cumulative=True,
                               ax=axes[1])
    ax.set_ylabel("Cumulative %")
    ax.set_xlabel("Elevation")
    for elev in ela:
        ax.axvline(elev, color='red')

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename)
