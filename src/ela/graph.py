from itertools import chain
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
from slice import Slice


def prep_slices_for_histogram(slices: Iterable[Slice]) -> np.ndarray:
    # must be large so int() doesn't "quantize" the reconstruction
    NUM_SAMPLES = 1000000
    reconstruction = chain(*[[s.mid]*int(s.weight*NUM_SAMPLES) for s in slices])
    return np.fromiter(reconstruction, dtype=float)


def plot_elevation_histogram(values: np.ndarray, ela: float, interval: float, filename: str,
                             title: str = "Vertical Distribution of Surface Area") -> None:
    fig, axes = plt.subplots(nrows=2, sharex=True)
    fig.set_size_inches(w=8, h=8)
    # traditional histogram
    ax: plt.Axes = sb.histplot(values,
                               binwidth=interval, stat="percent",
                               ax=axes[0])
    ax.set_ylabel("% of Total")
    ax.set_xlabel("Elevation")
    ax.axvline(ela, color='red')
    # cumulative histogram
    ax: plt.Axes = sb.histplot(values,
                               binwidth=interval, stat="percent", cumulative=True,
                               ax=axes[1])
    ax.set_ylabel("Cumulative %")
    ax.set_xlabel("Elevation")
    ax.axvline(ela, color='red')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename)
