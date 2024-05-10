import math
from typing import Iterable

import numpy as np

from .geoprocessing import Geoprocessor
from .search import SearchTest, SearchValue
from .slice import Slice


def absolute_difference(target: float, tolerance: float = 0.001) -> SearchTest:
    """
    Return a function that will test if a trial number is close to `target`
    within a a certain absolute difference.
    """
    return lambda trial: abs(trial-target) <= tolerance


def percent_difference(target: float, tolerance: float = 1.0) -> SearchTest:
    """
    Return a function that will test if a trial number is close to `target`
    within a certain percent difference.
    """
    return lambda trial: abs((trial-target)/target) <= (tolerance/100.0)


def br_from_slices(slices: Iterable[Slice]) -> SearchValue:
    """
    Return a function that calculates a AABR for a trial ELA using 
    pregenerated slices.
    """
    def ratio(trial_ela: float) -> float:
        acc = 0.0
        abl = 0.0
        for s in slices:
            n = (s.mid-trial_ela)*s.weight
            if n >= 0:
                acc += n
            else:
                abl += n

        return acc/abs(abl) if abl != 0 else math.inf

    return ratio


def aar_from_slices(slices: Iterable[Slice]) -> SearchValue:
    """
    Get a function that calculates an AAR for a trial ELA using
    pre-generated slices.
    """
    def ratio(trial_ela: float) -> float:
        return sum((s.weight for s in slices if s.mid >= trial_ela))
    return ratio

# get AAR by directly using surface area calcs at that elevation


def aar_direct(surface: str, be: Geoprocessor) -> SearchValue:
    """
    Get a function that calculates an AAR for a trial ELA by directly
    querying the surface area above the trial ELA.
    """
    _, max_elev = be.min_max_elevations(surface)
    total_area = be.surface_area(surface, max_elev)

    def ratio(trial_ela: float) -> float:
        # this surface_area() was originally the only one using "ABOVE"
        accumulation_area = total_area - be.surface_area(surface, trial_ela)
        return accumulation_area/total_area

    return ratio

# uses raster histogram to calculate ratios based on 2d surface area


def aar_from_elev_hist(elevations: np.ndarray, counts: np.ndarray) -> SearchValue:
    """
    Get a function that calculates an AAR based on a trial ELA using a 2D
    elevation histogram.
    """
    # remove last item from elevation bins since length of
    # bins and counts must match here
    if len(elevations) > len(counts):
        elevations = elevations[:-1]

    total_area = counts.sum()

    def ratio(trial_ela: float) -> float:
        upper_area = counts[elevations >= trial_ela].sum()
        return upper_area/total_area

    return ratio


def br_from_elev_hist(elevations: np.ndarray, counts: np.ndarray) -> SearchValue:
    """
    Get a function that calculates BR based on a trial ELA using a 2D
    elevation histogram.
    """
    # remove last item from elevation bins since length of
    # bins and counts must match here
    if len(elevations) > len(counts):
        elevations = elevations[:-1]

    def ratio(trial_ela: float) -> float:
        # don't need to divide the counts by the total area because the
        # total area figure just ends up canceling out in the numerator
        # and denominator
        upper_area = (counts[elevations >= trial_ela] *
                      (elevations[elevations >= trial_ela]-trial_ela)).sum()
        lower_area = (counts[elevations < trial_ela] *
                      (trial_ela-elevations[elevations < trial_ela])).sum()
        return upper_area/lower_area

    return ratio
