import math
from typing import Callable

SearchValue = Callable[[float], float]
SearchTest = Callable[[float], bool]


def binary_search(target: float, low: float, high: float,
                  value: SearchValue, found: SearchTest) -> float:
    """
    Perform a binary search for a `target` ratio at elevations between `low`
    and `high`. Trial ratio values are provided by the `value` function and
    are tested against `target` using the `found` function.
    """
    candidate: float = math.nan
    iterations: int = 0
    end_loop: bool = False
    while not end_loop:
        iterations += 1
        mid = (high+low)/2
        candidate = value(mid)
        # print(f"{iterations: 3d}{candidate: 14.7f}{mid: 10.2f}")
        if found(candidate) or high-low < 1e-6:
            end_loop = True
        elif candidate > target:
            low = mid
        else:
            high = mid
    return mid


def linear_search(target: float, low: float, high: float,
                  value: SearchValue, found: SearchTest) -> float:
    """
    Perform a linear search for a `target` ratio at elevations between `low`
    and `high`. Trial ratio values are provided by the `value` function and
    are tested against `target` using the `found` function.
    """
    MAX_INTERVAL: int = 50
    iterations: int = 0
    elev = low
    while elev <= high:
        iterations += 1
        candidate = value(elev)
        diff = abs(target-candidate)
        interval = max(0.25, min(diff*MAX_INTERVAL, MAX_INTERVAL))  # needs work
        # print(f"{iterations:3d}{candidate: 14.7f}{elev:10.2f}{interval:10.4f}")
        if found(candidate):
            break
        elev += interval
    return elev
