import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import repeat
from typing import Iterable

from ela import backend


@dataclass
class Slice:
    top: float
    mid: float
    bottom: float
    area: float = 0
    weight: float = 0


def _determine_worker_count() -> int:
    cpus = os.cpu_count()
    if cpus:
        workers: int = round(cpus * 2/3)  # may not report correctly in docker containers
    else:
        workers = 4
    return workers


def create_slices(surface: str, interval: float) -> Iterable[Slice]:
    min_elev, max_elev = backend.min_max_elevations(surface)
    total_area = backend.surface_area(surface, max_elev)
    elev = min_elev+interval
    prev_area = 0
    slices: list[Slice] = []
    while elev <= max_elev + interval:
        area = backend.surface_area(surface, elev)
        slice_area = area-prev_area
        slices.append(Slice(
            top=elev,
            mid=elev-(interval/2),
            bottom=elev-interval,
            area=slice_area,
            weight=slice_area/total_area))
        elev += interval
        prev_area = area

    return slices


def slice_surface_area(data: tuple[Slice, str, float]) -> Slice:
    """
    Calculate the surface area of the Slice, which has values for
    top, mid, and bottom, at a minimum. Returns a new complete Slice.
    """
    s, surface, total_area = data
    area_from_top = backend.surface_area(surface, s.top)
    area_from_bottom = backend.surface_area(surface, s.bottom)
    slice_area = area_from_top-area_from_bottom
    s.area = slice_area
    s.weight = slice_area/total_area
    return s


def create_slices_threadpool(surface: str, interval: float) -> Iterable[Slice]:
    min_elev, max_elev = backend.min_max_elevations(surface)
    total_area = backend.surface_area(surface, max_elev)
    slices: list[Slice] = []

    elev = min_elev+interval
    while elev <= max_elev + interval:
        slices.append(Slice(top=elev, mid=elev-(interval/2), bottom=elev-interval))
        elev += interval

    WORKERS = _determine_worker_count()
    print(f"using {WORKERS} threads")
    data = zip(slices, repeat(surface, len(slices)), repeat(total_area, len(slices)))
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        complete_slices = pool.map(slice_surface_area, data, chunksize=4)

    complete_slices = sorted(complete_slices, key=lambda s: s.bottom)
    return complete_slices
