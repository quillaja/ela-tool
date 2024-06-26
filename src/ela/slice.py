import concurrent.futures
import multiprocessing
import os
from dataclasses import dataclass
from functools import partial

from .geoprocessing import Geoprocessor


@dataclass
class Slice:
    """
    Contains the elevation of the top, middle, and bottom of a terrain
    surface 'slice'. Area is usually the 3D surface area of the terrain within
    the slice, and the weight is the porportion of the slices area compared to
    the entire surface area.
    """
    top: float
    mid: float
    bottom: float
    area: float = 0
    weight: float = 0


def _partial_slices(min_elev: float, max_elev: float, interval: float) -> list[Slice]:
    slices: list[Slice] = []

    elev = min_elev+interval
    while elev <= max_elev + interval:
        slices.append(Slice(top=elev, mid=elev-(interval/2), bottom=elev-interval))
        elev += interval

    return slices


def _determine_worker_count() -> int:
    """Returns 2/3 the cpu count, or 4 if the system cpu count is unavailable."""
    cpus = os.cpu_count()
    if cpus:
        workers: int = round(cpus * 2/3)  # may not report correctly in docker containers
    else:
        workers = 4
    return workers


def create_slices(surface: str, interval: float, be: Geoprocessor) -> list[Slice]:
    """
    Create slices of surface with vertical height interval. Works
    in a single threaded manner.

    :param surface: a path to a DEM
    :param interval: a vertical height to make each slice
    :param be: a Geoprocessor interface to get info from `surface`
    :return: the complete slices
    """
    min_elev, max_elev = be.min_max_elevations(surface)
    total_area = be.surface_area(surface, max_elev)
    slices = _partial_slices(min_elev, max_elev, interval)

    prev_area = 0
    for s in slices:
        area = be.surface_area(surface, s.top)
        slice_area = area-prev_area
        s.area = slice_area
        s.weight = slice_area/total_area
        prev_area = area

    return slices


def _slice_surface_area(s: Slice, surface: str, total_area: float, be: Geoprocessor) -> Slice:
    """
    Calculate the surface area of the Slice, which has values for
    top, mid, and bottom, at a minimum. Returns a new complete Slice.
    """
    area_from_top = be.surface_area(surface, s.top)
    area_from_bottom = be.surface_area(surface, s.bottom)
    slice_area = area_from_top-area_from_bottom
    s.area = slice_area
    s.weight = slice_area/total_area
    return s


def create_slices_threadpool(surface: str, interval: float, be: Geoprocessor) -> list[Slice]:
    """
    Create slices of surface with vertical height interval. Works
    in a multithreaded manner using concurrent.futures.ThreadPoolExecutor.

    :param surface: a path to a DEM
    :param interval: a vertical height to make each slice
    :param be: a Geoprocessor interface to get info from `surface`
    :return: the complete slices
    """
    min_elev, max_elev = be.min_max_elevations(surface)
    total_area = be.surface_area(surface, max_elev)
    incomplete_slices = _partial_slices(min_elev, max_elev, interval)

    WORKERS = _determine_worker_count()
    work = partial(_slice_surface_area, surface=surface, total_area=total_area, be=be)
    # print(f"using {WORKERS} threads")
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        complete_slices = pool.map(work, incomplete_slices, chunksize=4)

    complete_slices = sorted(complete_slices, key=lambda s: s.bottom)
    return complete_slices


def create_slices_processpool(surface: str, interval: float, be: Geoprocessor) -> list[Slice]:
    """
    Create slices of surface with vertical height interval. Works
    in a multithreaded manner using multiprocessing.Pool.

    :param surface: a path to a DEM
    :param interval: a vertical height to make each slice
    :param be: a Geoprocessor interface to get info from `surface`
    :return: the complete slices
    """
    min_elev, max_elev = be.min_max_elevations(surface)
    total_area = be.surface_area(surface, max_elev)
    incomplete_slices = _partial_slices(min_elev, max_elev, interval)

    WORKERS = _determine_worker_count()
    work = partial(_slice_surface_area, surface=surface, total_area=total_area, be=be)
    # print(f"using {WORKERS} threads")
    with multiprocessing.Pool(processes=WORKERS) as pool:
        complete_slices = pool.map(work, incomplete_slices, chunksize=4)

    complete_slices = sorted(complete_slices, key=lambda s: s.bottom)
    return complete_slices
