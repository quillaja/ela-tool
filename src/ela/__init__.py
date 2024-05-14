"""
The ela module provides classes and data structures for calculating the
Equilibrium Line Altitude (ELA) for glaciers based on a DEM of their surface.

`Geoprocessor` should be implemented by packages that wish to interface with
software like ArcGIS Pro or QGIS to provide the core geospatial processing.
"""
import csv
import typing

from .geoprocessing import Geoprocessor, HistData
from .methods import (AABR, AABR2D, AAR, AAR2D, ELA, AABROriginal, AAROriginal,
                      ELAMethod)
from .slice import Slice
from .slice import create_slices_threadpool as make_slices


def to_csv(results_csv: str, elas: typing.Iterable[ELA]) -> None:
    with open(results_csv, "w", encoding="utf-8") as file:
        w = csv.DictWriter(file, fieldnames=ELA._fields, lineterminator="\n")
        w.writeheader()
        for r in elas:
            w.writerow(r._asdict())
