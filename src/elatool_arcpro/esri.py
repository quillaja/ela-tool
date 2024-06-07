import math
import re
from functools import cache

import arcpy
import ela
import numpy as np


class EsriGeoprocessor(ela.Geoprocessor):
    """A Geoprocessor() that interfaces with ArcGIS Pro via `arcpy`."""

    def message(self, msg: str) -> None:
        arcpy.AddMessage(msg)

    @cache
    def surface_area(self, dem: str, elevation: float) -> float:
        # https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/surface-volume.htm
        side = "BELOW"
        result = arcpy.SurfaceVolume_3d(dem, None, side, elevation).getMessage(1)
        # print(result)
        match = re.search(r"3D Area=\s*(?P<area>\d*[.,]?\d*)", result)
        area = match.group("area")
        return float(area)

    def min_max_elevations(self, dem: str) -> tuple[float, float]:
        # https://pro.arcgis.com/en/pro-app/latest/arcpy/classes/raster-object.htm
        r = arcpy.Raster(dem)
        return (r.minimum, r.maximum)

    def array(self, dem: str) -> np.ndarray:
        # https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/rastertonumpyarray-function.htm
        return arcpy.RasterToNumPyArray(dem, nodata_to_value=math.nan)

    def histogram(self, dem: str) -> ela.HistData:
        values: np.ndarray = self.array(dem)
        values = values.flatten()
        values = values[~np.isnan(values)]
        values.sort()
        low, high = math.floor(values.min()), math.ceil(values.max())
        bins = (x for x in range(low, high+1))
        counts, bins = np.histogram(values, bins=np.fromiter(bins, dtype=float))
        return ela.HistData(values, counts, bins)
