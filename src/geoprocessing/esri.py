import math
import re
from functools import cache
from pathlib import Path
from typing import Literal

import arcpy
import numpy as np
from interface import Geoprocessor

# arcpy.CheckOutExtension("3D")  # do i need this?


class EsriGeoprocessor(Geoprocessor):

    @cache
    def surface_area(self, surface: str, elevation: float) -> float:
        side = "BELOW"
        result = arcpy.SurfaceVolume_3d(surface, None, side, elevation).getMessage(1)
        # print(result)
        match = re.search(r"3D Area=\s*(?P<area>\d*[.,]?\d*)", result)
        area = match.group("area")
        return float(area)

    def min_max_elevations(self, dem: str) -> tuple[float, float]:
        r = arcpy.Raster(dem)
        return (r.minimum, r.maximum)

    # def cellsize(self, dem: str) -> tuple[float, float]:
    #     r = arcpy.Raster(dem)
        # crs: arcpy.SpatialReference = r.spatialReference
        # print(f"{crs.XYResolution = !s}")
        # print(f"{crs.abbreviation = !s}")
        # print(f"{crs.alias = !s}")
        # print(f"{crs.factoryCode = !s}")
        # print(f"{crs.name = !s}")
        # print(f"{crs.type = !s}")
        # print(f"{crs.linearUnitCode = !s}")
        # print(f"{crs.linearUnitName = !s}")
        # print(f"{crs.metersPerUnit = !s}")
        # print(f"{crs.projectionName = !s}")
        # print(f"{crs.datumName = !s}")
        # print(f"{crs.spheroidName = !s}")
        # if crs.VCS:
        #     print(f"{crs.VCS.datumName = !s}")
        #     print(f"{crs.VCS.factoryCode = !s}")
        #     print(f"{crs.VCS.linearUnitName = !s}")
        #     print(f"{crs.VCS.metersPerUnit = !s}")
        #     print(f"{crs.VCS.name = !s}")
        # else:
        #     print("no VCS")
        # return (r.meanCellWidth, r.meanCellHeight)

    # def avg_slope(self, dem: str) -> float:
    #     slope_raster = "memory\\"+Path(dem).stem+"_slope"
    #     arcpy.Slope_3d(dem, slope_raster, "DEGREE")
    #     slope = arcpy.Raster(slope_raster)
    #     arcpy.Delete_management(slope_raster)
    #     return slope.mean

    def histogram(self, dem: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        values: np.ndarray = arcpy.RasterToNumPyArray(dem, nodata_to_value=math.nan)
        values = values.flatten()
        values = values[~np.isnan(values)]
        values.sort()
        low, high = math.floor(values.min()), math.ceil(values.max())
        # print("hist:", low, high, high-low)
        bins = (x for x in range(low, high+1))
        counts, bins = np.histogram(values, bins=np.fromiter(bins, dtype=float))
        return (values, counts, bins)
