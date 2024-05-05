import math
import re
from functools import cache

import arcpy
import numpy as np

from ela import Geoprocessor, HistData

# arcpy.CheckOutExtension("3D")  # do i need this?


class EsriGeoprocessor(Geoprocessor):
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

    # def cellsize(self, dem: str) -> tuple[float, float]:
    #     r = arcpy.Raster(dem)
        # #https://pro.arcgis.com/en/pro-app/latest/arcpy/classes/spatialreference.htm
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

    def array(self, dem: str) -> np.ndarray:
        # https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/rastertonumpyarray-function.htm
        return arcpy.RasterToNumPyArray(dem, nodata_to_value=math.nan)

    def histogram(self, dem: str) -> HistData:
        values: np.ndarray = self.array(dem)
        values = values.flatten()
        values = values[~np.isnan(values)]
        values.sort()
        low, high = math.floor(values.min()), math.ceil(values.max())
        bins = (x for x in range(low, high+1))
        counts, bins = np.histogram(values, bins=np.fromiter(bins, dtype=float))
        return HistData(values, counts, bins)

    def create_contours(self, dem: str, polylines: str, contours: list[float]) -> None:
        # https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/contour-list.htm
        arcpy.ContourList_3d(dem, polylines, contours)
