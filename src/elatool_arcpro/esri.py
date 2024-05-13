import math
import re
from functools import cache
from pathlib import Path
from typing import Iterable, Union

import arcpy
import ela
import numpy as np

# arcpy.CheckOutExtension("3D")  # do i need this?


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

    def histogram(self, dem: str) -> ela.HistData:
        values: np.ndarray = self.array(dem)
        values = values.flatten()
        values = values[~np.isnan(values)]
        values.sort()
        low, high = math.floor(values.min()), math.ceil(values.max())
        bins = (x for x in range(low, high+1))
        counts, bins = np.histogram(values, bins=np.fromiter(bins, dtype=float))
        return ela.HistData(values, counts, bins)


def elas_to_feature_class(elas: Iterable[ela.ELA], out_fc: str, crs: Union[int, arcpy.SpatialReference]) -> None:
    """
    Create or overwrite a new feature class at `out_fc` with the ELA
    information and contour lines contained within `elas`. The new feature class
    will be created with the spatial reference given in `crs`.

    :param elas: The ELA data to make into a feature class.
    :param out_fc: The feature class to create.
    :param crs: spatial reference for the new ELA feature class.
    """
    if arcpy.Exists(out_fc):
        arcpy.Delete_management(out_fc)

    ela_fields = [
        ["ela_id", "LONG"],
        ["ELA", "DOUBLE"],
        ["ratio", "DOUBLE"],
        ["interval", "DOUBLE"],
        ["method", "TEXT"],
        ["dem", "TEXT"]]
    ela_crs = crs if isinstance(
        crs, arcpy.SpatialReference) else arcpy.SpatialReference(crs)

    arcpy.CreateFeatureclass_management(
        out_path=str(Path(out_fc).parent),
        out_name=str(Path(out_fc).name),
        geometry_type="POLYLINE",
        spatial_reference=ela_crs)
    arcpy.AddFields_management(in_table=out_fc, field_description=ela_fields)

    # in order to add fields with ELA information, individual contours must be created
    # in memory, then their geometry can be written to a separate geodatabase along
    # with the field information. ContourList_3d() can write multiple contours at
    # once, but will not add a fields.
    temp_contour = r"memory\temp_contour"
    for id, r in enumerate(elas):
        arcpy.Delete_management(temp_contour)
        arcpy.ContourList_3d(r.dem, temp_contour, r.ela)
        with (arcpy.da.SearchCursor(temp_contour, ["SHAPE@"]) as source,
              arcpy.da.InsertCursor(out_fc, ["SHAPE@"]+[f[0] for f in ela_fields]) as dest):
            for source_row in source:
                dest.insertRow((source_row[0], id, r.ela, r.ratio, r.interval, r.method, r.dem))
